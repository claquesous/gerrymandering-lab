//! Builds a single self-contained HTML report: national summary, a
//! tolerance curve for the line rule, and an SVG map plus district table for
//! every state. Pass/fail coloring, the tolerance slider and tooltips run in
//! a small inline script (`report_template.html`) over the embedded data.

use crate::fetch;
use crate::load::{self, district_code, state_abbr, state_name, Feature};
use crate::scorecard::{self, Scorecard};
use anyhow::{bail, Result};
use geo::{BoundingRect, InteriorPoint, MultiPolygon, Simplify};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

const TEMPLATE: &str = include_str!("report_template.html");
/// Longer side of each state map, in SVG pixels.
const MAP_PX: f64 = 640.0;
/// Geometry is simplified to this many pixels.
const SIMPLIFY_PX: f64 = 0.6;

struct View {
    xmin: f64,
    ymax: f64,
    scale: f64,
    w: f64,
    h: f64,
}

impl View {
    fn fit(features: &[&Feature]) -> View {
        let (mut xmin, mut ymin, mut xmax, mut ymax) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
        for f in features {
            if let Some(r) = f.geom.bounding_rect() {
                xmin = xmin.min(r.min().x);
                ymin = ymin.min(r.min().y);
                xmax = xmax.max(r.max().x);
                ymax = ymax.max(r.max().y);
            }
        }
        let scale = MAP_PX / (xmax - xmin).max(ymax - ymin);
        View { xmin, ymax, scale, w: (xmax - xmin) * scale, h: (ymax - ymin) * scale }
    }

    fn px(&self, x: f64, y: f64) -> (f64, f64) {
        ((x - self.xmin) * self.scale, (self.ymax - y) * self.scale)
    }

    fn path(&self, geom: &MultiPolygon<f64>) -> String {
        let simple: MultiPolygon<f64> = geom.simplify(SIMPLIFY_PX / self.scale);
        let mut d = String::new();
        for poly in &simple.0 {
            for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
                if ring.0.len() < 4 {
                    continue;
                }
                for (i, c) in ring.0[..ring.0.len() - 1].iter().enumerate() {
                    let (x, y) = self.px(c.x, c.y);
                    let _ = write!(d, "{}{:.1} {:.1}", if i == 0 { "M" } else { "L" }, x, y);
                }
                d.push('Z');
            }
        }
        d
    }
}

fn svg_map(view: &View, card: &Scorecard, districts: &[&Feature], counties: &[&Feature]) -> String {
    let split: HashSet<&str> = card
        .districts
        .iter()
        .filter(|r| r.geoid[..2] == districts[0].state_fp)
        .flat_map(|r| r.partial_county_geoids.split_whitespace())
        .collect();
    let mut s = String::new();
    let _ = write!(
        s,
        r#"<svg class="map" data-congress="{}" viewBox="-4 -4 {:.0} {:.0}" role="img" aria-label="{} districts, {}th Congress">"#,
        card.congress,
        view.w + 8.0,
        view.h + 8.0,
        state_name(&districts[0].state_fp),
        card.congress
    );
    s.push_str(r#"<g class="dst">"#);
    for d in districts {
        let _ = write!(s, r#"<path data-g="{}" d="{}"/>"#, d.geoid, view.path(&d.geom));
    }
    s.push_str(r#"</g><g class="cty">"#);
    for c in counties {
        let class = if split.contains(c.geoid.as_str()) { r#" class="split""# } else { "" };
        let _ = write!(s, r#"<path{class} d="{}"/>"#, view.path(&c.geom));
    }
    s.push_str(r#"</g><g class="dline">"#);
    for d in districts {
        let _ = write!(s, r#"<path d="{}"/>"#, view.path(&d.geom));
    }
    s.push_str(r#"</g><g class="wit">"#);
    for d in districts {
        // One line per crossing count; the page shows the one that matters
        // for the current setting.
        for (n, w) in card.witnesses.get(&d.geoid).into_iter().flatten().enumerate() {
            let [a, b] = w.xy;
            let (x1, y1) = view.px(a[0], a[1]);
            let (x2, y2) = view.px(b[0], b[1]);
            let _ = write!(
                s,
                r#"<g data-g="{}" data-n="{}"><line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}"/><circle cx="{x1:.1}" cy="{y1:.1}" r="3"/><circle cx="{x2:.1}" cy="{y2:.1}" r="3"/></g>"#,
                d.geoid,
                n + 1
            );
        }
    }
    s.push_str(r#"</g><g class="lbl">"#);
    for d in districts {
        let Some(r) = d.geom.bounding_rect() else { continue };
        // Label only districts big enough on screen to hold one.
        if r.width() * view.scale < 22.0 || r.height() * view.scale < 14.0 {
            continue;
        }
        let largest = d.geom.0.iter().max_by(|a, b| {
            use geo::Area;
            a.unsigned_area().total_cmp(&b.unsigned_area())
        });
        if let Some(p) = largest.and_then(|p| p.interior_point()) {
            let (x, y) = view.px(p.x(), p.y());
            let code = district_code(&d.geoid);
            let text = if code == "00" { "AL" } else { code.trim_start_matches('0') };
            let _ = write!(s, r#"<text x="{x:.1}" y="{y:.1}">{text}</text>"#);
        }
    }
    s.push_str("</g></svg>");
    s
}

fn rows_json(card: &Scorecard) -> Value {
    let rows: Vec<Value> = card
        .districts
        .iter()
        .map(|r| {
            json!({
                "g": r.geoid,
                "l": r.label,
                "c": r.county_pass,
                "w": r.whole_counties,
                "pc": r.partial_counties,
                "pl": r.partial_county_list,
                // Depth profile (km) that decides the line rule for any
                // number of allowed crossings; empty for single-district states.
                "d": parse_depths(&r.crossing_depths_km),
                // For each crossing count, the districts its worst line crosses.
                "wx": card.witnesses.get(&r.geoid).map(|ws| ws.iter().map(|w| &w.crossed).collect::<Vec<_>>()),
                "f": (r.violating_pair_fraction * 10000.0).round() / 10000.0,
            })
        })
        .collect();
    let states: Vec<Value> = card
        .states
        .iter()
        .map(|s| json!({"st": s.state, "name": s.name, "n": s.districts, "cty": s.counties, "split": s.split_counties, "extra": s.extra_county_pieces}))
        .collect();
    json!({
        "congress": card.congress,
        "min_cross_km": card.min_cross_m / 1000.0,
        "samples": card.samples,
        "resolution": card.resolution,
        "rows": rows,
        "states": states,
    })
}

/// Crossing depth profiles (km) from `counties`, if it has been run, for the
/// county baseline on the tolerance chart.
fn county_depths(out: &Path) -> Result<Option<Vec<Vec<f64>>>> {
    let path = out.join("county_line_rule.csv");
    if !path.exists() {
        return Ok(None);
    }
    let mut rdr = csv::Reader::from_path(path)?;
    let col = rdr
        .headers()?
        .iter()
        .position(|h| h == "crossing_depths_km")
        .ok_or_else(|| anyhow::anyhow!("county_line_rule.csv is from an older version; rerun `counties`"))?;
    let mut v = Vec::new();
    for row in rdr.records() {
        v.push(parse_depths(&row?[col]));
    }
    Ok(Some(v))
}

fn parse_depths(s: &str) -> Vec<f64> {
    s.split_whitespace().filter_map(|x| x.parse::<f64>().ok()).map(|x| (x * 1000.0).round() / 1000.0).collect()
}

pub fn build(data: &Path, out: &Path) -> Result<PathBuf> {
    let cards: Vec<Scorecard> = [119, 118].into_iter().filter_map(|c| scorecard::read(out, c).ok()).collect();
    if cards.is_empty() {
        bail!("no scorecards in {} (run `score` first)", out.display());
    }
    let geoms: Vec<Vec<Feature>> = cards
        .iter()
        .map(|c| load::load_districts(&fetch::district_sources(data, c.congress)?.0))
        .collect::<Result<_>>()?;
    let counties = load::load_counties(&fetch::county_source(data))?;
    let counties_by_state = load::by_state(&counties);

    // Keyed by full name so sections appear alphabetically.
    let mut states: BTreeMap<&str, String> = BTreeMap::new();
    let all_fp: BTreeMap<&str, ()> = geoms.iter().flatten().map(|f| (f.state_fp.as_str(), ())).collect();
    for fp in all_fp.keys() {
        let abbr = state_abbr(fp).unwrap_or("??");
        let first: Vec<&Feature> = geoms[0].iter().filter(|f| f.state_fp == *fp).collect();
        let view = View::fit(&first);
        let cty = counties_by_state.get(*fp).cloned().unwrap_or_default();
        let mut maps = String::new();
        for (card, g) in cards.iter().zip(&geoms) {
            let ds: Vec<&Feature> = g.iter().filter(|f| f.state_fp == *fp).collect();
            if !ds.is_empty() {
                maps.push_str(&svg_map(&view, card, &ds, &cty));
            }
        }
        let section = format!(
            r#"<section class="state" id="st-{abbr}" data-st="{abbr}"><details><summary><span class="sname">{}</span><span class="ssum"></span></summary><div class="sbody"><div class="mapwrap">{maps}</div><div class="dtwrap"><table class="dt"></table></div></div></details></section>"#,
            state_name(fp)
        );
        states.insert(state_name(fp), section);
    }

    let data_json = json!({
        "cards": cards.iter().map(rows_json).collect::<Vec<_>>(),
        "county_depths": county_depths(out)?,
    });
    let html = TEMPLATE
        .replace("__DATA__", &data_json.to_string().replace("</", "<\\/"))
        .replace("__STATES__", &states.into_values().collect::<String>());
    let path = out.join("report.html");
    std::fs::write(&path, html)?;
    Ok(path)
}

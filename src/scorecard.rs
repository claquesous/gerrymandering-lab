//! Joins the two rules into one row per district and one per state, and
//! writes them as CSV and JSON.

use crate::county_split::{self, DistrictCounties, StateCounties};
use crate::line_rule::{self, LineResult, Params};
use crate::load::{self, district_code, district_label, state_abbr, Feature};
use anyhow::{anyhow, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistrictRow {
    pub geoid: String,
    pub label: String,
    pub state: String,
    pub at_large: bool,
    pub county_pass: bool,
    pub whole_counties: usize,
    pub partial_counties: usize,
    /// "Name (share%)" for each partial county.
    pub partial_county_list: String,
    /// Space-separated GEOIDs of the partial counties.
    pub partial_county_geoids: String,
    pub line_pass: bool,
    pub max_other_crossed: usize,
    pub violating_pair_fraction: f64,
    pub pairs_tested: usize,
    /// Longest run (km) through the second-deepest other district over the
    /// segments checked; the line rule fails at or above `min_cross_m`.
    pub second_crossing_km: f64,
    /// Districts crossed by the worst confirmed segment, e.g. "MD-2 12.3km; MD-7 1.1km".
    pub witness_crossed: String,
    pub witness_from_lon: Option<f64>,
    pub witness_from_lat: Option<f64>,
    pub witness_to_lon: Option<f64>,
    pub witness_to_lat: Option<f64>,
    pub both_pass: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateRow {
    pub state: String,
    pub districts: usize,
    pub county_pass: usize,
    pub line_pass: usize,
    pub both_pass: usize,
    pub counties: usize,
    pub split_counties: usize,
    pub extra_county_pieces: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Scorecard {
    pub congress: u32,
    pub resolution: usize,
    pub samples: usize,
    pub min_cross_m: f64,
    pub districts: Vec<DistrictRow>,
    pub states: Vec<StateRow>,
    /// Projected witness segments keyed by district GEOID, for the report.
    pub witness_xy: HashMap<String, [[f64; 2]; 2]>,
}

pub fn build(congress: u32, districts: &[Feature], parts: &[load::CountyPart], params: &Params) -> Result<Scorecard> {
    let counties = county_split::score(parts);

    let states = load::by_state(districts);
    let mut lines: HashMap<String, LineResult> = HashMap::new();
    let mut state_list: Vec<(&String, &Vec<&Feature>)> = states.iter().collect();
    // Big states first so the thread pool stays busy at the end.
    state_list.sort_by_key(|(_, ds)| std::cmp::Reverse(ds.len()));
    let results: Vec<Vec<(String, LineResult)>> = state_list
        .par_iter()
        .filter(|(_, ds)| ds.len() > 1)
        .map(|(st, ds)| {
            let t = std::time::Instant::now();
            let r = line_rule::score_state(ds, params);
            eprintln!("  {} ({} districts) line rule in {:.1}s", state_abbr(st).unwrap_or(st), ds.len(), t.elapsed().as_secs_f64());
            r
        })
        .collect();
    for r in results {
        lines.extend(r);
    }

    let mut rows = Vec::new();
    let mut witness_xy = HashMap::new();
    for d in districts {
        let at_large = district_code(&d.geoid) == "00";
        let c: &DistrictCounties = counties
            .districts
            .get(&d.geoid)
            .ok_or_else(|| anyhow!("{} is missing from the county relationship file", d.geoid))?;
        let l = lines.get(&d.geoid);
        let line_pass = l.map_or(true, |l| l.pass);
        let w = l.and_then(|l| l.witness.as_ref());
        if let Some(w) = w {
            witness_xy.insert(d.geoid.clone(), [w.from_xy, w.to_xy]);
        }
        rows.push(DistrictRow {
            geoid: d.geoid.clone(),
            label: district_label(&d.geoid),
            state: state_abbr(&d.state_fp).unwrap_or("??").to_string(),
            at_large,
            county_pass: c.pass,
            whole_counties: c.whole_counties,
            partial_counties: c.partial.len(),
            partial_county_list: c
                .partial
                .iter()
                .map(|p| format!("{} ({:.0}%)", p.name, p.share * 100.0))
                .collect::<Vec<_>>()
                .join("; "),
            partial_county_geoids: c.partial.iter().map(|p| p.geoid.as_str()).collect::<Vec<_>>().join(" "),
            line_pass,
            max_other_crossed: l.map_or(0, |l| l.max_other_crossed),
            violating_pair_fraction: l.map_or(0.0, |l| l.violating_pair_fraction),
            pairs_tested: l.map_or(0, |l| l.pairs_tested),
            second_crossing_km: l.map_or(0.0, |l| l.second_crossing_m / 1000.0),
            witness_crossed: w.map_or(String::new(), |w| {
                w.crossed.iter().map(|c| format!("{} {:.1}km", district_label(&c.geoid), c.length_m / 1000.0)).collect::<Vec<_>>().join("; ")
            }),
            witness_from_lon: w.map(|w| w.from[0]),
            witness_from_lat: w.map(|w| w.from[1]),
            witness_to_lon: w.map(|w| w.to[0]),
            witness_to_lat: w.map(|w| w.to[1]),
            both_pass: c.pass && line_pass,
        });
    }

    let mut state_rows: Vec<StateRow> = states
        .keys()
        .map(|st| {
            let ds: Vec<&DistrictRow> = rows.iter().filter(|r| &r.geoid[..2] == st).collect();
            let sc: StateCounties = counties.states.get(st).cloned().unwrap_or_default();
            StateRow {
                state: state_abbr(st).unwrap_or("??").to_string(),
                districts: ds.len(),
                county_pass: ds.iter().filter(|r| r.county_pass).count(),
                line_pass: ds.iter().filter(|r| r.line_pass).count(),
                both_pass: ds.iter().filter(|r| r.both_pass).count(),
                counties: sc.counties,
                split_counties: sc.split_counties,
                extra_county_pieces: sc.extra_pieces,
            }
        })
        .collect();
    state_rows.sort_by(|a, b| a.state.cmp(&b.state));

    Ok(Scorecard {
        congress,
        resolution: params.resolution,
        samples: params.samples,
        min_cross_m: params.min_cross_m,
        districts: rows,
        states: state_rows,
        witness_xy,
    })
}

pub fn write(card: &Scorecard, out_dir: &Path) -> Result<()> {
    std::fs::create_dir_all(out_dir)?;
    let n = card.congress;
    let mut w = csv::Writer::from_path(out_dir.join(format!("scorecard_cd{n}.csv")))?;
    for r in &card.districts {
        w.serialize(r)?;
    }
    w.flush()?;
    let mut w = csv::Writer::from_path(out_dir.join(format!("states_cd{n}.csv")))?;
    for r in &card.states {
        w.serialize(r)?;
    }
    w.flush()?;
    std::fs::write(out_dir.join(format!("scorecard_cd{n}.json")), serde_json::to_string_pretty(card)?)?;
    Ok(())
}

pub fn read(out_dir: &Path, congress: u32) -> Result<Scorecard> {
    let path = out_dir.join(format!("scorecard_cd{congress}.json"));
    let text = std::fs::read_to_string(&path).map_err(|e| anyhow!("{}: {e} (run `score --congress {congress}` first)", path.display()))?;
    Ok(serde_json::from_str(&text)?)
}

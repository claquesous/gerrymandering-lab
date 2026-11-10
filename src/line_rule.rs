//! Line rule: a straight line between any two points of a district may not
//! cross two or more other districts of the same state. Water and other
//! states are neutral.
//!
//! Checking every pair of points is impossible, so the check runs in two
//! stages:
//! 1. Search: sample points on the district's boundary (including its convex
//!    hull, where violations tend to be) and in its interior. Walk every pair
//!    through a raster of the state and note which other districts the
//!    segment may reach at least `min_cross_m` into.
//! 2. Confirm: re-check the worst candidates, plus an even spread of the
//!    rest, by clipping the segment against the real polygons. Only confirmed
//!    crossings decide pass/fail. The confirmed share of the sample scales
//!    the violating-pair fraction.
//!
//! A crossing is measured by its *depth*: how far the segment gets from the
//! other district's border while inside it. Length would be wrong: a line
//! between two points on one straight-looking border can run for kilometers
//! a few meters inside the neighbor (projected parallels are slightly
//! curved), which isn't crossing it in any useful sense.

use crate::load::Feature;
use crate::project::Albers;
use crate::raster::{Raster, DEPTH_UNIT, NONE};
use geo::{BooleanOps, BoundingRect, ConvexHull, Coord, Euclidean, Length, LineString, MultiLineString, MultiPolygon};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Raster cells along the longer side of each state's bounding box.
    pub resolution: usize,
    /// Sample points per district (pairs tested is about samples^2 / 2).
    pub samples: usize,
    /// A segment must reach at least this far into another district (away
    /// from its border) for it to count as crossing it.
    pub min_cross_m: f64,
    /// Candidate pairs per district re-checked against the real polygons.
    pub confirm: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params { resolution: 4096, samples: 500, min_cross_m: 250.0, confirm: 64 }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct Witness {
    /// Endpoints as [lon, lat].
    pub from: [f64; 2],
    pub to: [f64; 2],
    /// Endpoints in projected meters, for drawing.
    pub from_xy: [f64; 2],
    pub to_xy: [f64; 2],
    /// The other districts the segment crosses, longest crossing first.
    pub crossed: Vec<Crossing>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Crossing {
    pub geoid: String,
    /// How far into the district the segment reaches.
    pub depth_m: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct LineResult {
    pub pass: bool,
    /// The most other districts any tested segment crosses (confirmed).
    pub max_other_crossed: usize,
    /// Over the segments checked exactly, how far the worst segment reaches
    /// into the *second*-deepest other district it crosses. The district fails when this is at
    /// least `min_cross_m`, so it shows how close to the line a district is.
    pub second_crossing_m: f64,
    /// Estimated share of point pairs whose segment breaks the rule.
    pub violating_pair_fraction: f64,
    pub pairs_tested: usize,
    pub witness: Option<Witness>,
}

/// Scores every district of one state. `districts` must all be in the same
/// state. Any polygons work, so this also scores counties.
pub fn score_state(districts: &[&Feature], params: &Params) -> Vec<(String, LineResult)> {
    assert!(districts.len() < NONE as usize, "too many polygons for a u16 raster");
    let polys: Vec<&MultiPolygon<f64>> = districts.iter().map(|d| &d.geom).collect();
    let raster = Raster::new(&polys, params.resolution);
    let proj = Albers::for_state(&districts[0].state_fp);
    (0..districts.len())
        .into_par_iter()
        .map(|k| (districts[k].geoid.clone(), score_district(k, districts, &raster, &proj, params)))
        .collect()
}

fn score_district(k: usize, districts: &[&Feature], raster: &Raster, proj: &Albers, params: &Params) -> LineResult {
    let pts = sample_points(&districts[k].geom, k as u16, raster, params.samples);
    // Raster depth that a true crossing of `min_cross_m` must show. Cell
    // depths are coarse, so allow a cell of slack (this only widens the
    // candidate set; the exact check decides).
    let unit = DEPTH_UNIT as f64;
    let min_run = ((unit * (params.min_cross_m / raster.cell - 1.0)).floor() as u32).max(DEPTH_UNIT as u32);

    // Stage 1: raster search. Keep pairs that run through two or more other
    // districts, with the raster depth of the second-deepest crossing.
    let mut candidates: Vec<(u32, u32, Vec<u16>, u32)> = Vec::new();
    let mut pairs = 0usize;
    for i in 0..pts.len() {
        for j in i + 1..pts.len() {
            pairs += 1;
            let t = raster_crossings(raster, pts[i], pts[j], k as u16);
            if t.qualifying(min_run).nth(1).is_some() {
                candidates.push((i as u32, j as u32, t.qualifying(min_run).collect(), t.second()));
            }
        }
    }

    // Stage 2: exact check of the worst candidates, plus an even spread to
    // estimate how many raster candidates are real.
    // Worst by districts crossed, and worst by depth of the second crossing
    // (which sets `second_crossing_m`).
    let mut by_count: Vec<usize> = (0..candidates.len()).collect();
    by_count.sort_by_key(|&c| std::cmp::Reverse((candidates[c].2.len(), candidates[c].3)));
    let mut by_depth = by_count.clone();
    by_depth.sort_by_key(|&c| std::cmp::Reverse(candidates[c].3));
    let mut chosen: Vec<usize> = by_count.into_iter().take(params.confirm / 4).collect();
    chosen.extend(by_depth.into_iter().take(params.confirm / 4));
    let spread: Vec<usize> = if candidates.is_empty() {
        vec![]
    } else {
        let n = (params.confirm - chosen.len()).min(candidates.len());
        (0..n).map(|s| s * candidates.len() / n).collect()
    };
    let mut exact: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
    for &c in chosen.iter().chain(&spread) {
        let (i, j, ref ids, _) = candidates[c];
        exact.entry(c).or_insert_with(|| exact_crossings(pts[i as usize], pts[j as usize], ids, districts, params.min_cross_m / 2.0));
    }
    let confirmed = |c: &Vec<(usize, f64)>| c.iter().filter(|(_, len)| *len >= params.min_cross_m).count();
    let second = |c: &Vec<(usize, f64)>| c.get(1).map_or(0.0, |(_, len)| *len);
    let spread_confirmed = spread.iter().filter(|c| confirmed(&exact[c]) >= 2).count();
    let second_crossing_m = exact.values().map(second).fold(0.0, f64::max);

    // The witness is the violating segment with the deepest second crossing,
    // so it still demonstrates a failure at any tolerance below
    // `second_crossing_m`.
    let mut best: Option<(f64, usize)> = None; // (second depth, candidate)
    let mut max_count = 0;
    for (&c, crossings) in &exact {
        let n = confirmed(crossings);
        if n >= 2 {
            max_count = max_count.max(n);
            if best.map_or(true, |(bd, _)| second(crossings) > bd) {
                best = Some((second(crossings), c));
            }
        }
    }

    let raw_fraction = candidates.len() as f64 / pairs.max(1) as f64;
    let rate = if spread.is_empty() { 0.0 } else { spread_confirmed as f64 / spread.len() as f64 };
    match best {
        Some((_, c)) => {
            let (i, j) = (candidates[c].0, candidates[c].1);
            let (a, b) = (pts[i as usize], pts[j as usize]);
            let (alon, alat) = proj.inverse(a.x, a.y);
            let (blon, blat) = proj.inverse(b.x, b.y);
            LineResult {
                pass: false,
                max_other_crossed: max_count,
                second_crossing_m,
                // At least one confirmed violation, even if the spread found none.
                violating_pair_fraction: (raw_fraction * rate).max(1.0 / pairs as f64),
                pairs_tested: pairs,
                witness: Some(Witness {
                    from: [alon, alat],
                    to: [blon, blat],
                    from_xy: [a.x, a.y],
                    to_xy: [b.x, b.y],
                    crossed: exact[&c]
                        .iter()
                        .filter(|(_, len)| *len >= params.min_cross_m)
                        .map(|&(d, depth)| Crossing { geoid: districts[d].geoid.clone(), depth_m: depth })
                        .collect(),
                }),
            }
        }
        None => LineResult {
            pass: true,
            max_other_crossed: max_single_crossing(&pts, raster, k as u16, min_run),
            second_crossing_m,
            violating_pair_fraction: 0.0,
            pairs_tested: pairs,
            witness: None,
        },
    }
}

/// The deepest raster cell a segment reaches in each other polygon, for up
/// to `TALLY_MAX` distinct polygons (a segment inside one polygon rarely
/// touches more).
struct Tally {
    ids: [u16; TALLY_MAX],
    depth: [u32; TALLY_MAX],
    len: usize,
}

const TALLY_MAX: usize = 32;

impl Tally {
    /// Polygons the segment reaches at least `min_depth` raster units into.
    fn qualifying(&self, min_depth: u32) -> impl Iterator<Item = u16> + '_ {
        (0..self.len).filter(move |&i| self.depth[i] >= min_depth).map(|i| self.ids[i])
    }

    /// Raster depth of the second-deepest crossing.
    fn second(&self) -> u32 {
        let (mut top, mut second) = (0, 0);
        for &n in &self.depth[..self.len] {
            if n > top {
                second = top;
                top = n;
            } else if n > second {
                second = n;
            }
        }
        second
    }
}

fn raster_crossings(raster: &Raster, a: Coord<f64>, b: Coord<f64>, own: u16) -> Tally {
    let mut t = Tally { ids: [0; TALLY_MAX], depth: [0; TALLY_MAX], len: 0 };
    let mut last = (NONE, usize::MAX);
    raster.walk(a, b, |v, d| {
        if v == NONE || v == own {
            return;
        }
        // Consecutive cells are usually the same polygon; skip the search.
        let slot = if last.0 == v {
            last.1
        } else {
            match t.ids[..t.len].iter().position(|&id| id == v) {
                Some(i) => i,
                None if t.len < TALLY_MAX => {
                    t.ids[t.len] = v;
                    t.len += 1;
                    t.len - 1
                }
                None => return,
            }
        };
        t.depth[slot] = t.depth[slot].max(d as u32);
        last = (v, slot);
    });
    t
}

/// Exact depth (meters) segment a-b reaches into each district in `ids`,
/// deepest first. Depths below `floor` may come back as 0.
fn exact_crossings(a: Coord<f64>, b: Coord<f64>, ids: &[u16], districts: &[&Feature], floor: f64) -> Vec<(usize, f64)> {
    let seg = MultiLineString(vec![LineString(vec![a, b])]);
    let mut v: Vec<(usize, f64)> =
        ids.iter().map(|&d| (d as usize, crossing_depth(&districts[d as usize].geom, &seg, floor))).collect();
    v.sort_by(|x, y| y.1.total_cmp(&x.1));
    v
}

/// How far `seg` gets from the border of `geom` while inside it.
fn crossing_depth(geom: &MultiPolygon<f64>, seg: &MultiLineString<f64>, floor: f64) -> f64 {
    let mut best: f64 = 0.0;
    for piece in &geom.clip(seg, false).0 {
        for line in piece.lines() {
            let len = ((line.end.x - line.start.x).powi(2) + (line.end.y - line.start.y).powi(2)).sqrt();
            // Both ends of a clipped piece lie on the border, so no point of
            // it is deeper than half its length.
            if len / 2.0 <= best.max(floor) {
                continue;
            }
            // Depth changes no faster than position along the line, so this
            // spacing under-reads by at most half a step.
            let n = ((len / (floor / 2.0).max(1.0)).ceil() as usize).clamp(2, 64);
            for k in 1..n {
                let t = k as f64 / n as f64;
                let p = Coord { x: line.start.x + t * (line.end.x - line.start.x), y: line.start.y + t * (line.end.y - line.start.y) };
                best = best.max(distance_to_border(geom, p));
            }
        }
    }
    best
}

fn distance_to_border(geom: &MultiPolygon<f64>, p: Coord<f64>) -> f64 {
    let mut best = f64::INFINITY;
    for poly in &geom.0 {
        for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
            for l in ring.lines() {
                let (dx, dy) = (l.end.x - l.start.x, l.end.y - l.start.y);
                let len2 = dx * dx + dy * dy;
                let t = if len2 > 0.0 { (((p.x - l.start.x) * dx + (p.y - l.start.y) * dy) / len2).clamp(0.0, 1.0) } else { 0.0 };
                let (cx, cy) = (l.start.x + t * dx - p.x, l.start.y + t * dy - p.y);
                best = best.min(cx * cx + cy * cy);
            }
        }
    }
    best.sqrt()
}

/// For a passing district: 1 if any sampled segment crosses another
/// district, else 0. Computed on the raster only.
fn max_single_crossing(pts: &[Coord<f64>], raster: &Raster, own: u16, min_run: u32) -> usize {
    let step = (pts.len() / 60).max(1);
    for i in (0..pts.len()).step_by(step) {
        for j in (i + 1..pts.len()).step_by(step) {
            if raster_crossings(raster, pts[i], pts[j], own).qualifying(min_run).next().is_some() {
                return 1;
            }
        }
    }
    0
}

/// Up to `n` points of the district: its convex-hull vertices, points spaced
/// evenly along its boundary, and an even spread of interior raster cells.
/// Boundary points lie exactly on the true boundary, so a segment between
/// them never starts or ends inside a neighbor.
pub fn sample_points(geom: &MultiPolygon<f64>, own: u16, raster: &Raster, n: usize) -> Vec<Coord<f64>> {
    let mut pts: Vec<Coord<f64>> = Vec::new();

    let hull = geom.convex_hull();
    let hull_pts: Vec<Coord<f64>> = hull.exterior().coords().copied().collect();
    let hull_take = (n / 5).max(1);
    pts.extend(spread(&hull_pts[..hull_pts.len().saturating_sub(1)], hull_take));

    // Even spacing along all exterior rings by arc length.
    let boundary_take = n * 7 / 15;
    let rings: Vec<&LineString<f64>> = geom.0.iter().map(|p| p.exterior()).collect();
    let total: f64 = rings.iter().map(|r| Euclidean.length(*r)).sum();
    if total > 0.0 {
        let spacing = total / boundary_take as f64;
        let mut next = spacing / 2.0;
        let mut walked = 0.0;
        for r in &rings {
            for seg in r.lines() {
                let len = ((seg.end.x - seg.start.x).powi(2) + (seg.end.y - seg.start.y).powi(2)).sqrt();
                while next <= walked + len && len > 0.0 {
                    let t = (next - walked) / len;
                    pts.push(Coord { x: seg.start.x + t * (seg.end.x - seg.start.x), y: seg.start.y + t * (seg.end.y - seg.start.y) });
                    next += spacing;
                }
                walked += len;
            }
        }
    }

    // Interior cells, evenly spread in scan order. Only the polygon's own
    // bounding box is scanned.
    let mut cells: Vec<usize> = Vec::new();
    if let Some(r) = geom.bounding_rect() {
        let (lo, hi) = (raster.to_grid(r.min()), raster.to_grid(r.max()));
        let (i0, j0) = (lo.x.floor().max(0.0) as usize, lo.y.floor().max(0.0) as usize);
        let (i1, j1) = ((hi.x.ceil() as usize + 1).min(raster.w), (hi.y.ceil() as usize + 1).min(raster.h));
        for j in j0..j1 {
            for i in i0..i1 {
                if raster.cells[j * raster.w + i] == own {
                    cells.push(j * raster.w + i);
                }
            }
        }
    }
    let interior_take = n.saturating_sub(pts.len());
    for idx in spread(&cells, interior_take) {
        pts.push(raster.center(idx % raster.w, idx / raster.w));
    }
    pts
}

/// Up to `n` items taken at even intervals.
fn spread<T: Copy>(items: &[T], n: usize) -> Vec<T> {
    if items.len() <= n {
        return items.to_vec();
    }
    (0..n).map(|i| items[i * items.len() / n]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::{polygon, MultiPolygon, Polygon};

    fn feature(geoid: &str, p: Polygon<f64>) -> Feature {
        Feature { geoid: geoid.into(), state_fp: "99".into(), geom: MultiPolygon(vec![p]) }
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon<f64> {
        polygon![(x: x0, y: y0), (x: x1, y: y0), (x: x1, y: y1), (x: x0, y: y1)]
    }

    fn params() -> Params {
        Params { resolution: 400, samples: 120, min_cross_m: 250.0, confirm: 32 }
    }

    #[test]
    fn grid_of_squares_passes() {
        let s = 10_000.0;
        let fs: Vec<Feature> = (0..9)
            .map(|i| feature(&format!("99{i:02}"), rect((i % 3) as f64 * s, (i / 3) as f64 * s, (i % 3 + 1) as f64 * s, (i / 3 + 1) as f64 * s)))
            .collect();
        let refs: Vec<&Feature> = fs.iter().collect();
        for (id, r) in score_state(&refs, &params()) {
            assert!(r.pass, "{id} failed: {r:?}");
        }
    }

    #[test]
    fn u_shape_around_two_districts_fails() {
        // District 01 is a U whose arms wrap around districts 02 and 03,
        // stacked in the gap. A line across the top of the U, from arm to
        // arm, crosses only 03; a line from the bottom-left of the U to the top
        // of the right arm crosses both 02 and 03.
        let u = polygon![
            (x: 0.0, y: 0.0), (x: 30_000.0, y: 0.0), (x: 30_000.0, y: 30_000.0), (x: 20_000.0, y: 30_000.0),
            (x: 20_000.0, y: 10_000.0), (x: 10_000.0, y: 10_000.0), (x: 10_000.0, y: 30_000.0), (x: 0.0, y: 30_000.0),
        ];
        let fs = vec![
            feature("9901", u),
            feature("9902", rect(10_000.0, 10_000.0, 20_000.0, 20_000.0)),
            feature("9903", rect(10_000.0, 20_000.0, 20_000.0, 30_000.0)),
        ];
        let refs: Vec<&Feature> = fs.iter().collect();
        let results = score_state(&refs, &params());
        let u = &results.iter().find(|(id, _)| id == "9901").unwrap().1;
        assert!(!u.pass);
        assert_eq!(u.max_other_crossed, 2);
        let w = u.witness.as_ref().unwrap();
        let mut ids: Vec<&str> = w.crossed.iter().map(|c| c.geoid.as_str()).collect();
        ids.sort();
        assert_eq!(ids, ["9902", "9903"]);
        assert!(u.second_crossing_m >= 250.0);
        assert!(u.violating_pair_fraction > 0.0);
        for id in ["9902", "9903"] {
            assert!(results.iter().find(|(i, _)| i == id).unwrap().1.pass);
        }
    }

    #[test]
    fn hugging_a_border_is_not_crossing() {
        // District 01's southern border bows 20 m north in the middle, so the
        // straight line between its two bottom corners runs 15 km inside each
        // of 02 and 03, but never more than 20 m from their border.
        let top = polygon![(x: 0.0, y: 0.0), (x: 15_000.0, y: 20.0), (x: 30_000.0, y: 0.0), (x: 30_000.0, y: 10_000.0), (x: 0.0, y: 10_000.0)];
        let left = polygon![(x: 0.0, y: -10_000.0), (x: 15_000.0, y: -10_000.0), (x: 15_000.0, y: 20.0), (x: 0.0, y: 0.0)];
        let right = polygon![(x: 15_000.0, y: -10_000.0), (x: 30_000.0, y: -10_000.0), (x: 30_000.0, y: 0.0), (x: 15_000.0, y: 20.0)];
        let fs = vec![feature("9901", top), feature("9902", left), feature("9903", right)];
        let refs: Vec<&Feature> = fs.iter().collect();
        let results = score_state(&refs, &params());
        assert!(results.iter().all(|(_, r)| r.pass), "{results:?}");
    }

    #[test]
    fn two_district_state_crosses_at_most_one() {
        // An L-shaped district wrapped around a square: lines from arm to arm
        // cross the square, but with only one other district the rule can't fail.
        let l = polygon![
            (x: 0.0, y: 0.0), (x: 20_000.0, y: 0.0), (x: 20_000.0, y: 10_000.0),
            (x: 10_000.0, y: 10_000.0), (x: 10_000.0, y: 20_000.0), (x: 0.0, y: 20_000.0),
        ];
        let fs = vec![feature("9901", l), feature("9902", rect(10_000.0, 10_000.0, 20_000.0, 20_000.0))];
        let refs: Vec<&Feature> = fs.iter().collect();
        let results = score_state(&refs, &params());
        let l = &results[0].1;
        assert!(l.pass);
        assert_eq!(l.max_other_crossed, 1);
    }
}

//! Line rule: a straight line between any two points of a district may not
//! cross two or more other districts of the same state. Water and other
//! states are neutral.
//!
//! Checking every pair of points is impossible, so the check runs in two
//! stages:
//! 1. Search: sample points on the district's boundary (including its convex
//!    hull, where violations tend to be) and in its interior. Walk every pair
//!    through a raster of the state and note which other districts the
//!    segment runs through for at least `min_cross_m`.
//! 2. Confirm: re-check the worst candidates, plus an even spread of the
//!    rest, by clipping the segment against the real polygons. Only confirmed
//!    crossings decide pass/fail. The confirmed share of the sample scales
//!    the violating-pair fraction.

use crate::load::Feature;
use crate::project::Albers;
use crate::raster::{Raster, NONE};
use geo::{BooleanOps, ConvexHull, Coord, Euclidean, Length, LineString, MultiLineString, MultiPolygon};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Raster cells along the longer side of each state's bounding box.
    pub resolution: usize,
    /// Sample points per district (pairs tested is about samples^2 / 2).
    pub samples: usize,
    /// A segment must run at least this far inside another district for it
    /// to count as crossing it.
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
    pub length_m: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct LineResult {
    pub pass: bool,
    /// The most other districts any tested segment crosses (confirmed).
    pub max_other_crossed: usize,
    /// Over the segments checked exactly, the longest run through the
    /// *second*-deepest other district. The district fails when this is at
    /// least `min_cross_m`, so it shows how close to the line a district is.
    pub second_crossing_m: f64,
    /// Estimated share of point pairs whose segment breaks the rule.
    pub violating_pair_fraction: f64,
    pub pairs_tested: usize,
    pub witness: Option<Witness>,
}

/// Scores every district of one state. `districts` must all be in the same
/// state; there can be at most 64 of them (California has 52).
pub fn score_state(districts: &[&Feature], params: &Params) -> Vec<(String, LineResult)> {
    assert!(districts.len() <= 64, "district bitmask holds 64 districts");
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
    let min_run = ((params.min_cross_m / raster.cell).ceil() as u32).max(2);

    // Stage 1: raster search. Keep pairs that run through two or more other
    // districts, with the raster depth of the second-deepest crossing.
    let mut candidates: Vec<(u32, u32, u64, u32)> = Vec::new();
    let mut pairs = 0usize;
    for i in 0..pts.len() {
        for j in i + 1..pts.len() {
            pairs += 1;
            let (mask, second) = raster_crossings(raster, pts[i], pts[j], k as u16, min_run);
            if mask.count_ones() >= 2 {
                candidates.push((i as u32, j as u32, mask, second));
            }
        }
    }

    // Stage 2: exact check of the worst candidates, plus an even spread to
    // estimate how many raster candidates are real.
    // Worst by districts crossed, and worst by depth of the second crossing
    // (which sets `second_crossing_m`).
    let mut by_count: Vec<usize> = (0..candidates.len()).collect();
    by_count.sort_by_key(|&c| std::cmp::Reverse((candidates[c].2.count_ones(), candidates[c].3)));
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
        let (i, j, mask, _) = candidates[c];
        exact.entry(c).or_insert_with(|| exact_crossings(pts[i as usize], pts[j as usize], mask, districts));
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
            let (i, j, _, _) = candidates[c];
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
                        .map(|&(d, len)| Crossing { geoid: districts[d].geoid.clone(), length_m: len })
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

/// Counts the cells segment a-b spends in each other district. Returns the
/// mask of districts it spends at least `min_run` cells in, and the cell
/// count of the second-deepest crossing.
fn raster_crossings(raster: &Raster, a: Coord<f64>, b: Coord<f64>, own: u16, min_run: u32) -> (u64, u32) {
    let mut cells = [0u32; 64];
    raster.walk(a, b, |v| {
        if v != NONE && v != own {
            cells[v as usize] += 1;
        }
    });
    let (mut mask, mut top, mut second) = (0u64, 0u32, 0u32);
    for (d, &n) in cells.iter().enumerate() {
        if n >= min_run {
            mask |= 1 << d;
        }
        if n > top {
            second = top;
            top = n;
        } else if n > second {
            second = n;
        }
    }
    (mask, second)
}

/// Exact length (meters) of segment a-b inside each district in `mask`,
/// longest first.
fn exact_crossings(a: Coord<f64>, b: Coord<f64>, mask: u64, districts: &[&Feature]) -> Vec<(usize, f64)> {
    let seg = MultiLineString(vec![LineString(vec![a, b])]);
    let mut v: Vec<(usize, f64)> = (0..districts.len())
        .filter(|&d| mask & (1 << d) != 0)
        .map(|d| (d, Euclidean.length(&districts[d].geom.clip(&seg, false))))
        .collect();
    v.sort_by(|x, y| y.1.total_cmp(&x.1));
    v
}

/// For a passing district: 1 if any sampled segment crosses another
/// district, else 0. Computed on the raster only.
fn max_single_crossing(pts: &[Coord<f64>], raster: &Raster, own: u16, min_run: u32) -> usize {
    let step = (pts.len() / 60).max(1);
    for i in (0..pts.len()).step_by(step) {
        for j in (i + 1..pts.len()).step_by(step) {
            if raster_crossings(raster, pts[i], pts[j], own, min_run).0 != 0 {
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

    // Interior cells, evenly spread in scan order.
    let cells: Vec<usize> = raster.cells.iter().enumerate().filter(|(_, &v)| v == own).map(|(i, _)| i).collect();
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

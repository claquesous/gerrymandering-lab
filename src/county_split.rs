//! County rule: a district may include at most one partial county.
//!
//! A county is *split* when its land lies in two or more districts, and every
//! district holding a piece of a split county includes it as a partial
//! county. The source data is the Census district-to-county relationship
//! table, which is built from census blocks, so land areas are exact and no
//! geometric tolerance is needed. Pieces that are only water are ignored.

use crate::load::CountyPart;
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};

#[derive(Debug, Clone, Serialize)]
pub struct PartialCounty {
    pub geoid: String,
    pub name: String,
    /// Share of the county's land that lies in this district, 0..1.
    pub share: f64,
    /// Number of districts the county is split between.
    pub pieces: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct DistrictCounties {
    pub whole_counties: usize,
    pub partial: Vec<PartialCounty>,
    pub pass: bool,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct StateCounties {
    pub counties: usize,
    pub split_counties: usize,
    /// Sum over split counties of (pieces - 1).
    pub extra_pieces: usize,
}

pub struct CountyResult {
    pub districts: HashMap<String, DistrictCounties>,
    pub states: HashMap<String, StateCounties>,
}

pub fn score(parts: &[CountyPart]) -> CountyResult {
    let land: Vec<&CountyPart> = parts.iter().filter(|p| p.part_land_m2 > 0.0).collect();

    let mut pieces: HashMap<&str, usize> = HashMap::new();
    for p in &land {
        *pieces.entry(p.county.as_str()).or_default() += 1;
    }

    let mut districts: BTreeMap<String, DistrictCounties> = BTreeMap::new();
    for p in &land {
        let d = districts
            .entry(p.district.clone())
            .or_insert_with(|| DistrictCounties { whole_counties: 0, partial: vec![], pass: true });
        let n = pieces[p.county.as_str()];
        if n > 1 {
            d.partial.push(PartialCounty {
                geoid: p.county.clone(),
                name: p.county_name.clone(),
                share: p.part_land_m2 / p.county_land_m2,
                pieces: n,
            });
        } else {
            d.whole_counties += 1;
        }
    }
    for d in districts.values_mut() {
        d.partial.sort_by(|a, b| a.geoid.cmp(&b.geoid));
        d.pass = d.partial.len() <= 1;
    }

    let mut states: HashMap<String, StateCounties> = HashMap::new();
    for (county, &n) in &pieces {
        let s = states.entry(county[..2].to_string()).or_default();
        s.counties += 1;
        if n > 1 {
            s.split_counties += 1;
            s.extra_pieces += n - 1;
        }
    }

    CountyResult { districts: districts.into_iter().collect(), states }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(d: &str, c: &str, county_land: f64, part_land: f64) -> CountyPart {
        CountyPart {
            district: d.into(),
            county: c.into(),
            county_name: format!("County {c}"),
            county_land_m2: county_land,
            part_land_m2: part_land,
        }
    }

    #[test]
    fn counts_partial_counties() {
        // County A is whole in district 01. County B is split between 01 and
        // 02. County C is split between 02 and 03. County D has a water-only
        // piece in 03, which doesn't count.
        let parts = vec![
            part("9901", "99001", 100.0, 100.0),
            part("9901", "99002", 100.0, 40.0),
            part("9902", "99002", 100.0, 60.0),
            part("9902", "99003", 100.0, 50.0),
            part("9903", "99003", 100.0, 50.0),
            part("9903", "99004", 100.0, 100.0),
            part("9901", "99004", 100.0, 0.0),
        ];
        let r = score(&parts);
        let d1 = &r.districts["9901"];
        assert_eq!((d1.whole_counties, d1.partial.len(), d1.pass), (1, 1, true));
        assert!((d1.partial[0].share - 0.4).abs() < 1e-12);
        let d2 = &r.districts["9902"];
        assert_eq!((d2.whole_counties, d2.partial.len(), d2.pass), (0, 2, false));
        let d3 = &r.districts["9903"];
        assert_eq!((d3.whole_counties, d3.partial.len(), d3.pass), (1, 1, true));
        let s = &r.states["99"];
        assert_eq!((s.counties, s.split_counties, s.extra_pieces), (4, 2, 2));
    }

    #[test]
    fn district_inside_one_big_county_passes() {
        let parts = vec![
            part("9901", "99001", 300.0, 100.0),
            part("9902", "99001", 300.0, 100.0),
            part("9903", "99001", 300.0, 100.0),
        ];
        let r = score(&parts);
        assert!(r.districts.values().all(|d| d.pass && d.partial.len() == 1));
        assert_eq!(r.states["99"].extra_pieces, 2);
    }
}

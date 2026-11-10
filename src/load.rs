//! Reads Census cartographic-boundary shapefiles and relationship tables.
//! Geometries come back already projected to meters (see `project`).

use crate::project::Albers;
use anyhow::{anyhow, Context, Result};
use geo::{Coord, LineString, MultiPolygon, Polygon};
use shapefile::dbase::{FieldValue, Record};
use shapefile::{PolygonRing, Shape};
use std::collections::HashMap;
use std::path::Path;

/// One polygon feature (a district or a county).
#[derive(Debug, Clone)]
pub struct Feature {
    /// "0650" for CA-50, "06037" for Los Angeles County.
    pub geoid: String,
    pub state_fp: String,
    pub geom: MultiPolygon<f64>,
}

/// Two-letter postal code for a state FIPS code, or `None` for territories
/// and DC, which have no voting House members.
pub fn state_abbr(fp: &str) -> Option<&'static str> {
    Some(match fp {
        "01" => "AL", "02" => "AK", "04" => "AZ", "05" => "AR", "06" => "CA", "08" => "CO",
        "09" => "CT", "10" => "DE", "12" => "FL", "13" => "GA", "15" => "HI", "16" => "ID",
        "17" => "IL", "18" => "IN", "19" => "IA", "20" => "KS", "21" => "KY", "22" => "LA",
        "23" => "ME", "24" => "MD", "25" => "MA", "26" => "MI", "27" => "MN", "28" => "MS",
        "29" => "MO", "30" => "MT", "31" => "NE", "32" => "NV", "33" => "NH", "34" => "NJ",
        "35" => "NM", "36" => "NY", "37" => "NC", "38" => "ND", "39" => "OH", "40" => "OK",
        "41" => "OR", "42" => "PA", "44" => "RI", "45" => "SC", "46" => "SD", "47" => "TN",
        "48" => "TX", "49" => "UT", "50" => "VT", "51" => "VA", "53" => "WA", "54" => "WV",
        "55" => "WI", "56" => "WY",
        _ => return None,
    })
}

/// Full state name for a FIPS code.
pub fn state_name(fp: &str) -> &'static str {
    match fp {
        "01" => "Alabama", "02" => "Alaska", "04" => "Arizona", "05" => "Arkansas", "06" => "California",
        "08" => "Colorado", "09" => "Connecticut", "10" => "Delaware", "12" => "Florida", "13" => "Georgia",
        "15" => "Hawaii", "16" => "Idaho", "17" => "Illinois", "18" => "Indiana", "19" => "Iowa",
        "20" => "Kansas", "21" => "Kentucky", "22" => "Louisiana", "23" => "Maine", "24" => "Maryland",
        "25" => "Massachusetts", "26" => "Michigan", "27" => "Minnesota", "28" => "Mississippi",
        "29" => "Missouri", "30" => "Montana", "31" => "Nebraska", "32" => "Nevada", "33" => "New Hampshire",
        "34" => "New Jersey", "35" => "New Mexico", "36" => "New York", "37" => "North Carolina",
        "38" => "North Dakota", "39" => "Ohio", "40" => "Oklahoma", "41" => "Oregon", "42" => "Pennsylvania",
        "44" => "Rhode Island", "45" => "South Carolina", "46" => "South Dakota", "47" => "Tennessee",
        "48" => "Texas", "49" => "Utah", "50" => "Vermont", "51" => "Virginia", "53" => "Washington",
        "54" => "West Virginia", "55" => "Wisconsin", "56" => "Wyoming",
        _ => "Unknown",
    }
}

/// District code within its state: "00" for at-large, otherwise "01".."52".
pub fn district_code(geoid: &str) -> &str {
    &geoid[2..]
}

/// Short label such as "CA-50" or "WY-AL".
pub fn district_label(geoid: &str) -> String {
    let st = state_abbr(&geoid[..2]).unwrap_or("??");
    match district_code(geoid) {
        "00" => format!("{st}-AL"),
        c => format!("{st}-{}", c.trim_start_matches('0')),
    }
}

fn text(rec: &Record, field: &str) -> Result<String> {
    match rec.get(field) {
        Some(FieldValue::Character(Some(s))) => Ok(s.trim().to_string()),
        Some(FieldValue::Character(None)) => Ok(String::new()),
        other => Err(anyhow!("field {field}: unexpected value {other:?}")),
    }
}

fn project_polygon(shape: &shapefile::Polygon, proj: &Albers) -> MultiPolygon<f64> {
    let ring = |pts: &[shapefile::Point]| -> LineString<f64> {
        pts.iter()
            .map(|p| {
                let (x, y) = proj.forward(p.x, p.y);
                Coord { x, y }
            })
            .collect()
    };
    // Shapefiles list each outer ring followed by its holes.
    let mut polys: Vec<Polygon<f64>> = Vec::new();
    for r in shape.rings() {
        match r {
            PolygonRing::Outer(pts) => polys.push(Polygon::new(ring(pts), vec![])),
            PolygonRing::Inner(pts) => match polys.last_mut() {
                Some(p) => p.interiors_push(ring(pts)),
                None => polys.push(Polygon::new(ring(pts), vec![])),
            },
        }
    }
    MultiPolygon(polys)
}

/// Reads every polygon in `shp`, keeping only the 50 states.
fn read_features(shp: &Path) -> Result<Vec<Feature>> {
    let mut reader = shapefile::Reader::from_path(shp).with_context(|| format!("opening {}", shp.display()))?;
    let mut out = Vec::new();
    for item in reader.iter_shapes_and_records() {
        let (shape, rec) = item?;
        let state_fp = text(&rec, "STATEFP")?;
        if state_abbr(&state_fp).is_none() {
            continue;
        }
        let poly = match shape {
            Shape::Polygon(p) => p,
            Shape::NullShape => continue,
            other => return Err(anyhow!("{}: unexpected shape type {}", shp.display(), other.shapetype())),
        };
        out.push(Feature {
            geoid: text(&rec, "GEOID")?,
            geom: project_polygon(&poly, &Albers::for_state(&state_fp)),
            state_fp,
        });
    }
    Ok(out)
}

/// Finds the single `.shp` file inside a fetched dataset directory.
pub fn find_shp(dir: &Path) -> Result<std::path::PathBuf> {
    std::fs::read_dir(dir)
        .with_context(|| format!("reading {} (run `fetch` first)", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .find(|p| p.extension().is_some_and(|x| x == "shp"))
        .ok_or_else(|| anyhow!("no .shp in {}", dir.display()))
}

/// Congressional districts, excluding the "ZZ" water/undefined pseudo-districts.
pub fn load_districts(dir: &Path) -> Result<Vec<Feature>> {
    let mut v = read_features(&find_shp(dir)?)?;
    v.retain(|f| district_code(&f.geoid) != "ZZ");
    v.sort_by(|a, b| a.geoid.cmp(&b.geoid));
    Ok(v)
}

pub fn load_counties(dir: &Path) -> Result<Vec<Feature>> {
    read_features(&find_shp(dir)?)
}

/// One row of a Census district-to-county relationship table: the piece of
/// `county` that lies in `district`.
#[derive(Debug, Clone)]
pub struct CountyPart {
    pub district: String,
    pub county: String,
    pub county_name: String,
    pub county_land_m2: f64,
    pub part_land_m2: f64,
}

/// Reads a pipe-delimited `tab20_cdXXX20_county20_natl.txt` file.
pub fn load_county_parts(path: &Path) -> Result<Vec<CountyPart>> {
    let raw = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let raw = raw.trim_start_matches('\u{feff}');
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'|').from_reader(raw.as_bytes());
    let headers = rdr.headers()?.clone();
    let col = |prefix: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| h.starts_with(prefix))
            .ok_or_else(|| anyhow!("{}: no column starting with {prefix}", path.display()))
    };
    let (cd, county, cname, cland, pland) = (
        col("GEOID_CD")?,
        col("GEOID_COUNTY")?,
        col("NAMELSAD_COUNTY")?,
        col("AREALAND_COUNTY")?,
        col("AREALAND_PART")?,
    );
    let num = |s: &str| s.parse::<f64>().unwrap_or(0.0);
    let mut out = Vec::new();
    for row in rdr.records() {
        let row = row?;
        let district = row[cd].to_string();
        if district.len() != 4 || state_abbr(&district[..2]).is_none() || district.ends_with("ZZ") {
            continue;
        }
        out.push(CountyPart {
            district,
            county: row[county].to_string(),
            county_name: row[cname].to_string(),
            county_land_m2: num(&row[cland]),
            part_land_m2: num(&row[pland]),
        });
    }
    Ok(out)
}

/// Groups features by state FIPS code.
pub fn by_state(features: &[Feature]) -> HashMap<String, Vec<&Feature>> {
    let mut m: HashMap<String, Vec<&Feature>> = HashMap::new();
    for f in features {
        m.entry(f.state_fp.clone()).or_default().push(f);
    }
    m
}

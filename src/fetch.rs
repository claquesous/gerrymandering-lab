//! Downloads Census boundary and relationship files. Uses the system `curl`
//! and `unzip` so the crate stays free of TLS/HTTP dependencies.

use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

const GENZ: &str = "https://www2.census.gov/geo/tiger";
const REL: &str = "https://www2.census.gov/geo/docs/maps-data/data/rel2020/cd-sld";

pub struct Dataset {
    pub name: &'static str,
    pub url: String,
    pub zipped: bool,
}

pub fn datasets() -> Vec<Dataset> {
    let shp = |year: u32, name: &'static str| Dataset {
        name,
        url: format!("{GENZ}/GENZ{year}/shp/{name}.zip"),
        zipped: true,
    };
    let rel = |name: &'static str| Dataset { name, url: format!("{REL}/{name}.txt"), zipped: false };
    vec![
        shp(2024, "cb_2024_us_cd119_500k"),
        shp(2023, "cb_2023_us_cd118_500k"),
        shp(2020, "cb_2020_us_county_500k"),
        rel("tab20_cd11920_county20_natl"),
        rel("tab20_cd11820_county20_natl"),
    ]
}

/// Where a dataset lives once fetched: a directory for shapefiles, a file
/// for relationship tables.
pub fn local_path(raw_dir: &Path, ds: &Dataset) -> PathBuf {
    if ds.zipped {
        raw_dir.join(ds.name)
    } else {
        raw_dir.join(format!("{}.txt", ds.name))
    }
}

pub fn fetch_all(raw_dir: &Path, force: bool) -> Result<()> {
    std::fs::create_dir_all(raw_dir)?;
    for ds in datasets() {
        let dest = local_path(raw_dir, &ds);
        if dest.exists() && !force {
            println!("have   {}", ds.name);
            continue;
        }
        println!("fetch  {}", ds.url);
        let download = if ds.zipped { raw_dir.join(format!("{}.zip", ds.name)) } else { dest.clone() };
        run(Command::new("curl").args(["-fsSL", "--retry", "3", "-o"]).arg(&download).arg(&ds.url))
            .with_context(|| format!("downloading {}", ds.url))?;
        if ds.zipped {
            std::fs::create_dir_all(&dest)?;
            run(Command::new("unzip").args(["-o", "-q"]).arg(&download).arg("-d").arg(&dest))
                .with_context(|| format!("unzipping {}", download.display()))?;
            std::fs::remove_file(&download)?;
        }
    }
    Ok(())
}

/// The district shapefile directory and county relationship table for a
/// Congress.
pub fn district_sources(data: &Path, congress: u32) -> Result<(PathBuf, PathBuf)> {
    let (shp, rel) = match congress {
        119 => ("cb_2024_us_cd119_500k", "tab20_cd11920_county20_natl.txt"),
        118 => ("cb_2023_us_cd118_500k", "tab20_cd11820_county20_natl.txt"),
        _ => bail!("only the 118th and 119th Congresses are supported"),
    };
    Ok((data.join(shp), data.join(rel)))
}

/// County shapes used to draw the report's maps (2020 counties, matching the
/// relationship tables).
pub fn county_source(data: &Path) -> PathBuf {
    data.join("cb_2020_us_county_500k")
}

fn run(cmd: &mut Command) -> Result<()> {
    let status = cmd.status()?;
    if !status.success() {
        bail!("{:?} exited with {status}", cmd);
    }
    Ok(())
}

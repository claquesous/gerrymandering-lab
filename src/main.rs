mod county_split;
mod fetch;
mod line_rule;
mod load;
mod project;
mod raster;
mod report;
mod scorecard;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

/// Scores US congressional districts against two anti-gerrymandering rules:
/// (1) at most one partial county per district, and (2) no straight line
/// inside a district crosses two or more other districts.
#[derive(Parser)]
struct Cli {
    /// Where downloaded Census files live.
    #[arg(long, default_value = "data/raw", global = true)]
    data: PathBuf,
    /// Where scorecards and the report are written.
    #[arg(long, default_value = "results", global = true)]
    out: PathBuf,
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Download district, county and relationship files from the Census Bureau.
    Fetch {
        /// Download again even if the files are already present.
        #[arg(long)]
        force: bool,
    },
    /// Score one Congress's districts and write results/scorecard_cdNNN.*
    Score {
        /// 118 or 119.
        #[arg(long, default_value_t = 119)]
        congress: u32,
        /// Raster cells along the longer side of each state.
        #[arg(long, default_value_t = 4096)]
        resolution: usize,
        /// Sample points per district for the line rule.
        #[arg(long, default_value_t = 500)]
        samples: usize,
        /// Minimum distance (meters) a line must run inside another district
        /// to count as crossing it.
        #[arg(long, default_value_t = 250.0)]
        min_cross_m: f64,
        /// Only score these states (postal codes, comma-separated).
        #[arg(long, value_delimiter = ',')]
        states: Vec<String>,
    },
    /// Apply the line rule to counties instead of districts, as a baseline.
    /// Writes results/county_line_rule.csv.
    Counties {
        #[arg(long, default_value_t = 4096)]
        resolution: usize,
        #[arg(long, default_value_t = 500)]
        samples: usize,
        #[arg(long, default_value_t = 250.0)]
        min_cross_m: f64,
    },
    /// Build results/report.html from the scorecards present.
    Report,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Fetch { force } => fetch::fetch_all(&cli.data, force),
        Cmd::Score { congress, resolution, samples, min_cross_m, states } => {
            let (shp, rel) = fetch::district_sources(&cli.data, congress)?;
            let t = std::time::Instant::now();
            let mut districts = load::load_districts(&shp)?;
            let mut parts = load::load_county_parts(&rel)?;
            if !states.is_empty() {
                let keep = |fp: &str| load::state_abbr(fp).is_some_and(|a| states.iter().any(|s| s.eq_ignore_ascii_case(a)));
                districts.retain(|d| keep(&d.state_fp));
                parts.retain(|p| keep(&p.district[..2]));
            }
            eprintln!("loaded {} districts in {:.1}s", districts.len(), t.elapsed().as_secs_f64());
            let params = line_rule::Params { resolution, samples, min_cross_m, ..Default::default() };
            let card = scorecard::build(congress, &districts, &parts, &params)?;
            scorecard::write(&card, &cli.out)?;
            let n = card.districts.len();
            let count = |f: fn(&scorecard::DistrictRow) -> bool| card.districts.iter().filter(|r| f(r)).count();
            println!(
                "{}th Congress: {n} districts | county rule pass {} | line rule pass {} | both {}  ({:.1}s)",
                congress,
                count(|r| r.county_pass),
                count(|r| r.line_pass),
                count(|r| r.both_pass),
                t.elapsed().as_secs_f64()
            );
            Ok(())
        }
        Cmd::Counties { resolution, samples, min_cross_m } => {
            let params = line_rule::Params { resolution, samples, min_cross_m, ..Default::default() };
            let path = scorecard::county_line_rule(&cli.data, &cli.out, &params)?;
            println!("wrote {}", path.display());
            Ok(())
        }
        Cmd::Report => {
            let path = report::build(&cli.data, &cli.out)?;
            println!("wrote {}", path.display());
            Ok(())
        }
    }
}

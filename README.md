# District Rule Scorecard

Scores existing US congressional districts against two anti-gerrymandering rules:

1. **County rule:** a district may include at most one *partial* county, meaning a county whose land is split with another district.
2. **Line rule:** a straight line between any two points in a district may not cross two or more *other* districts in the same state. Water and other states are neutral.

It covers the 118th and 119th Congresses and doesn't draw new maps.

## Usage

Needs Rust 1.85+ (pinned to 1.99.1 in `.tool-versions`), plus `curl` and `unzip`.

```sh
cargo run --release -- fetch                  # ~40 MB of Census files into data/raw/
cargo run --release -- score --congress 119   # results/scorecard_cd119.{csv,json}, states_cd119.csv
cargo run --release -- score --congress 118
cargo run --release -- counties               # line rule applied to counties (baseline), ~20 s
cargo run --release -- report                 # results/report.html (self-contained)
```

`score` options:

| Option | Default | Meaning |
|---|---|---|
| `--min-cross-m` | 250 | How far a line must reach into another district, measured from its border, to count as crossing it |
| `--samples` | 500 | Sample points per district (about 125k segments each) |
| `--resolution` | 4096 | Raster cells along each state's longer side |
| `--states MD,IL` | all | Score only these states |

A full national run takes about 8 seconds.

## Data

| Use | Source |
|---|---|
| District shapes | Census cartographic boundary files `cb_2024_us_cd119_500k`, `cb_2023_us_cd118_500k`. These are clipped to the shoreline, so water is a gap. |
| County splits | Census relationship tables `tab20_cd119(118)20_county20_natl.txt`, which give the exact land area of every county piece in every district |
| Map counties | `cb_2020_us_county_500k` (2020 counties, matching the relationship tables) |

The 119th Congress files reflect maps in force in January 2025. Maps redrawn after that aren't included.

## Method

- **County rule:** exact. It uses the relationship table's land areas, which are built from census blocks, so no geometric tolerance is involved. Pieces that are only water are ignored.
- **Line rule** (`src/line_rule.rs`), all in an Albers equal-area projection:
  1. Sample points on each district's convex hull, evenly along its true boundary, and spread through its interior.
  2. Walk every pair of points through a raster of the state, noting the deepest cell reached in each other district. Each cell's depth (distance to the nearest border) comes from a chamfer distance transform.
  3. Re-check the worst candidates, plus an even spread of the rest, by clipping the segment against the real polygons.
  4. Only confirmed crossings fail a district.
- **Crossings are measured by depth, not length.** Depth is how far the line gets from the other district's border. A line between two points on one border, such as a county line that follows a parallel, can run for kilometers a few meters inside the neighbor, because parallels are slightly curved in the projection. That shouldn't count as crossing.
- **County baseline:** `counties` applies the same rule to every county, with the other counties of its state in place of other districts. At 250 m, 58% of counties fail (mostly by under 2 km, along rivers and coastlines), against 89% of districts. At 5 km the figures are 4% and 41%.
- `second_crossing_km` is the depth of the worst line's *second* crossing. A district fails at any tolerance below it, which is how the report's tolerance slider works without rerunning anything.
- `violating_pair_fraction` is the raster share of failing pairs, scaled by the share of that sample confirmed exactly.

## Layout

```
src/fetch.rs          downloads (curl/unzip)
src/load.rs           shapefile + relationship-table readers
src/project.rs        Albers equal-area projection
src/county_split.rs   county rule
src/raster.rs         per-state raster + segment walker
src/line_rule.rs      line rule
src/scorecard.rs      CSV/JSON output
src/report.rs         HTML report (template: src/report_template.html)
```

`cargo test` runs unit tests on synthetic geometry: a U-shaped district wrapped around two others fails, a grid of squares passes, a line hugging a slightly bowed border doesn't count as crossing, and projection area and round-trip checks pass.

//! A per-state grid of district indices, and a walker that lists the cells a
//! straight segment passes through. This is the fast search half of the
//! line rule; `line_rule` confirms its findings against the real polygons.

use geo::{Coord, MultiPolygon};

/// Cell value for water, another state, or anything outside every district.
pub const NONE: u16 = u16::MAX;

pub struct Raster {
    /// World coordinates of the grid's lower-left corner.
    pub x0: f64,
    pub y0: f64,
    /// Cell edge length in meters.
    pub cell: f64,
    pub w: usize,
    pub h: usize,
    pub cells: Vec<u16>,
    /// Approximate distance from each cell to the nearest cell of a
    /// different value, in `DEPTH_UNIT`s per cell (3-4 chamfer distance).
    pub depth: Vec<u16>,
}

/// Depth units per cell width.
pub const DEPTH_UNIT: u16 = 3;

impl Raster {
    /// Rasterizes `polys` (index i gets value i) so the longer side of their
    /// combined bounding box spans about `resolution` cells. A cell takes
    /// the value of the polygon containing its center.
    pub fn new(polys: &[&MultiPolygon<f64>], resolution: usize) -> Self {
        let (mut xmin, mut ymin, mut xmax, mut ymax) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
        for mp in polys {
            for p in &mp.0 {
                for c in p.exterior().coords() {
                    xmin = xmin.min(c.x);
                    ymin = ymin.min(c.y);
                    xmax = xmax.max(c.x);
                    ymax = ymax.max(c.y);
                }
            }
        }
        let cell = (xmax - xmin).max(ymax - ymin) / resolution as f64;
        let w = ((xmax - xmin) / cell).ceil() as usize + 1;
        let h = ((ymax - ymin) / cell).ceil() as usize + 1;
        let mut r = Raster { x0: xmin, y0: ymin, cell, w, h, cells: vec![NONE; w * h], depth: vec![] };
        for (i, mp) in polys.iter().enumerate() {
            r.fill(mp, i as u16);
        }
        r.depth = r.distance_transform();
        r
    }

    /// Two-pass 3-4 chamfer distance from every cell to the nearest border
    /// between different values. Border cells get one cell (3 units).
    fn distance_transform(&self) -> Vec<u16> {
        let (w, h) = (self.w, self.h);
        let mut d = vec![u16::MAX; w * h];
        for j in 0..h {
            for i in 0..w {
                let v = self.cells[j * w + i];
                let border = i == 0
                    || j == 0
                    || i + 1 == w
                    || j + 1 == h
                    || self.cells[j * w + i - 1] != v
                    || self.cells[j * w + i + 1] != v
                    || self.cells[(j - 1) * w + i] != v
                    || self.cells[(j + 1) * w + i] != v;
                if border {
                    d[j * w + i] = DEPTH_UNIT;
                }
            }
        }
        let relax = |d: &mut [u16], at: usize, from: usize, cost: u16| {
            let c = d[from].saturating_add(cost);
            if c < d[at] {
                d[at] = c;
            }
        };
        for j in 1..h {
            for i in 1..w - 1 {
                let at = j * w + i;
                relax(&mut d, at, at - 1, 3);
                relax(&mut d, at, at - w, 3);
                relax(&mut d, at, at - w - 1, 4);
                relax(&mut d, at, at - w + 1, 4);
            }
        }
        for j in (0..h - 1).rev() {
            for i in (1..w - 1).rev() {
                let at = j * w + i;
                relax(&mut d, at, at + 1, 3);
                relax(&mut d, at, at + w, 3);
                relax(&mut d, at, at + w + 1, 4);
                relax(&mut d, at, at + w - 1, 4);
            }
        }
        d
    }

    /// Scanline fill with the even-odd rule over every ring of `mp`, so holes
    /// come out right without tracking which ring is which.
    fn fill(&mut self, mp: &MultiPolygon<f64>, value: u16) {
        let mut rows: Vec<Vec<f64>> = vec![Vec::new(); self.h];
        for p in &mp.0 {
            for ring in std::iter::once(p.exterior()).chain(p.interiors()) {
                for seg in ring.lines() {
                    let a = self.to_grid(seg.start);
                    let b = self.to_grid(seg.end);
                    if a.y == b.y {
                        continue;
                    }
                    let (lo, hi) = if a.y < b.y { (a, b) } else { (b, a) };
                    // Rows whose center line y = j + 0.5 lies in [lo.y, hi.y).
                    let j0 = (lo.y - 0.5).ceil().max(0.0) as usize;
                    let j1 = ((hi.y - 0.5).ceil() as usize).min(self.h);
                    for j in j0..j1 {
                        let yc = j as f64 + 0.5;
                        rows[j].push(lo.x + (yc - lo.y) * (hi.x - lo.x) / (hi.y - lo.y));
                    }
                }
            }
        }
        for (j, xs) in rows.iter_mut().enumerate() {
            xs.sort_by(|a, b| a.total_cmp(b));
            for pair in xs.chunks_exact(2) {
                // Cells whose center x = i + 0.5 lies in [pair[0], pair[1]).
                let i0 = (pair[0] - 0.5).ceil().max(0.0) as usize;
                let i1 = ((pair[1] - 0.5).ceil().max(0.0) as usize).min(self.w);
                for i in i0..i1 {
                    self.cells[j * self.w + i] = value;
                }
            }
        }
    }

    /// World coordinates to continuous grid coordinates.
    pub fn to_grid(&self, c: Coord<f64>) -> Coord<f64> {
        Coord { x: (c.x - self.x0) / self.cell, y: (c.y - self.y0) / self.cell }
    }

    /// World coordinates of the center of cell (i, j).
    pub fn center(&self, i: usize, j: usize) -> Coord<f64> {
        Coord { x: self.x0 + (i as f64 + 0.5) * self.cell, y: self.y0 + (j as f64 + 0.5) * self.cell }
    }

    #[cfg(test)]
    pub fn get(&self, i: i64, j: i64) -> u16 {
        self.get_with_depth(i, j).0
    }

    /// Cell value and its depth (see `depth`).
    pub fn get_with_depth(&self, i: i64, j: i64) -> (u16, u16) {
        if i < 0 || j < 0 || i as usize >= self.w || j as usize >= self.h {
            (NONE, 0)
        } else {
            let at = j as usize * self.w + i as usize;
            (self.cells[at], self.depth[at])
        }
    }

    /// Calls `f` with the value and depth of every cell the segment a-b
    /// passes through, in order (Amanatides & Woo voxel traversal).
    pub fn walk(&self, a: Coord<f64>, b: Coord<f64>, mut f: impl FnMut(u16, u16)) {
        let (a, b) = (self.to_grid(a), self.to_grid(b));
        let (mut i, mut j) = (a.x.floor() as i64, a.y.floor() as i64);
        let (ie, je) = (b.x.floor() as i64, b.y.floor() as i64);
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let (si, sj) = (dx.signum() as i64, dy.signum() as i64);
        let next = |p: f64, cell: i64, step: i64, d: f64| {
            if step == 0 {
                f64::INFINITY
            } else {
                let edge = if step > 0 { cell + 1 } else { cell } as f64;
                (edge - p) / d
            }
        };
        let mut tx = next(a.x, i, si, dx);
        let mut ty = next(a.y, j, sj, dy);
        let (ddx, ddy) = (if dx != 0.0 { 1.0 / dx.abs() } else { f64::INFINITY }, if dy != 0.0 { 1.0 / dy.abs() } else { f64::INFINITY });
        let steps = (ie - i).abs() + (je - j).abs();
        let (v, d) = self.get_with_depth(i, j);
        f(v, d);
        for _ in 0..steps {
            if tx < ty {
                i += si;
                tx += ddx;
            } else {
                j += sj;
                ty += ddy;
            }
            let (v, d) = self.get_with_depth(i, j);
            f(v, d);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::{polygon, MultiPolygon};

    fn square(x0: f64, y0: f64, s: f64) -> MultiPolygon<f64> {
        MultiPolygon(vec![polygon![(x: x0, y: y0), (x: x0 + s, y: y0), (x: x0 + s, y: y0 + s), (x: x0, y: y0 + s)]])
    }

    #[test]
    fn fills_squares_and_holes() {
        let a = square(0.0, 0.0, 10.0);
        let mut b = square(10.0, 0.0, 10.0);
        b.0[0].interiors_push(vec![(13.0, 3.0), (17.0, 3.0), (17.0, 7.0), (13.0, 7.0)]);
        let r = Raster::new(&[&a, &b], 20);
        assert_eq!(r.cell, 1.0);
        assert_eq!(r.get(2, 2), 0);
        assert_eq!(r.get(12, 2), 1);
        assert_eq!(r.get(15, 5), NONE);
        let count = |v: u16| r.cells.iter().filter(|&&c| c == v).count();
        assert_eq!(count(0), 100);
        assert_eq!(count(1), 100 - 16);
        // Depth grows toward the middle of square 0 (cells 0..10).
        let depth = |i: usize, j: usize| r.depth[j * r.w + i];
        assert_eq!(depth(0, 0), DEPTH_UNIT);
        assert_eq!(depth(9, 5), DEPTH_UNIT);
        assert_eq!(depth(4, 4), 5 * DEPTH_UNIT);
    }

    #[test]
    fn walk_visits_connected_cells() {
        let a = square(0.0, 0.0, 10.0);
        let r = Raster::new(&[&a], 10);
        let mut n = 0;
        r.walk(Coord { x: 0.5, y: 0.5 }, Coord { x: 9.5, y: 3.5 }, |_, _| n += 1);
        // 9 steps in x plus 3 in y, plus the start cell.
        assert_eq!(n, 13);
        let mut seen = Vec::new();
        r.walk(Coord { x: 9.5, y: 9.5 }, Coord { x: 9.5, y: 0.5 }, |v, _| seen.push(v));
        assert_eq!(seen, vec![0; 10]);
    }
}

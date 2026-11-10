//! Albers equal-area conic projection on the GRS80 ellipsoid (Snyder, "Map
//! Projections: A Working Manual", pp. 101-102). Census shapefiles are NAD83
//! lon/lat, which is GRS80 for our purposes. Areas and "straight lines" are
//! both measured in the projected plane, in meters.

const A: f64 = 6_378_137.0;
const E2: f64 = 0.006_694_380_022_90;

#[derive(Debug, Clone, Copy)]
pub struct Albers {
    lon0: f64,
    n: f64,
    c: f64,
    rho0: f64,
}

fn q(phi: f64) -> f64 {
    let e = E2.sqrt();
    let s = phi.sin();
    (1.0 - E2) * (s / (1.0 - E2 * s * s) - (1.0 / (2.0 * e)) * ((1.0 - e * s) / (1.0 + e * s)).ln())
}

fn m(phi: f64) -> f64 {
    let s = phi.sin();
    phi.cos() / (1.0 - E2 * s * s).sqrt()
}

impl Albers {
    /// All parameters in degrees.
    pub fn new(lat1: f64, lat2: f64, lat0: f64, lon0: f64) -> Self {
        let (p1, p2, p0) = (lat1.to_radians(), lat2.to_radians(), lat0.to_radians());
        let (m1, m2) = (m(p1), m(p2));
        let (q1, q2, q0) = (q(p1), q(p2), q(p0));
        let n = (m1 * m1 - m2 * m2) / (q2 - q1);
        let c = m1 * m1 + n * q1;
        let rho0 = A * (c - n * q0).sqrt() / n;
        Albers { lon0: lon0.to_radians(), n, c, rho0 }
    }

    /// EPSG:5070, NAD83 / Conus Albers.
    pub fn conus() -> Self {
        Self::new(29.5, 45.5, 23.0, -96.0)
    }

    /// EPSG:3338, NAD83 / Alaska Albers.
    pub fn alaska() -> Self {
        Self::new(55.0, 65.0, 50.0, -154.0)
    }

    /// ESRI:102007, Hawaii Albers Equal Area Conic.
    pub fn hawaii() -> Self {
        Self::new(8.0, 18.0, 3.0, -157.0)
    }

    /// Projection to use for a state, by two-digit FIPS code.
    pub fn for_state(state_fp: &str) -> Self {
        match state_fp {
            "02" => Self::alaska(),
            "15" => Self::hawaii(),
            _ => Self::conus(),
        }
    }

    pub fn forward(&self, lon: f64, lat: f64) -> (f64, f64) {
        let rho = A * (self.c - self.n * q(lat.to_radians())).sqrt() / self.n;
        let mut dlon = lon.to_radians() - self.lon0;
        // Keep Aleutian longitudes east of 180 on the same sheet.
        if dlon > std::f64::consts::PI {
            dlon -= 2.0 * std::f64::consts::PI;
        } else if dlon < -std::f64::consts::PI {
            dlon += 2.0 * std::f64::consts::PI;
        }
        let theta = self.n * dlon;
        (rho * theta.sin(), self.rho0 - rho * theta.cos())
    }

    pub fn inverse(&self, x: f64, y: f64) -> (f64, f64) {
        let dy = self.rho0 - y;
        let rho = (x * x + dy * dy).sqrt();
        let theta = x.atan2(dy);
        let qv = (self.c - rho * rho * self.n * self.n / (A * A)) / self.n;
        let e = E2.sqrt();
        let mut phi = (qv / 2.0).asin();
        for _ in 0..15 {
            let s = phi.sin();
            let one = 1.0 - E2 * s * s;
            let delta = one * one / (2.0 * phi.cos())
                * (qv / (1.0 - E2) - s / one + (1.0 / (2.0 * e)) * ((1.0 - e * s) / (1.0 + e * s)).ln());
            phi += delta;
            if delta.abs() < 1e-12 {
                break;
            }
        }
        ((self.lon0 + theta / self.n).to_degrees(), phi.to_degrees())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn origin_maps_to_zero() {
        let (x, y) = Albers::conus().forward(-96.0, 23.0);
        assert!(x.abs() < 1e-6 && y.abs() < 1e-6);
    }

    #[test]
    fn preserves_ellipsoidal_area() {
        // Exact area of a lon/lat quadrangle on the ellipsoid is
        // (a^2 / 2) * dlon * (q(lat2) - q(lat1)). Densify its parallels
        // (which project to arcs) and compare the shoelace area.
        let p = Albers::conus();
        let (lon1, lon2, lat1, lat2): (f64, f64, f64, f64) = (-105.0, -102.0, 37.0, 41.0);
        let exact = A * A / 2.0 * (lon2 - lon1).to_radians() * (q(lat2.to_radians()) - q(lat1.to_radians()));
        let mut ring = Vec::new();
        let steps = 2000;
        for i in 0..=steps {
            ring.push(p.forward(lon1 + (lon2 - lon1) * i as f64 / steps as f64, lat1));
        }
        for i in 0..=steps {
            ring.push(p.forward(lon2 - (lon2 - lon1) * i as f64 / steps as f64, lat2));
        }
        let mut area = 0.0;
        for i in 0..ring.len() {
            let (x1, y1) = ring[i];
            let (x2, y2) = ring[(i + 1) % ring.len()];
            area += x1 * y2 - x2 * y1;
        }
        let area = area.abs() / 2.0;
        assert!((area - exact).abs() / exact < 1e-6, "{area} vs {exact}");
    }


    #[test]
    fn round_trip() {
        for p in [Albers::conus(), Albers::alaska(), Albers::hawaii()] {
            for &(lon, lat) in &[(-120.5, 37.2), (-70.1, 44.9), (-155.5, 19.6), (-150.0, 61.2)] {
                let (x, y) = p.forward(lon, lat);
                let (lon2, lat2) = p.inverse(x, y);
                assert!((lon - lon2).abs() < 1e-9 && (lat - lat2).abs() < 1e-9);
            }
        }
    }
}

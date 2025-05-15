#![allow(dead_code)]
use mag_calib::SphericalDownSampler;
use plotly::{
    color::Rgb,
    common::{ColorScale, ColorScalePalette, Line, Marker, Mode},
    HeatMap, Plot, Scatter3D,
};
use std::f64::consts::{PI, TAU};

pub fn scatter3d(points_vec: Vec<&[[f64; 3]]>, labels: Vec<String>) -> String {
    let mut plot = Plot::new();
    for (i, points) in points_vec.iter().enumerate() {
        let mut x = Vec::with_capacity(points.len());
        let mut y = Vec::with_capacity(points.len());
        let mut z = Vec::with_capacity(points.len());

        for &[px, py, pz] in points.iter() {
            x.push(px);
            y.push(py);
            z.push(pz);
        }

        let trace = Scatter3D::new(x, y, z)
            .mode(Mode::Markers)
            .marker(Marker::new().size(3))
            .name(labels[i].clone());
        plot.add_trace(trace);
    }
    return plot.to_html();
}

pub fn heatmap_buckets(downsampler: &SphericalDownSampler, max_count: Option<usize>) -> String {
    let lon_div = downsampler.lon_div;
    let height_div = downsampler.height_div;

    // Initialize 2D array of zeros
    let mut z = vec![vec![0.0; lon_div]; height_div]; // z[lat][lon]

    // Fill heatmap values
    for (&(lon, lat), cell) in &downsampler.buckets {
        if let Some(max_count) = max_count {
            z[lat][lon] = std::cmp::min(cell.points.len(), max_count) as f64;
        } else {
            z[lat][lon] = cell.points.len() as f64;
        }
    }

    let heatmap = HeatMap::new_z(z).color_scale(ColorScale::Palette(ColorScalePalette::Jet));

    let mut plot = Plot::new();
    plot.add_trace(heatmap);
    plot.to_html()
}

fn create_sphere_graticules(
    lon_div: usize,
    height_div: usize,
    radius: f64,
) -> Vec<Box<Scatter3D<f64, f64, f64>>> {
    let mut traces: Vec<Box<Scatter3D<f64, f64, f64>>> = Vec::new();

    for lat_idx in 1..height_div {
        // Divide by Z instead of phi
        let z_frac = lat_idx as f64 / height_div as f64;
        let z_val = radius * (1.0 - 2.0 * z_frac); // from +r to -r
        let r_xy = (radius.powi(2) - z_val.powi(2)).sqrt(); // = radius * sin(phi)

        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut z = Vec::new();

        for step in 0..=lon_div {
            let theta = TAU * step as f64 / lon_div as f64;
            x.push(r_xy * theta.cos());
            y.push(r_xy * theta.sin());
            z.push(z_val);
        }

        traces.push(
            Scatter3D::new(x, y, z)
                .mode(Mode::Lines)
                .line(Line::new().color("rgb(0,0,0)").width(5.0))
                .name(format!("lat {}", lat_idx)),
        );
    }

    // Longitude lines (constant theta)
    for lon_idx in 0..lon_div {
        let theta = TAU * lon_idx as f64 / lon_div as f64;
        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut z = Vec::new();

        for step in 0..=height_div {
            let phi = PI * step as f64 / height_div as f64;
            let xp = radius * phi.sin() * theta.cos();
            let yp = radius * phi.sin() * theta.sin();
            let zp = radius * phi.cos();
            x.push(xp);
            y.push(yp);
            z.push(zp);
        }

        traces.push(
            Scatter3D::new(x, y, z)
                .mode(Mode::Lines)
                .line(Line::new().color("rgb(0,0,0)").width(5.0))
                .name(format!("lon {}", lon_idx)),
        );
    }

    traces
}

pub fn scatter3d_buckets_with_graticules(
    downsampler: &SphericalDownSampler,
    graticules_radius: f64,
) -> String {
    let mut x = Vec::with_capacity(downsampler.buckets.len());
    let mut y = Vec::with_capacity(downsampler.buckets.len());
    let mut z = Vec::with_capacity(downsampler.buckets.len());
    let mut c: Vec<Rgb> = Vec::with_capacity(downsampler.buckets.len());

    let colors: [Rgb; 4] = [
        Rgb::new(31, 119, 180), // Blue
        Rgb::new(255, 127, 14), // Orange
        Rgb::new(44, 160, 44),  // Green
        Rgb::new(214, 39, 40),  // Red
    ];

    for ((lat, lon), cell) in downsampler.buckets.iter() {
        let color_index = (lon % 2) + 2 * (lat % 2);

        for p in &cell.points {
            x.push(p[0]);
            y.push(p[1]);
            z.push(p[2]);
            c.push(colors[color_index]);
        }
    }

    let trace = Scatter3D::new(x, y, z)
        .mode(Mode::Markers)
        .marker(Marker::new().size(3).color_array(c));

    let mut plot = Plot::new();
    plot.add_trace(trace);

    // Add graticules
    let graticules = create_sphere_graticules(
        downsampler.lon_div,
        downsampler.height_div,
        graticules_radius,
    );
    for trace in graticules {
        plot.add_trace(trace);
    }

    plot.to_html()
}

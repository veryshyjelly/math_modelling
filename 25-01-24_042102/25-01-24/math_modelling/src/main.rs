use std::iter::zip;

use ndarray::prelude::*;
use plotlib::{page::Page, repr::Plot, style::LineStyle, view::ContinuousView};

fn main() {
    const K: usize = 1000;
    let t = Array1::linspace(0.0, 1.0, K);
    let (a, b, c, d) = (15.0, 0.1, 0.01, 10.0);
    let (n0, p0) = (2000.0, 100.0);
    let dt = 1.0 / (K as f64);
    let dn_by_dt = |nt: f64, pt: f64| nt * a - b * nt * pt;
    let dp_by_dt = |nt: f64, pt: f64| c * nt * pt - d * pt;

    let (n, p) = solve_dual_forward(dn_by_dt, dp_by_dt, dt, n0, p0, K);

    let data_n: Vec<(f64, f64)> = zip(t.clone(), n).collect();
    let data_p: Vec<(f64, f64)> = zip(t, p).collect();

    let s1: Plot = Plot::new(data_n).line_style(LineStyle::new().colour("#35C788"));
    let s2: Plot = Plot::new(data_p).line_style(LineStyle::new().colour("#DD3355"));

    let v = ContinuousView::new()
        .add(s1)
        .add(s2)
        .x_range(0.0, 1.5)
        .y_range(-10.0, n0)
        .x_label("Time");
    Page::single(&v).save("prey_predator.svg").unwrap();
}

fn solve_dual_forward<N, P>(
    dn_by_dt: N,
    dp_by_dt: P,
    dt: f64,
    n0: f64,
    p0: f64,
    k: usize,
) -> (Array1<f64>, Array1<f64>)
where
    N: Fn(f64, f64) -> f64,
    P: Fn(f64, f64) -> f64,
{
    let mut n = Array1::zeros(k);
    let mut p = Array1::zeros(k);
    n[0] = n0;
    p[0] = p0;

    for i in 1..k {
        n[i] = n[i - 1] + dt * dn_by_dt(n[i - 1], p[i - 1]);
        p[i] = p[i - 1] + dt * dp_by_dt(n[i], p[i - 1]);
    }

    (n, p)
}

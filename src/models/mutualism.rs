use crate::models::ChartDrawer;
use crate::ode_solvers::ODESolver2;
use ndarray::Array1;

pub fn mutualism1(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (alpha1, beta1, alpha2, beta2) = (15., 0.1, 10.0, 0.01);
    let (n0, p0, tn, n_steps) = (2000., 100., 1., 10000);
    let dn_by_dt = |_: f64, n: f64, p: f64| alpha1 * n * (1. + beta1 * p);
    let dp_by_dt = |_: f64, n: f64, p: f64| alpha2 * p * (1. + beta2 * n);

    let t = Array1::linspace(0., tn, n_steps);
    let (nt, pt) = ode_solver.solve(dn_by_dt, dp_by_dt, tn / n_steps as f64, n_steps, &t, n0, p0);

    chart_drawer(
        "plots/mutualism1.png",
        "Mutualism 1",
        vec![(t.clone(), nt, "n"), (t, pt, "p")],
    );
}

pub fn mutualism2(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (alpha1, beta1, k1, alpha2, beta2, k2) = (15., 0.1, 300., 10.0, 0.01, 400.);
    let (n0, p0, tn, n_steps) = (2000., 100., 1., 10000);
    let dn_by_dt = |_: f64, n: f64, p: f64| alpha1 * n * (1. - (n - beta1 * p) / k1);
    let dp_by_dt = |_: f64, n: f64, p: f64| alpha2 * p * (1. - (p - beta2 * n) / k2);

    let t = Array1::linspace(0., tn, n_steps);
    let (nt, pt) = ode_solver.solve(dn_by_dt, dp_by_dt, tn / n_steps as f64, n_steps, &t, n0, p0);

    chart_drawer(
        "plots/mutualism2.png",
        "Mutualism 2",
        vec![(t.clone(), nt, "n"), (t, pt, "p")],
    );
}
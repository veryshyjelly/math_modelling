use super::ChartDrawer;
use crate::ode_solvers::ODESolver2;
use ndarray::Array1;

pub fn lotka_volterra(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (alpha1, beta1, alpha2, beta2) = (15., 0.1, 10.0, 0.01);
    let (n0, p0, tn, n_steps) = (2000., 100., 1., 10000);
    let dn_by_dt = |_: f64, n: f64, p: f64| alpha1 * n - beta1 * n * p;
    let dp_by_dt = |_: f64, n: f64, p: f64| -alpha2 * p + beta2 * n * p;

    let t = Array1::linspace(0., tn, n_steps);
    let (nt, pt) = ode_solver.solve(dn_by_dt, dp_by_dt, tn / n_steps as f64, n_steps, &t, n0, p0);

    chart_drawer(
        "plots/lotka_volterra_model.png",
        "Predator Prey - Lotka Volterra",
        vec![(t.clone(), nt, "preys"), (t, pt, "predator")],
    );
}
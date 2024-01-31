use super::ChartDrawer;
use crate::ode_solvers::ODESolver1;
use ndarray::Array1;

pub fn seasonal_capacity(ode_solver: ODESolver1, chart_drawer: ChartDrawer) {
    let (alpha, beta, gamma, k) = (0.8, 0.5, 1.5, 10.);
    let (n0, tn, n_steps) = (2., 10., 1000000);
    let dn_by_dt = |t: f64, n: f64| alpha * n * (1. - n * (1. + beta * f64::cos(gamma * t)) / k);

    let t = Array1::linspace(0., tn, n_steps);
    let nt = ode_solver.solve(dn_by_dt, tn / n_steps as f64, n_steps, &t, n0);

    chart_drawer(
        "plots/seasonal_capacity_model.png",
        "Seasonal Capacity Model",
        vec![(t, nt, "population")],
    );
}
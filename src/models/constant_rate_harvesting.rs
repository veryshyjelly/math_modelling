use crate::ode_solvers::ODESolver1;
use ndarray::Array1;

pub fn constant_rate_harvesting(
    ode_solver: ODESolver1,
    chart_drawer: fn(&str, &str, Vec<(Array1<f64>, Array1<f64>)>),
) {
    let (alpha, h, k) = (1., 2., 10.);
    let (y0, tn, n_steps) = (2., 1.6, 10000);
    let dn_by_dt = |_: f64, n: f64| alpha * n * (1. - n / k) - h;

    let t = Array1::linspace(0., tn, n_steps);
    let nt = ode_solver.solve(dn_by_dt, tn / n_steps as f64, n_steps, &t, y0);

    chart_drawer(
        "plots/constant_rate_harvesting.png",
        "Constant Rate Harvesting",
        vec![(t, nt)],
    );
}
use crate::ode_solvers::ODESolver1;
use ndarray::Array1;
use std::f64::consts::E;

pub fn generalist_vs_specialist_predator(
    ode_solver: ODESolver1,
    chart_drawer: fn(&str, &str, Vec<(Array1<f64>, Array1<f64>)>),
) {
    let (r, k) = (0.5, 10.);
    let (a, b) = (1.5, 2.5);
    let (y0, tn, n_steps) = (2., 10., 10000);
    let dn_by_dt_generalist =
        |_: f64, n: f64| r * n * (1. - n / k) - b * (1. - E.powf(-(n * n) / (a * a)));
    let dn_by_dt_specialist = |_: f64, n: f64| r * n * (1. - n / k) - b * n / (a + n);

    let t = Array1::linspace(0., tn, n_steps);
    let nt_generalist = ode_solver.solve(dn_by_dt_generalist, tn / n_steps as f64, n_steps, &t, y0);
    let nt_specialist = ode_solver.solve(dn_by_dt_specialist, tn / n_steps as f64, n_steps, &t, y0);

    chart_drawer(
        "plots/generalist_vs_specialist_predator.png",
        "Generalist vs Specialist Predator",
        vec![(t.clone(), nt_generalist), (t, nt_specialist)],
    );
}
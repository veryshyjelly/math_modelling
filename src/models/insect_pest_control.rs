use super::ChartDrawer;
use crate::ode_solvers::ODESolver2;
use ndarray::Array1;

pub fn pest_control1(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (a, b, k) = (1., 2., 10.);
    let (N0, n0, tn, n_steps) = (2., 0.2, 1.6, 10000);
    let dN_by_dt = |_: f64, N: f64, n: f64| (a * N / (N + n) - b) * N - k * N * (N + n);
    let dn_by_dt = |_: f64, _: f64, n: f64| -b * n;

    let t = Array1::linspace(0., tn, n_steps);
    let (Nt, nt) = ode_solver.solve(dN_by_dt, dn_by_dt, tn / n_steps as f64, n_steps, &t, N0, n0);

    chart_drawer(
        "plots/insect_pest_control.png",
        "Insect Pest Control",
        vec![(t.clone(), Nt, "pest"), (t, nt, "insect")],
    );
}

pub fn pest_control2(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (a, b, gamma, k) = (1., 2., 3., 10.);
    let (N0, n0, tn, n_steps) = (2., 0.2, 1.6, 10000);
    let dN_by_dt = |_: f64, N: f64, n: f64| (a * N / (N + n) - b) * N - k * N * (N + n);
    let dn_by_dt = |_: f64, N: f64, n: f64| gamma * N - b * n;

    let t = Array1::linspace(0., tn, n_steps);
    let (Nt, nt) = ode_solver.solve(dN_by_dt, dn_by_dt, tn / n_steps as f64, n_steps, &t, N0, n0);

    chart_drawer(
        "plots/insect_pest_control.png",
        "Insect Pest Control",
        vec![(t.clone(), Nt, "pest"), (t, nt, "insect")],
    );
}
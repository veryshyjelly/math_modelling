use super::ChartDrawer;
use crate::ode_solvers::ODESolver1;
use ndarray::Array1;

pub fn gompertz(ode_solver: ODESolver1, chart_drawer: ChartDrawer) {
    let (alpha, k) = (0.8, 1.);
    let (n0, tn, n_steps) = (2., 2., 1000000);
    let dn_by_dt = |_: f64, n: f64| -alpha * n * f64::ln(n / k);

    let t = Array1::linspace(0., tn, n_steps);
    let nt = ode_solver.solve(dn_by_dt, tn / n_steps as f64, n_steps, &t, n0);

    chart_drawer(
        "plots/gompertz_model.png",
        "Gompertz Model",
        vec![(t, nt, "population")],
    );
}
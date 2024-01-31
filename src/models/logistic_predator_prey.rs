use crate::models::ChartDrawer;
use crate::ode_solvers::ODESolver2;
use ndarray::Array1;

pub fn logistic_predator_prey(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (alpha1, beta, alpha2, gamma, k) = (15., 0.1, 10.0, 0.01, 200.);
    let (n0, p0, tn, n_steps) = (2000., 100., 1., 10000);
    let dn_by_dt = |_: f64, n: f64, p: f64| alpha1 * n * (1. - n / k - beta * p);
    let dp_by_dt = |_: f64, n: f64, p: f64| -alpha2 * p * (1. - gamma * n);

    let t = Array1::linspace(0., tn, n_steps);
    let (nt, pt) = ode_solver.solve(dn_by_dt, dp_by_dt, tn / n_steps as f64, n_steps, &t, n0, p0);

    chart_drawer(
        "plots/logistic_predator_prey.png",
        "Logistic Predator Prey",
        vec![(t.clone(), nt, "preys"), (t, pt, "predator")],
    );
}
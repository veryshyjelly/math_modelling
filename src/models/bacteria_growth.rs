use super::ChartDrawer;
use crate::ode_solvers::ODESolver1;
use ndarray::Array1;

pub fn bacteria_growth(ode_solver: ODESolver1, chart_drawer: ChartDrawer) {
    let r = 1.5;
    let (n0, tn, n_steps) = (0.1, 1., 100000);
    let dn_by_dt = |_: f64, n: f64| r * n.powf(2. / 3.);

    let t = Array1::linspace(0., tn, n_steps);
    let nt = ode_solver.solve(dn_by_dt, tn / n_steps as f64, n_steps, &t, n0);

    chart_drawer(
        "plots/bacteria_growth.png",
        "Bacteria Growth in Petri Dish",
        vec![(t, nt, "bacteria")],
    );
}
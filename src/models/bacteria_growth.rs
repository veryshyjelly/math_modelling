use crate::ode_solvers::OdeSolver;
use ndarray::Array1;

pub fn bacteria_growth(
    ode_solver: OdeSolver,
    chart_drawer: fn((f64, f64), (f64, f64), &str, &str, Vec<(Array1<f64>, Array1<f64>)>),
) {
    let r = 1.5;
    let (y0, tn, n_steps) = (0.1, 1., 100000);
    let dn_by_dt = |_: f64, n: f64| r * n.powf(2. / 3.);

    let t = Array1::linspace(0., tn, n_steps);
    let nt = ode_solver.solve(dn_by_dt, tn / n_steps as f64, n_steps, &t, y0);

    let max_n = nt.iter().rfold(0., |x, y| f64::max(x, *y));

    chart_drawer(
        (0., tn),
        (0., max_n),
        "plots/bacteria_growth.png",
        "Bacteria Growth in Petri Dish",
        vec![(t, nt)],
    );
}
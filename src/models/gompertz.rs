use crate::ode_solvers::OdeSolver;
use ndarray::Array1;

pub fn gompertz(
    ode_solver: OdeSolver,
    chart_drawer: fn((f64, f64), (f64, f64), &str, &str, Vec<(Array1<f64>, Array1<f64>)>),
) {
    let (alpha, k) = (0.8, 1.);
    let (y0, tn, n_steps) = (2., 2., 1000000);
    let dn_by_dt = |_: f64, n: f64| -alpha * n * f64::ln(n / k);

    let t = Array1::linspace(0., tn, n_steps);
    let nt = ode_solver.solve(dn_by_dt, tn / n_steps as f64, n_steps, &t, y0);

    let max_n = nt.iter().rfold(0., |x, y| f64::max(x, *y));

    chart_drawer(
        (0., tn),
        (0., max_n),
        "plots/gompertz_model.png",
        "Gompertz Model",
        vec![(t, nt)],
    );
}
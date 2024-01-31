use crate::models::ChartDrawer;
use crate::ode_solvers::ODESolver2;
use ndarray::Array1;

pub fn rabies_pest1(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (r, beta, u, k) = (0.1, 0.1, 10., 100.);
    let (s0, i0, tn, n_steps) = (100., 10., 1., 10000);
    let ds_by_dt = |_: f64, s: f64, i: f64| r * (s + i) * (1. - s / k) - beta * s * i;
    let di_by_dt = |_: f64, s: f64, i: f64| beta * s * i - u * i;

    let t = Array1::linspace(0., tn, n_steps);
    let (st, it) = ode_solver.solve(ds_by_dt, di_by_dt, tn / n_steps as f64, n_steps, &t, s0, i0);

    chart_drawer(
        "plots/rabies_pest1.png",
        "Rabies Pest 1",
        vec![(t.clone(), st, "susceptible"), (t, it, "infective")],
    )
}

pub fn rabies_pest2(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (r, beta, u, k, c) = (0.1, 0.1, 10., 100., 10.);
    let (s0, i0, tn, n_steps) = (100., 10., 1., 10000);
    let ds_by_dt = |_: f64, s: f64, i: f64| r * (s + i) * (1. - s / k) - beta * s * i - c * s;
    let di_by_dt = |_: f64, s: f64, i: f64| beta * s * i - u * i - c * i;

    let t = Array1::linspace(0., tn, n_steps);
    let (st, it) = ode_solver.solve(ds_by_dt, di_by_dt, tn / n_steps as f64, n_steps, &t, s0, i0);

    chart_drawer(
        "plots/rabies_pest2.png",
        "Rabies Pest 2",
        vec![(t.clone(), st, "susceptible"), (t, it, "infective")],
    )
}

pub fn rabies_pest3(ode_solver: ODESolver2, chart_drawer: ChartDrawer) {
    let (r, beta, u, k, v) = (0.1, 0.1, 10., 100., 10.);
    let (s0, i0, tn, n_steps) = (100., 10., 1., 10000);
    let ds_by_dt = |_: f64, s: f64, i: f64| r * (s + i) * (1. - s / k) - beta * (1. - v) * s * i;
    let di_by_dt = |_: f64, s: f64, i: f64| beta * (1. - v) * s * i - u * i;

    let t = Array1::linspace(0., tn, n_steps);
    let (st, it) = ode_solver.solve(ds_by_dt, di_by_dt, tn / n_steps as f64, n_steps, &t, s0, i0);

    chart_drawer(
        "plots/rabies_pest3.png",
        "Rabies Pest 3",
        vec![(t.clone(), st, "susceptible"), (t, it, "infective")],
    )
}
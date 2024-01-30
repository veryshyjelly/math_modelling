use crate::ode_solvers::ODESolver2;
use ndarray::Array1;

pub fn insect_pest_control(
    ode_solver: ODESolver2,
    chart_drawer: fn(&str, &str, Vec<(Array1<f64>, Array1<f64>)>),
) {
}
use crate::chart::draw_bitmap_line_chart;
use crate::models::bacteria_growth::bacteria_growth;
use crate::models::demographic::demographic;
use crate::models::gompertz::gompertz;
use crate::ode_solvers::{OdeSolver, SolverMethod};

mod chart;
mod models;
mod ode_solvers;

fn main() {
    let solver = OdeSolver::new(SolverMethod::ForwardEuler);

    // bacteria_growth(solver, draw_bitmap_line_chart);
    // gompertz(solver, draw_bitmap_line_chart)
    demographic(solver, draw_bitmap_line_chart);
}
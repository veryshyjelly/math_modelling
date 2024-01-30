use crate::chart::draw_bitmap_line_chart;
use crate::models::generalist_vs_specialist::generalist_vs_specialist_predator;
use crate::models::{
    bacteria_growth::bacteria_growth, constant_rate_harvesting::constant_rate_harvesting,
    demographic::demographic, gompertz::gompertz, optimal_harvesting::optimal_harvesting,
    seasonal_capacity::seasonal_capacity,
};
use crate::ode_solvers::{ODESolver1, SolverMethod};

mod chart;
mod models;
mod ode_solvers;

fn main() {
    let solver = ODESolver1::new(SolverMethod::Three8th);

    // bacteria_growth(solver, draw_bitmap_line_chart);
    // gompertz(solver, draw_bitmap_line_chart);
    // demographic(solver, draw_bitmap_line_chart);
    // seasonal_capacity(solver, draw_bitmap_line_chart);
    // constant_rate_harvesting(solver, draw_bitmap_line_chart);
    // optimal_harvesting(solver, draw_bitmap_line_chart);
    generalist_vs_specialist_predator(solver, draw_bitmap_line_chart);
}
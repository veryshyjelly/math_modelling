use chart::draw_bitmap_line_chart;
use ode_solvers::{ODESolver1, ODESolver2, SolverMethod};

mod chart;
mod models;
mod ode_solvers;

fn main() {
    // let solver = ODESolver1::new(SolverMethod::Three8th);
    // models::bacteria_growth(solver, draw_bitmap_line_chart);
    // models::gompertz(solver, draw_bitmap_line_chart);
    // models::demographic(solver, draw_bitmap_line_chart);
    // models::seasonal_capacity(solver, draw_bitmap_line_chart);
    // models::constant_rate_harvesting(solver, draw_bitmap_line_chart);
    // models::optimal_harvesting(solver, draw_bitmap_line_chart);
    // models::generalist_vs_specialist_predator(solver, draw_bitmap_line_chart);

    let solver = ODESolver2::new(SolverMethod::Classic4);
    // models::insect_pest_control(solver, draw_bitmap_line_chart);
    models::lotka_volterra(solver, draw_bitmap_line_chart);
}
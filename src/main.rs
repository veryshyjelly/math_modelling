use chart::draw_bitmap_line_chart;
use ode_solvers::{ODESolver1, ODESolver2, SolverMethod};

mod chart;
mod models;
mod ode_solvers;

fn main() {
    let solver1 = ODESolver1::new(SolverMethod::Three8th);
    let solver2 = ODESolver2::new(SolverMethod::ForwardEuler);

    models::bacteria_growth(solver1, draw_bitmap_line_chart);
    models::gompertz(solver1, draw_bitmap_line_chart);
    models::demographic(solver1, draw_bitmap_line_chart);
    models::seasonal_capacity(solver1, draw_bitmap_line_chart);
    models::constant_rate_harvesting(solver1, draw_bitmap_line_chart);
    models::optimal_harvesting(solver1, draw_bitmap_line_chart);
    models::generalist_vs_specialist_predator(solver1, draw_bitmap_line_chart);
    models::pest_control1(solver2, draw_bitmap_line_chart);
    models::pest_control2(solver2, draw_bitmap_line_chart);
    models::lotka_volterra(solver2, draw_bitmap_line_chart);
    models::logistic_predator_prey(solver2, draw_bitmap_line_chart);
    models::competition_model(solver2, draw_bitmap_line_chart);
    models::another_competition_model(solver2, draw_bitmap_line_chart);
    models::mutualism1(solver2, draw_bitmap_line_chart);
    models::mutualism2(solver2, draw_bitmap_line_chart);
    models::rabies_pest1(solver2, draw_bitmap_line_chart);
    models::rabies_pest2(solver2, draw_bitmap_line_chart);
    models::rabies_pest3(solver2, draw_bitmap_line_chart);
}
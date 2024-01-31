mod another_competition_model;
mod bacteria_growth;
mod competition_model;
mod constant_rate_harvesting;
mod demographic;
mod generalist_vs_specialist;
mod gompertz;
mod insect_pest_control;
mod logistic_predator_prey;
mod lotka_volterra;
mod mutualism;
mod optimal_harvesting;
mod rabies_pest;
mod seasonal_capacity;

pub use another_competition_model::another_competition_model;
pub use bacteria_growth::bacteria_growth;
pub use competition_model::competition_model;
pub use constant_rate_harvesting::constant_rate_harvesting;
pub use demographic::demographic;
pub use generalist_vs_specialist::generalist_vs_specialist_predator;
pub use gompertz::gompertz;
pub use insect_pest_control::{pest_control1, pest_control2};
pub use logistic_predator_prey::logistic_predator_prey;
pub use lotka_volterra::lotka_volterra;
pub use mutualism::{mutualism1, mutualism2};
pub use optimal_harvesting::optimal_harvesting;
pub use rabies_pest::{rabies_pest1, rabies_pest2, rabies_pest3};
pub use seasonal_capacity::seasonal_capacity;

use ndarray::Array1;

type ChartDrawer = fn(&str, &str, Vec<(Array1<f64>, Array1<f64>, &str)>);
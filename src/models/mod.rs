mod bacteria_growth;
mod constant_rate_harvesting;
mod demographic;
mod generalist_vs_specialist;
mod gompertz;
mod insect_pest_control;
mod lotka_volterra;
mod optimal_harvesting;
mod seasonal_capacity;

pub use bacteria_growth::bacteria_growth;
pub use constant_rate_harvesting::constant_rate_harvesting;
pub use demographic::demographic;
pub use generalist_vs_specialist::generalist_vs_specialist_predator;
pub use gompertz::gompertz;
pub use insect_pest_control::insect_pest_control;
pub use lotka_volterra::lotka_volterra;
use ndarray::Array1;
pub use optimal_harvesting::optimal_harvesting;
pub use seasonal_capacity::seasonal_capacity;

type ChartDrawer = fn(&str, &str, Vec<(Array1<f64>, Array1<f64>, &str)>);
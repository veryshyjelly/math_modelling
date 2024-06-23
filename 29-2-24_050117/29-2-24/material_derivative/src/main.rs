use ndarray::Array1;
use plotters::{
    backend::BitMapBackend,
    chart::ChartBuilder,
    drawing::IntoDrawingArea,
    series::LineSeries,
    style::{RED, WHITE},
};

mod chart;

fn main() {
    let (x_iter, t_iter) = (40, 10);
    let mut uxt: Vec<Vec<(f64, f64)>> = vec![];
    let (x_max, t_max) = (2., 1.);
    let x = Array1::linspace(0., x_max, x_iter);
    let t = Array1::linspace(0., t_max, t_iter);

    let ux0 = |x: f64| 2. * x + 5.;
    let u0t = |t: f64| t + 5.;

    let (delta_t, delta_x) = (t[1], x[1]);

    uxt.push(x.iter().map(|&xi| (xi, ux0(xi))).collect());
    println!("{:?}", uxt[0]);
    for i in 0..t_iter - 1 {
        let mut utxi = vec![(0., u0t(t[i + 1]))];
        for j in 0..x_iter - 2 {
            let mut val =
                uxt[i][j + 1].1 * (1. + (delta_t / delta_x) * (uxt[i][j + 1].1 - uxt[i][j + 2].1));
            if val > 1000. {
                val = 1000.;
            } else if val < -1000.0 {
                val = -1000.0;
            }
            utxi.push((x[j + 1], val))
        }
        utxi.push((x_max, uxt[i][x_iter - 1].1));
        println!("utxi {:?}", utxi);
        uxt.push(utxi);
        println!("calculated {}", i);
    }
    println!("calculation completed");

    let root = BitMapBackend::gif("material_derivative.gif", (800, 600), 200)
        .unwrap()
        .into_drawing_area();

    println!("num frames = {}", uxt.len());

    for i in 0..uxt.len() {
        root.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .caption(
                format!("Material Derivate (n_iter = {})", i),
                ("sans-serif", 20),
            )
            .build_cartesian_2d(0.0..x_max, -50.5..50.)
            .unwrap();
        chart.configure_mesh().draw().unwrap();
        chart
            .draw_series(LineSeries::new(uxt.remove(0), &RED))
            .unwrap();
        println!("drawn {}", i);
        root.present().unwrap();
    }

    println!("Result has been saved");
}

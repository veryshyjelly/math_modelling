use ndarray::{Array1, Array2};
use plotters::{
    backend::{BitMapBackend, SVGBackend},
    chart::ChartBuilder,
    drawing::IntoDrawingArea,
    series::LineSeries,
    style::{RED, WHITE},
};

mod chart;

fn main() {
    let (x_shape, y_shape) = (10, 10);
    let (x, y) = (
        Array1::linspace(0., 1., x_shape),
        Array1::linspace(0., 1., y_shape),
    );
    let mut u = Array2::<f64>::zeros((x_shape, y_shape));
    let mut v = Array2::<f64>::zeros((x_shape, y_shape));
    for i in 0..x_shape {
        for j in 0..y_shape {
            u[(i, j)] = 0.5 + x[i] + y[j];
            v[(i, j)] = 0.4 + x[i] + y[j];
        }
    }
    let (delta_x, delta_y) = (x[1], y[1]);

    let mut u_plot = vec![];
    let mut v_plot = vec![];

    let delta_t = 0.05;
    let t_steps = 20;
    let t = Array1::linspace(0., t_steps as f64 * delta_t, t_steps);

    for time in 0..t_steps {
        println!("u = {:?} \n v = {:?}", u, v);
        let mut u_next = u.clone();
        let mut v_next = v.clone();
        for i in 1..x_shape - 1 {
            for j in 1..y_shape - 1 {
                u_next[(i, j)] = u[(i, j)]
                    - delta_t
                        * (u[(i, j)] * (u[(i, j)] - u[(i - 1, j)]) / delta_x
                            + v[(i, j)] * (u[(i, j)] - u[(i, j - 1)]) / delta_y);
                v_next[(i, j)] = v[(i, j)]
                    - delta_t
                        * (u[(i, j)] * (v[(i, j)] - v[(i - 1, j)]) / delta_x
                            + v[(i, j)] * (v[(i, j)] - v[(i, j - 1)]) / delta_y);
                u_next[(0, j)] = 0.5 + y[j];
                u_next[(x_shape - 1, j)] = 1.5 + y[j];
                v_next[(0, j)] = 0.4 + y[j];
                v_next[(x_shape - 1, j)] = 1.4 + y[j];
            }
            u_next[(i, 0)] = 0.5 + x[i];
            u_next[(i, y_shape - 1)] = 0.5 + x[i] + 1.;
            v_next[(i, 0)] = 0.4 + x[i];
            v_next[(i, y_shape - 1)] = 0.4 + x[i] + 1.;
        }
        u = u_next;
        v = v_next;
        u_plot.push(u[(x_shape - 1, y_shape - 1)]);
        v_plot.push(v[(x_shape / 2, y_shape / 2)]);
        println!("time = {time}");
    }

    println!("{:?}\n{:?}", u_plot, v_plot);

    chart::draw_bitmap_line_chart(
        "2d_flow.svg",
        "",
        vec![
            (&t, Array1::from(u_plot), "u"),
            (&t, Array1::from(v_plot), "v"),
        ],
    );
}

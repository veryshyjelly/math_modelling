use std::time;

use chart::draw_bitmap_line_chart;
use ndarray::Array1;

mod chart;
mod ode_solvers;

fn main() {
    sub_model1();
    sub_model2();
}


fn sub_model1() {
    let start = time::Instant::now();
    let (s0, e0, i0, r0) = (1000., 700., 100., 0.);
    let (tn, n_steps) = (1e4, 1e5 as usize);
    let t = Array1::linspace(0., tn, n_steps);

    let (r, b, a) = (3e-4, 9e-1, 0.5);
    
    let ds_by_dt = |_: f64, S: f64, _e: f64, I: f64, _r: f64| - r * S * I;
    let de_by_dt = |_: f64, S: f64, E: f64, I: f64, _r: f64| r*S*I - b*E;
    let di_by_dt = |_: f64, _s: f64, E: f64, I: f64, _r: f64| b*E - a * I;
    let dr_by_dt = |_: f64, _S: f64, _e: f64, I: f64, _r: f64| a * I;

    let (Sn, En, In, Rn) = forward_euler4(ds_by_dt, de_by_dt, di_by_dt, dr_by_dt, t[1], n_steps, &t, s0, e0, i0, r0);

    println!("S = {} E = {} I = {}", Sn.last().unwrap(), En.last().unwrap(), In.last().unwrap());

    draw_bitmap_line_chart(
        "plots/sub_model1.svg",
        "Sub-Model 1",
        vec![
            (&t, Sn, "S"),
            (&t, En, "E"),
            (&t, In, "I"),
            (&t, Rn, "R"),
        ],
    );

    println!("time taken = {:?}", time::Instant::now() - start);
    
}

fn sub_model2() {
    let start = time::Instant::now();
    let (s0, e0, i0, r0) = (1000., 700., 100., 0.);
    let (tn, n_steps) = (50., 1e5 as usize);
    let t = Array1::linspace(0., tn, n_steps);

    let (r, b, a, c) = (3e-4, 9e-1, 0.5, 0.01);
    
    let ds_by_dt = |_: f64, S: f64, _e: f64, I: f64, _r: f64| - r * S * I;
    let de_by_dt = |_: f64, S: f64, E: f64, I: f64, _r: f64| r*S*I - b*E;
    let di_by_dt = |_: f64, _s: f64, E: f64, I: f64, _r: f64| b*E - a * I;
    let dr_by_dt = |_: f64, _S: f64, _e: f64, I: f64, R: f64| a * I - c *R;

    let (Sn, En, In, Rn) = forward_euler4(ds_by_dt, de_by_dt, di_by_dt, dr_by_dt, t[1], n_steps, &t, s0, e0, i0, r0);

    println!("S = {} E = {} I = {}", Sn.last().unwrap(), En.last().unwrap(), In.last().unwrap());

    draw_bitmap_line_chart(
        "plots/sub_model2.svg",
        "Sub-Model 2",
        vec![
            (&t, Sn, "S"),
            (&t, En, "E"),
            (&t, In, "I"),
            (&t, Rn, "R"),
        ],
    );

    println!("time taken = {:?}", time::Instant::now() - start);
    
}

pub fn forward_euler4<F, G, H, I>(
    f1: F,
    f2: G,
    f3: H,
    f4: I,
    h: f64,
    n_steps: usize,
    t: &Array1<f64>,
    y1_0: f64,
    y2_0: f64,
    y3_0: f64,
    y4_0: f64,
) -> (Array1<f64>, Array1<f64>, Array1<f64>, Array1<f64>)
where
    F: Fn(f64, f64, f64, f64, f64) -> f64,
    G: Fn(f64, f64, f64, f64, f64) -> f64,
    H: Fn(f64, f64, f64, f64, f64) -> f64,
    I: Fn(f64, f64, f64, f64, f64) -> f64,
{
    let mut y1 = Array1::zeros(n_steps + 1);
    let mut y2 = Array1::zeros(n_steps + 1);
    let mut y3 = Array1::zeros(n_steps + 1);
    let mut y4 = Array1::zeros(n_steps + 1);
    y1[0] = y1_0;
    y2[0] = y2_0;
    y3[0] = y3_0;
    y4[0] = y4_0;
    for i in 0..n_steps {
        let k1 = h * f1(t[i], y1[i], y2[i], y3[i], y4[i]);
        y1[i + 1] = y1[i] + k1;
        let l1 = h * f2(t[i], y1[i + 1], y2[i], y3[i], y4[i]);
        y2[i + 1] = y2[i] + l1;
        let m1 = h * f3(t[i], y1[i + 1], y2[i + 1], y3[i], y4[i]);
        y3[i + 1] = y3[i] + m1;
        let n1 = h * f4(t[i], y1[i + 1], y2[i + 1], y3[i], y4[i]);
        y4[i + 1] = y4[i] + n1;
    }
    (y1, y2, y3, y4)
}

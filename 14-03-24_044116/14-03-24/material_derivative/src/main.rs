use ndarray::Array1;

mod chart;

fn main() {
    let (x_iter, t_iter) = (40, 10);
    let mut uxt1: Vec<Vec<(f64, f64)>> = vec![];
    let mut uxt2: Vec<Vec<(f64, f64)>> = vec![];
    let mut uxt3: Vec<Vec<(f64, f64)>> = vec![];
    // let mut uxt: Vec<Vec<(f64, f64)>> = vec![];

    let (x_max, t_max) = (2., 3.);
    let x = Array1::linspace(0., x_max, x_iter);
    let t = Array1::linspace(0., t_max, t_iter);

    let ux0 = |x: f64| 2. * x + 5.;
    let u0t = |t: f64| t + 5.;

    let (delta_t, delta_x) = (t[1], x[1]);

    uxt1.push(x.iter().map(|&xi| (xi, ux0(xi))).collect());
    uxt2.push(x.iter().map(|&xi| (xi, ux0(xi))).collect());
    uxt3.push(x.iter().map(|&xi| (xi, ux0(xi))).collect());
    // println!("{:?}", uxt1[0]);

    for i in 0..t_iter - 1 {
        let mut uxti1 = vec![(0., u0t(t[i + 1]))];
        let mut uxti2 = vec![(0., u0t(t[i + 1]))];
        let mut uxti3 = vec![(0., u0t(t[i + 1]))];
        for j in 0..x_iter - 2 {
            let MU = 0.001;

            let mut val1 = uxt1[i][j + 1].1
                * (1.
                    + (delta_t / delta_x) * (uxt1[i][j + 1].1 - uxt1[i][j + 2].1)
                    + MU * (delta_t / delta_x.powi(2))
                        * (2. * uxt1[i][j + 1].1 - uxt1[i][j + 2].1 - uxt1[i][j].1));
            if val1 > 1000. {
                val1 = 1000.;
            } else if val1 < -1000.0 {
                val1 = -1000.0;
            }
            uxti1.push((x[j + 1], val1));

            let MU = 0.05;

            let mut val2 = uxt2[i][j + 1].1
                * (1.
                    + (delta_t / delta_x) * (uxt2[i][j + 1].1 - uxt2[i][j + 2].1)
                    + MU * (delta_t / delta_x.powi(2))
                        * (2. * uxt2[i][j + 1].1 - uxt2[i][j + 2].1 - uxt2[i][j].1));
            if val2 > 1000. {
                val2 = 1000.;
            } else if val2 < -1000.0 {
                val2 = -1000.0;
            }
            uxti2.push((x[j + 1], val2));

            let MU = 0.01;
            let mut val3 = uxt3[i][j + 1].1
                * (1.
                    + (delta_t / delta_x) * (uxt3[i][j + 1].1 - uxt3[i][j + 2].1)
                    + MU * (delta_t / delta_x.powi(2))
                        * (2. * uxt3[i][j + 1].1 - uxt3[i][j + 2].1 - uxt3[i][j].1));
            if val3 > 1000. {
                val3 = 1000.;
            } else if val3 < -1000.0 {
                val3 = -1000.0;
            }
            uxti3.push((x[j + 1], val3))
        }
        uxti1.push((x_max, uxt1[i][x_iter - 1].1));
        uxti2.push((x_max, uxt2[i][x_iter - 1].1));
        uxti3.push((x_max, uxt3[i][x_iter - 1].1));
        // println!("utxi {:?}", uxti1);
        uxt1.push(uxti1);
        uxt2.push(uxti2);
        uxt3.push(uxti3);
        println!("calculated {}", i);
    }

    let mut u_plot1 = vec![];
    let mut u_plot2 = vec![];
    let mut u_plot3 = vec![];
    for i in 0..t_iter - 1 {
        u_plot1.push(uxt1[i][x_iter / 2].1);
        u_plot2.push(uxt2[i][x_iter / 2].1);
        u_plot3.push(uxt3[i][x_iter / 2].1);
    }

    println!("{:?}", u_plot1);

    chart::draw_bitmap_line_chart(
        "1d_flow.svg",
        "",
        vec![
            (&t, Array1::from(u_plot1), "u = 0.001"),
            (&t, Array1::from(u_plot2), "u = 0.05"),
            (&t, Array1::from(u_plot3), "u = 0.01"),
        ],
    );
}

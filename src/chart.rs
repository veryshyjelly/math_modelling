use ndarray::Array1;
use plotters::backend::BitMapBackend;
use plotters::chart::{ChartBuilder, LabelAreaPosition};
use plotters::drawing::IntoDrawingArea;
use plotters::prelude::IntoFont;
use plotters::series::LineSeries;
use plotters::style::{RGBColor, BLUE, CYAN, GREEN, MAGENTA, RED, WHITE, YELLOW};

const COLORS: [&RGBColor; 6] = [&RED, &GREEN, &BLUE, &YELLOW, &CYAN, &MAGENTA];

pub fn draw_bitmap_line_chart(
    x_range: (f64, f64),
    y_range: (f64, f64),
    file_name: &str,
    caption: &str,
    lines: Vec<(Array1<f64>, Array1<f64>)>,
) {
    let root_drawing_area = BitMapBackend::new(file_name, (1024, 768)).into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root_drawing_area)
        .caption(caption, ("sans-serif", 40).into_font())
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(x_range.0..x_range.1, y_range.0..y_range.1)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    for (i, (x, y)) in lines.into_iter().enumerate() {
        chart
            .draw_series(LineSeries::new(x.into_iter().zip(y), COLORS[i % 6]))
            .unwrap();
    }
}
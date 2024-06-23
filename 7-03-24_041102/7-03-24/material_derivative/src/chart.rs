use ndarray::Array1;
use plotters::backend::SVGBackend;
use plotters::chart::{ChartBuilder, LabelAreaPosition};
use plotters::drawing::IntoDrawingArea;
use plotters::element::PathElement;
use plotters::prelude::IntoFont;
use plotters::series::LineSeries;
use plotters::style::{Color, RGBColor, BLACK, BLUE, CYAN, GREEN, MAGENTA, RED, WHITE, YELLOW};

const COLORS: [&RGBColor; 6] = [&RED, &GREEN, &BLUE, &YELLOW, &CYAN, &MAGENTA];

pub fn draw_bitmap_line_chart(
    file_name: &str,
    caption: &str,
    lines: Vec<(&Array1<f64>, Array1<f64>, &str)>,
) {
    let (x_min, x_max) = (
        lines
            .iter()
            .map(|(x, _, _)| x.iter().rfold(0., |x, y| f64::min(x, *y)))
            .rfold(0., |x, y| f64::min(x, y)),
        lines
            .iter()
            .map(|(x, _, _)| x.iter().rfold(0., |x, y| f64::max(x, *y)))
            .rfold(0., |x, y| f64::max(x, y)),
    );

    println!("x max = {} min = {}", x_max, x_min);

    let (y_min, y_max) = (
        lines
            .iter()
            .map(|(_, y, _)| y.iter().rfold(0., |x, y| f64::min(x, *y)))
            .rfold(0., |x, y| f64::min(x, y)),
        lines
            .iter()
            .map(|(_, y, _)| y.iter().rfold(0., |x, y| f64::max(x, *y)))
            .rfold(0., |x, y| f64::max(x, y)),
    );

    println!("y max = {} min = {}", y_max, y_min);

    let root_drawing_area = SVGBackend::new(file_name, (1024, 768)).into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root_drawing_area)
        .caption(caption, ("sans-serif", 40).into_font())
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    for (i, (x, y, label)) in lines.into_iter().enumerate() {
        println!("label: {}", label);
        let color = COLORS[i % 6];
        chart
            .draw_series(LineSeries::new(x.to_owned().into_iter().zip(y), color))
            .unwrap()
            .label(label)
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.clone()));
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    root_drawing_area.present().unwrap();
}

use ndarray::Array1;

#[derive(Copy, Clone)]
pub enum SolverMethod {
    ForwardEuler,
    ExplicitMidpoint,
    Heun2,
    Ralston2,
    Kutta3,
    Wray3,
    Ralston3,
    SSPRK3,
    Classic4,
    Three8th,
}

#[derive(Copy, Clone)]
pub struct ODESolver1 {
    method: SolverMethod,
}

impl ODESolver1 {
    pub fn new(method: SolverMethod) -> Self {
        Self { method }
    }

    pub fn solve<F>(self, f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        match self.method {
            SolverMethod::ForwardEuler => ODESolver1::forward_euler(f, h, n_steps, t, y_0),
            SolverMethod::ExplicitMidpoint => ODESolver1::explicit_midpoint(f, h, n_steps, t, y_0),
            SolverMethod::Heun2 => ODESolver1::heun_s2(f, h, n_steps, t, y_0),
            SolverMethod::Ralston2 => ODESolver1::ralston_s2(f, h, n_steps, t, y_0),
            SolverMethod::Kutta3 => ODESolver1::kutta_s3(f, h, n_steps, t, y_0),
            SolverMethod::Wray3 => ODESolver1::wray_s3(f, h, n_steps, t, y_0),
            SolverMethod::Ralston3 => ODESolver1::ralston_s3(f, h, n_steps, t, y_0),
            SolverMethod::SSPRK3 => ODESolver1::ssprk3(f, h, n_steps, t, y_0),
            SolverMethod::Classic4 => ODESolver1::classic4(f, h, n_steps, t, y_0),
            SolverMethod::Three8th => ODESolver1::three_8th(f, h, n_steps, t, y_0),
        }
    }

    pub fn forward_euler<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            y[i + 1] = y[i] + h * k1;
        }
        y
    }

    pub fn explicit_midpoint<F>(
        f: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y_0: f64,
    ) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + h / 2., y[i] + k1 / 2.);
            y[i + 1] = y[i] + h * k2;
        }
        y
    }

    pub fn heun_s2<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + h, y[i] + k1);
            y[i + 1] = y[i] + h * (k1 + k2) / 2.;
        }
        y
    }

    pub fn ralston_s2<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + 2. * h / 3., y[i] + 2. * k1 / 3.);
            y[i + 1] = y[i] + h * (k1 + 3. * k2) / 4.;
        }
        y
    }

    pub fn kutta_s3<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + h / 2., y[i] + k1 / 2.);
            let k3 = f(t[i] + h, y[i] + 2. * k2 - k1);
            y[i + 1] = y[i] + h * (k1 + 4. * k2 + k3) / 6.;
        }
        y
    }

    pub fn wray_s3<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + 8. * h / 15., y[i] + 8. * k1 / 15.);
            let k3 = f(t[i] + 2. * h / 3., y[i] + k1 / 4. + 5. * k2 / 12.);
            y[i + 1] = y[i] + h * (k1 + 3. * k3) / 4.;
        }
        y
    }

    pub fn ralston_s3<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + h / 2., y[i] + k1 / 2.);
            let k3 = f(t[i] + 3. * h / 4., y[i] + 3. * k2 / 4.);
            y[i + 1] = y[i] + h * (2. * k1 + 3. * k2 + 4. * k3) / 9.;
        }
        y
    }

    pub fn ssprk3<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + h, y[i] + k1);
            let k3 = f(t[i] + h / 2., y[i] + k1 / 4. + k2 / 4.);
            y[i + 1] = y[i] + h * (k1 + k2 + 4. * k3) / 6.;
        }
        y
    }

    pub fn classic4<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + h / 2., y[i] + k1 / 2.);
            let k3 = f(t[i] + h / 2., y[i] + k2 / 2.);
            let k4 = f(t[i] + h, y[i] + k3);
            y[i + 1] = y[i] + h * (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
        }
        y
    }

    pub fn three_8th<F>(f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut y = Array1::zeros(n_steps + 1);
        y[0] = y_0;
        for i in 0..n_steps {
            let k1 = f(t[i], y[i]);
            let k2 = f(t[i] + h / 3., y[i] + k1 / 3.);
            let k3 = f(t[i] + 2. * h / 3., y[i] + k2 - k1 / 3.);
            let k4 = f(t[i] + h, y[i] + k1 - k2 + k3);
            y[i + 1] = y[i] + h * (k1 + 3. * k2 + 3. * k3 + k4) / 8.;
        }
        y
    }
}

#[derive(Copy, Clone)]
pub struct ODESolver2 {
    method: SolverMethod,
}

impl ODESolver2 {
    pub fn new(method: SolverMethod) -> Self {
        Self { method }
    }

    pub fn solve<F>(
        self,
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        match self.method {
            SolverMethod::ForwardEuler => {
                ODESolver2::forward_euler(f1, f2, h, n_steps, t, y1_0, y2_0)
            }
            SolverMethod::ExplicitMidpoint => {
                ODESolver2::explicit_midpoint(f1, f2, h, n_steps, t, y1_0, y2_0)
            }
            SolverMethod::Heun2 => ODESolver2::heun_s2(f1, f2, h, n_steps, t, y1_0, y2_0),
            SolverMethod::Ralston2 => ODESolver2::ralston_s2(f1, f2, h, n_steps, t, y1_0, y2_0),
            SolverMethod::Kutta3 => ODESolver2::kutta_s3(f1, f2, h, n_steps, t, y1_0, y2_0),
            SolverMethod::Wray3 => ODESolver2::wray_s3(f1, f2, h, n_steps, t, y1_0, y2_0),
            SolverMethod::Ralston3 => ODESolver2::ralston_s3(f1, f2, h, n_steps, t, y1_0, y2_0),
            SolverMethod::SSPRK3 => ODESolver2::ssprk3(f1, f2, h, n_steps, t, y1_0, y2_0),
            SolverMethod::Classic4 => ODESolver2::classic4(f1, f2, h, n_steps, t, y1_0, y2_0),
            SolverMethod::Three8th => ODESolver2::three_8th(f1, f2, h, n_steps, t, y1_0, y2_0),
        }
    }

    pub fn forward_euler<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            y1[i + 1] = y1[i] + h * k1;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            y2[i + 1] = y2[i] + h * k1;
        }
        (y1, y2)
    }

    pub fn explicit_midpoint<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + h / 2., y1[i] + k1 / 2., y2[i]);
            y1[i + 1] = y1[i] + h * k2;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            let k2 = f2(t[i] + h / 2., y1[i + 1], y2[i] + k1 / 2.);
            y2[i + 1] = y2[i] + h * k2;
        }
        (y1, y2)
    }

    pub fn heun_s2<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + h, y1[i] + k1, y2[i]);
            y1[i + 1] = y1[i] + h * (k1 + k2) / 2.;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            let k2 = f2(t[i] + h, y1[i + 1], y2[i] + k1);
            y2[i + 1] = y2[i] + h * (k1 + k2) / 2.;
        }
        (y1, y2)
    }

    pub fn ralston_s2<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + 2. * h / 3., y1[i] + 2. * k1 / 3., y2[i]);
            y1[i + 1] = y1[i] + h * (k1 + 3. * k2) / 4.;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            let k2 = f2(t[i] + 2. * h / 3., y1[i + 1], y2[i] + 2. * k1 / 3.);
            y2[i + 1] = y2[i] + h * (k1 + 3. * k2) / 4.;
        }
        (y1, y2)
    }

    pub fn kutta_s3<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + h / 2., y1[i] + k1 / 2., y2[i]);
            let k3 = f1(t[i] + h, y1[i] + 2. * k2 - k1, y2[i]);
            y1[i + 1] = y1[i] + h * (k1 + 4. * k2 + k3) / 6.;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            let k2 = f2(t[i] + h / 2., y1[i + 1], y2[i] + k1 / 2.);
            let k3 = f2(t[i] + h, y1[i + 1], y2[i] + 2. * k2 - k1);
            y2[i + 1] = y2[i] + h * (k1 + 4. * k2 + k3) / 6.;
        }
        (y1, y2)
    }

    pub fn wray_s3<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + 8. * h / 15., y1[i] + 8. * k1 / 15., y2[i]);
            let k3 = f1(t[i] + 2. * h / 3., y1[i] + k1 / 4. + 5. * k2 / 12., y2[i]);
            y1[i + 1] = y1[i] + h * (k1 + 3. * k3) / 4.;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            let k2 = f2(t[i] + 8. * h / 15., y1[i + 1], y2[i] + 8. * k1 / 15.);
            let k3 = f2(
                t[i] + 2. * h / 3.,
                y1[i + 1],
                y2[i] + k1 / 4. + 5. * k2 / 12.,
            );
            y2[i + 1] = y2[i] + h * (k1 + 3. * k3) / 4.;
        }

        (y1, y2)
    }

    pub fn ralston_s3<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + h / 2., y1[i] + k1 / 2., y2[i]);
            let k3 = f1(t[i] + 3. * h / 4., y1[i] + 3. * k2 / 4., y2[i]);
            y1[i + 1] = y1[i] + h * (2. * k1 + 3. * k2 + 4. * k3) / 9.;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            let k2 = f2(t[i] + h / 2., y1[i + 1], y2[i] + k1 / 2.);
            let k3 = f2(t[i] + 3. * h / 4., y1[i + 1], y2[i] + 3. * k2 / 4.);
            y2[i + 1] = y2[i] + h * (2. * k1 + 3. * k2 + 4. * k3) / 9.;
        }

        (y1, y2)
    }

    pub fn ssprk3<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + h, y1[i] + k1, y2[i]);
            let k3 = f1(t[i] + h / 2., y1[i] + k1 / 4. + k2 / 4., y2[i]);
            y1[i + 1] = y1[i] + h * (k1 + k2 + 4. * k3) / 6.;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            let k2 = f2(t[i] + h, y1[i + 1], y2[i] + k1);
            let k3 = f2(t[i] + h / 2., y1[i + 1], y2[i] + k1 / 4. + k2 / 4.);
            y2[i + 1] = y2[i] + h * (k1 + k2 + 4. * k3) / 6.;
        }

        (y1, y2)
    }

    pub fn classic4<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + h / 2., y1[i] + k1 / 2., y2[i]);
            let k3 = f1(t[i] + h / 2., y1[i] + k2 / 2., y2[i]);
            let k4 = f1(t[i] + h, y1[i] + k3, y2[i]);
            y1[i + 1] = y1[i] + h * (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
            let k1 = f2(t[i], y1[i + 1], y2[i]);
            let k2 = f2(t[i] + h / 2., y1[i + 1], y2[i] + k1 / 2.);
            let k3 = f2(t[i] + h / 2., y1[i + 1], y2[i] + k2 / 2.);
            let k4 = f2(t[i] + h, y1[i + 1], y2[i] + k3);
            y2[i + 1] = y2[i] + h * (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
        }
        (y1, y2)
    }

    pub fn three_8th<F>(
        f1: F,
        f2: F,
        h: f64,
        n_steps: usize,
        t: &Array1<f64>,
        y1_0: f64,
        y2_0: f64,
    ) -> (Array1<f64>, Array1<f64>)
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let mut y1 = Array1::zeros(n_steps + 1);
        let mut y2 = Array1::zeros(n_steps + 1);
        y1[0] = y1_0;
        y2[0] = y2_0;
        for i in 0..n_steps {
            let k1 = f1(t[i], y1[i], y2[i]);
            let k2 = f1(t[i] + h / 3., y1[i] + k1 / 3., y2[i]);
            let k3 = f1(t[i] + 2. * h / 3., y1[i] + k2 - k1 / 3., y2[i]);
            let k4 = f1(t[i] + h, y1[i] + k1 - k2 + k3, y2[i]);
            y1[i + 1] = y1[i] + h * (k1 + 3. * k2 + 3. * k3 + k4) / 8.;
            let k1 = f2(t[i], y1[i], y2[i]);
            let k2 = f2(t[i] + h / 3., y1[i], y2[i] + k1 / 3.);
            let k3 = f2(t[i] + 2. * h / 3., y1[i], y2[i] + k2 - k1 / 3.);
            let k4 = f2(t[i] + h, y1[i], y2[i] + k1 - k2 + k3);
            y1[i + 1] = y1[i] + h * (k1 + 3. * k2 + 3. * k3 + k4) / 8.;
        }

        (y1, y2)
    }
}
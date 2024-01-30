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
pub struct OdeSolver {
    method: SolverMethod,
}

impl OdeSolver {
    pub fn new(method: SolverMethod) -> OdeSolver {
        OdeSolver { method }
    }

    pub fn solve<F>(self, f: F, h: f64, n_steps: usize, t: &Array1<f64>, y_0: f64) -> Array1<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        match self.method {
            SolverMethod::ForwardEuler => OdeSolver::forward_euler(f, h, n_steps, t, y_0),
            SolverMethod::ExplicitMidpoint => OdeSolver::explicit_midpoint(f, h, n_steps, t, y_0),
            SolverMethod::Heun2 => OdeSolver::heun_s2(f, h, n_steps, t, y_0),
            SolverMethod::Ralston2 => OdeSolver::ralston_s2(f, h, n_steps, t, y_0),
            SolverMethod::Kutta3 => OdeSolver::kutta_s3(f, h, n_steps, t, y_0),
            SolverMethod::Wray3 => OdeSolver::wray_s3(f, h, n_steps, t, y_0),
            SolverMethod::Ralston3 => OdeSolver::ralston_s3(f, h, n_steps, t, y_0),
            SolverMethod::SSPRK3 => OdeSolver::ssprk3(f, h, n_steps, t, y_0),
            SolverMethod::Classic4 => OdeSolver::classic4(f, h, n_steps, t, y_0),
            SolverMethod::Three8th => OdeSolver::three_8th(f, h, n_steps, t, y_0),
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
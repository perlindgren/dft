// naive dft
// https://en.wikipedia.org/wiki/Discrete_Fourier_transform

#![allow(non_snake_case)]
use num_complex::*;
use std::f32::consts::TAU;

// x being a real valued vector
pub fn dft(x: &[f32]) -> Vec<Complex<f32>> {
    let N = x.len();
    let Nf = N as f32;
    let mut X = vec![];
    for k in 0..N {
        let kf = k as f32;
        let mut Xk = c32(0.0, 0.0);
        for (n, x) in x.iter().enumerate() {
            let nf = n as f32;
            let Xf = *x * c32(0.0, -nf * kf * TAU / Nf).exp();
            Xk += Xf;
        }
        X.push(Xk)
    }
    X
}

// x being a real valued vector
pub fn idft_c(X: &[Complex<f32>]) -> Vec<Complex<f32>> {
    let N = X.len();
    let Nf = N as f32;
    let mut x = vec![];
    for k in 0..N {
        let kf = k as f32;
        let mut xk = c32(0.0, 0.0);
        for (n, X) in X.iter().enumerate() {
            let nf = n as f32;
            let x = *X * c32(0.0, nf * kf * TAU / Nf).exp();
            xk += x;
        }
        x.push(xk / Nf);
    }
    x
}

pub fn idft_r(X: &[Complex<f32>]) -> Vec<f32> {
    let N = X.len();
    let Nf = N as f32;
    let mut x = vec![];
    for k in 0..N {
        let kf = k as f32;
        let mut xk = 0.0;
        for (n, X) in X.iter().enumerate() {
            let nf = n as f32;
            let x = ((*X) * c32(0.0, nf * kf * TAU / Nf).exp()).re();
            xk += x;
        }
        x.push(xk / Nf);
    }
    x
}

pub fn gen_sin(fs: f32, t: f32, f: f32) -> Vec<f32> {
    let mut x = vec![];
    for i in 0..(fs * t) as usize {
        x.push((f * TAU * i as f32 / fs).sin());
    }
    x
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_dft_c() {
        let X = dft(&[1.0, 0.0, 2.0]);
        println!("--- dft ---\n{X:?}");

        let x = idft_c(X.as_ref());
        println!("--- idft ---\n{x:?}");
    }

    #[test]
    fn test_dft_r() {
        let X = dft(&[1.0, 0.0, 2.0]);
        println!("--- dft ---\n{X:?}");

        let x = idft_r(X.as_ref());
        println!("--- idft ---\n{x:?}");
    }

    #[test]
    fn test_dft_sine() {}
}

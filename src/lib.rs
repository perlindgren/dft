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
            println!("k {k}, n {n}, x {x}, Xf {Xf}");
            Xk += Xf;
        }
        X.push(Xk)
    }
    X
}

// x being a real valued vector
pub fn idft(X: &[Complex<f32>]) -> Vec<Complex<f32>> {
    let N = X.len();
    let Nf = N as f32;
    let mut x = vec![];
    for k in 0..N {
        println!("k {k}");
        let kf = k as f32;
        let mut xk = c32(0.0, 0.0);
        for (n, X) in X.iter().enumerate() {
            let nf = n as f32;
            let x = *X * c32(0.0, nf * kf * TAU / Nf).exp();
            println!("k {k}, n {n}, x {X}, Xf {x}");
            xk += x;
        }
        x.push(xk / Nf);
    }
    x
}

#[cfg(test)]
mod test {
    use std::process::id;

    use super::*;

    #[test]
    fn test_dft() {
        let X = dft(&[1.0, 0.0, 2.0]);
        println!("--- dft ---\n{X:?}");

        let x = idft(X.as_ref());
        println!("--- idft ---\n{x:?}");
    }
}

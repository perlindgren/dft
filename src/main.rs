use dft::*;

fn main() {
    let x = gen_sin(1000.0, 1.0, 100.0);

    let fft_x = dft(x.as_ref());
    let ix = idft_r(fft_x.as_ref());
    println!("{ix:?}");
}

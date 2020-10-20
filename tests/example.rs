#![allow(unused_imports)]
use indoc::indoc;
use pretty_assertions::{assert_eq, assert_ne};

use itertools_num::linspace;
use lbm_blood_1d::boundary_conditions::*;
use lbm_blood_1d::simulation::*;
use lbm_blood_1d::simulation_parsing::*;

macro_rules! assert_delta {
    ($x:expr, $y:expr, $d:expr) => {
        if !($x - $y < $d || $y - $x < $d) {
            panic!();
        }
    };
}

#[test]
fn test_standardize_inlet() {
    let data = SimulationParsing::read_inlet(&"vascularNetworks/adan56_inlet.dat".to_string());

    let inlet = Inlet::new(&data);

    println!(
        "{} {}",
        inlet.flow_at_time(0.8277777777777778 + (1.1 * 10.0)),
        -5.840427802383133e-06
    );

    println!(
        "{} {}",
        inlet.flow_at_time(1.09999),
        -5.239489231024034121e-07
    );

    println!(
        "{} {}",
        inlet.flow_at_time(0.0000001),
        -5.239489231023915536e-07
    );
    // );
    // println!(
    //     "{} {} {}",
    //     base[T % (base.len() - 1)],
    //     extended[T],
    //     -2.329872060883883496e-05
    // );
    //
    // let T = (((1.1 * 80.0) + 8.222222222222222987e-01) / dt).round() as usize;
    // assert_delta!(
    //     base[T % (base.len() - 1)],
    //     -6.037601030251277555e-06,
    //     0.0000001f64
    // );
    // println!(
    //     "{} {} {}",
    //     base[T % (base.len() - 1)],
    //     extended[T],
    //     -6.037601030251277555e-06
    // );
    //
    /* println!("{} {}", base.len(), extended.len()); */
}

#![allow(unused_imports)]
#![allow(non_snake_case)]
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

macro_rules! expand_names {
    ( $last:ident ) => {
        $last
    };

    ( $head:ident $(, $tail:ident)* ) => {
        (expand_names!($($tail),*), $head)
    };
}

macro_rules! zip_map {
    ( $name:expr, |($head:ident $(, $tail:ident)*)| $body:expr ) => {
        izip!($name, $head $(, $tail)+)
            .map(|expand_names!($head $(, $tail)+)| $body ).collect()
    };
}

macro_rules! izip {
    ( $name: expr, $last:ident ) => {
        $name.$last.iter()
    };

    ( $name: expr, $head:ident $(, $tail:ident)* ) => {
        izip!($name, $($tail),*).zip($name.$head.iter())
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

    println!("{} {}", inlet.flow_at_time(1.09999), -5.239489231024034121e-07);

    println!("{} {}", inlet.flow_at_time(0.0000001), -5.239489231023915536e-07);

    let sim = Simulation::new("vascularNetworks/bifur.json");
    sim.vessels.iter().for_each(|v| {
        let owo: Vec<f64> = zip_map!(v.cells, |(A, u, f0)| f0 + 0.0);
        println!("aaa {} {}", 5, owo[5]);
    });

    for v in sim.vessels {
        let mut owo = vec![0.0f64; v.x_dim];
        for i in 0..v.x_dim {
            owo[i] = v.cells.f0[i];
        }
        println!("bbb {} {}", 5, owo[5]);
    }

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

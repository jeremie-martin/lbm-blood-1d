// https://bheisler.github.io/criterion.rs/book/getting_started.html

#![allow(non_snake_case)]

use criterion::*;
use lbm_blood_1d::constants::*;
use lbm_blood_1d::lbm_algorithm::*;
use lbm_blood_1d::simulation::*;
use lbm_blood_1d::simulation_parsing::*;
use lbm_blood_1d::vessels::*;
use std::collections::VecDeque;
use std::ptr;

macro_rules! expand_names {
    ( $last:ident ) => {
        $last
    };

    ( $head:ident $(, $tail:ident)* ) => {
        (expand_names!($($tail),*), $head)
    };
}

macro_rules! expand_iter {
    ( $name: expr, $last:ident ) => {
        $name.$last.iter()
    };

    ( $name: expr, $head:ident $(, $tail:ident)* ) => {
        expand_iter!($name, $($tail),*).zip($name.$head.iter())
    };
}

macro_rules! zip_map {
    // n = 1
    ( $name:expr, |$head:ident| $body:expr ) => {
        $name.$head.iter().map(|$head| $body ).collect()
    };

    // n > 1
    ( $name:expr, |($head:ident $(, $tail:ident)*)| $body:expr ) => {
        expand_iter!($name, $head $(, $tail)*)
            .map(|expand_names!($head $(, $tail)*)| $body ).collect()
    };
}

macro_rules! expand_iter_mut {
    ( $name: expr, $last:ident ) => {
        (&$name.$last).into_iter()
    };

    ( $name:expr, mut $last:ident ) => {
        (&mut $name.$last).into_iter()
    };

    ( $name: expr, mut $head:ident, $($($tail:ident)+),+) => {
        expand_iter_mut!($name $(, $($tail)+)*).zip(&mut $name.$head)
    };

    ( $name: expr, $head:ident, $($($tail:ident)+),+) => {
        expand_iter_mut!($name $(, $($tail)+)*).zip(&$name.$head)
    };
}

macro_rules! expand_names_mut {
    ( $last:ident ) => {
        $last
    };

    ( mut $last:ident ) => {
        $last
    };

    ( mut $head:ident $(, $($tail:ident)+)* ) => {
        (expand_names_mut!($($($tail)+),*), $head)
    };

    ( $head:ident $(, $($tail:ident)+)* ) => {
        (expand_names_mut!($($($tail)+),*), $head)
    };
}

macro_rules! zip_for_each {
    // n = 1
    ( $name:expr, |$head:ident| $body:expr ) => {
        (&mut $name.$head).iter_mut().map(|$head| $body).collect()
    };

    // n > 1
    ( $name:expr, |($($($list:ident)+),+)| $body:expr ) => {
        expand_iter_mut!($name, $($($list)+),+)
        .for_each(|expand_names_mut!($($($list)+),+)| $body )
    };
}

macro_rules! zip_for_each_enumerate {
    // n = 1
    ( $name:expr, |($idx:ident, $head:ident)| $body:expr ) => {
        (&mut $name.$head).iter_mut().enumerate().map(|$idx, $head| $body).collect()
    };

    // n > 1
    ( $name:expr, |($idx:ident, $($($list:ident)+),+)| $body:expr ) => {
        expand_iter_mut!($name, $($($list)+),+).enumerate()
        .for_each(|($idx, expand_names2!($($($list)+),+))| $body )
    };
    // n > 1
    ( $name:expr, |($idx:ident, $($($list:ident)+),+)| $body:expr ) => {
        expand_iter_mut!($name, $($($list)+),+).enumerate()
        .for_each(|($idx, expand_names2!($($($list)+),+))| $body )
    };
}

macro_rules! zip_map_enumerate {
    // n = 1
    ( $name:expr, |($idx:ident, $head:ident)| $body:expr ) => {
        $name.$head.iter().enumerate().map(|($idx, $head)| $body ).collect()
    };

    // n > 1
    ( $name:expr, |($idx:ident, $head:ident $(, $tail:ident)*)| $body:expr ) => {
        expand_iter!($name, $head $(, $tail)*).enumerate()
            .map(|($idx, expand_names!($head $(, $tail)*))| $body ).collect()
    };
}

macro_rules! add_prefix {
    ( $prefix:expr, $last:ident ) => {
        $prefix.$last
    };

    ( $prefix:expr, $head:ident $(, $tail:ident)* ) => {
        $prefix.$head, add_prefix!($prefix, $($tail),*)
    };
}

macro_rules! izip2 {
    // @closure creates a tuple-flattening closure for .map() call. usage:
    // @closure partial_pattern => partial_tuple , rest , of , iterators
    // eg. izip!( @closure ((a, b), c) => (a, b, c) , dd , ee )
    ( @closure $p:pat => $tup:expr ) => {
        |$p| $tup
    };

    // The "b" identifier is a different identifier on each recursion level thanks to hygiene.
    ( @closure $p:pat => ( $($tup:tt)* ) , $_iter:expr $( , $tail:expr )* ) => {
        izip2!(@closure ($p, b) => ( $($tup)*, b ) $( , $tail )*)
    };

    // unary
    ($first:expr $(,)*) => {
        $first.into_iter()
    };

    // binary
    ($first:expr, $second:expr $(,)*) => {
        izip2!($first)
            .zip($second)
    };

    // n-ary where n > 2
    ( $first:expr $( , $rest:expr )* $(,)* ) => {
        izip!($first)
            $(
                .zip($rest)
            )*
            .map(
                izip2!(@closure a => (a) $( , $rest )*)
            )
    };
}

macro_rules! izip {
    // @closure creates a tuple-flattening closure for .map() call. usage:
    // @closure partial_pattern => partial_tuple , rest , of , iterators
    // eg. izip!( @closure ((a, b), c) => (a, b, c) , dd , ee )
    ( @closure $p:pat => $tup:ident) => {
        |$p| $tup
    };

    // The "b" identifier is a different identifier on each recursion level thanks to hygiene.
    ( @closure $p:pat => ( $($tup:tt)* ) , $_iter:ident $( , $tail:ident )* ) => {
        izip!(@closure ($p, b) => ( $($tup)*, b ) $( , $tail )*)
    };

    // unary
    ($first:expr $(,)*) => {
        $first.into_iter()
    };

    // binary
    ($first:expr, $second:expr $(,)*) => {
        izip!($first)
            .zip($second)
    };

    // n-ary where n > 2
    ($prefix:expr, $first:ident $( , $rest:ident )* $(,)* ) => {
        izip!($prefix.$first)
            $(
                .zip($prefix.$rest)
            )*
            .map(
                izip!(@closure a => (a) $( , $rest )*)
            )
    };
}

pub fn functional<T: Compute>(sim: &mut Simulation<T>) {
    sim.vessels.iter_mut().for_each(|v| {
        // for i in 0..v.x_dim {
        //     let P_derivative = (v.cells.beta[i]) / (2.0 * (v.cells.A[i] * v.cells.A0[i]).sqrt());
        //     let H_derivative = (v.cells.A[i] / v.consts.rho) * P_derivative;
        //
        //     let A_derivative = match i {
        //         0 => (v.cells.A[i + 1] - v.cells.A[i]) / v.consts.dx,
        //         i if (i == v.x_last) => (v.cells.A[i] - v.cells.A[i - 1]) / dx,
        //         _ => (v.cells.A[i + 1] - v.cells.A[i - 1]) / (2.0 * dx),
        //     };
        //
        //     v.cells.F[i] = A_derivative * ((v.consts.cs2) - (H_derivative));
        // }
        // v.cells.F = zip_map_enumerate!(v.cells, |(i, A0, beta)| {
        //     let P_derivative = (beta) / (2.0 * (v.cells.A[i] * A0).sqrt());
        //     let H_derivative = (v.cells.A[i] / v.consts.rho) * P_derivative;
        //
        //     let A_derivative = match i {
        //         0 => (v.cells.A[i + 1] - v.cells.A[i]) / v.consts.dx,
        //         i if (i == v.x_last) => (v.cells.A[i] - v.cells.A[i - 1]) / v.consts.dx,
        //         _ => (v.cells.A[i + 1] - v.cells.A[i - 1]) / (2.0 * v.consts.dx),
        //     };
        //
        //     A_derivative * ((v.consts.cs2) - (H_derivative))
        // });

        for i in 0..v.x_dim {
            v.cells.u[i] += ((v.consts.dt / 2.0) * v.cells.F[i]) / v.cells.A[i];
            let uc = v.cells.u[i] / v.consts.c;
            let uc2 = uc * uc;
            v.cells.f0[i] = (1.0 / 3.0) * v.cells.A[i] * (2.0 - 3.0 * uc2);
            v.cells.f1[i] = (1.0 / 6.0) * v.cells.A[i] * (1.0 - 3.0 * uc + 3.0 * uc2);
            v.cells.f2[i] = (1.0 / 6.0) * v.cells.A[i] * (1.0 + 3.0 * uc + 3.0 * uc2);
            v.cells.A[i] = v.cells.f0[i] + v.cells.f1[i] + v.cells.f2[i];
        }

        // v.cells.u = zip_map!(v.cells, |(A, u, F)| u + ((v.consts.dt / 2.0) * F) / A);
        // //
        // v.cells.f = zip_map!(v.cells, |(A, u)| {
        //     let uc = u / v.consts.c;
        //     let uc2 = uc * uc;
        //     (
        //         (1.0 / 3.0) * A * (2.0 - 3.0 * uc2),
        //         (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2),
        //         (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2),
        //     )
        // });
        //
        // v.cells.A = zip_map!(v.cells, |f| f.0 + f.1 + f.2);
    });
}

pub fn functional2<T: Compute>(sim: &mut Simulation<T>) {
    let dt = sim.consts.dt;
    let owo = sim.consts.c;
    let rho = sim.consts.rho;
    let dx = sim.consts.dx;
    let cs2 = sim.consts.cs2;

    for v in &mut sim.vessels {
        for i in 0..v.x_dim {
            let P_derivative = (v.cells.beta[i]) / (2.0 * (v.cells.A[i] * v.cells.A0[i]).sqrt());
            let H_derivative = (v.cells.A[i] / v.consts.rho) * P_derivative;

            let A_derivative = match i {
                0 => (v.cells.A[i + 1] - v.cells.A[i]) / v.consts.dx,
                i if (i == v.x_last) => (v.cells.A[i] - v.cells.A[i - 1]) / dx,
                _ => (v.cells.A[i + 1] - v.cells.A[i - 1]) / (2.0 * dx),
            };

            v.cells.F[i] = A_derivative * ((v.consts.cs2) - (H_derivative));
        }
        // for i in 1..v.x_last {
        //     let P_derivative = (v.cells.beta[i]) / (2.0 * (v.cells.A[i] * v.cells.A0[i]).sqrt());
        //     let H_derivative = (v.cells.A[i] / v.consts.rho) * P_derivative;
        //
        //     let A_derivative = (v.cells.A[i + 1] - v.cells.A[i - 1]) / (2.0 * v.consts.dx);
        //
        //     v.cells.F[i] = A_derivative * ((v.consts.cs2) - (H_derivative));
        // }
        // let i = 0;
        // let P_derivative = (v.cells.beta[i]) / (2.0 * (v.cells.A[i] * v.cells.A0[i]).sqrt());
        // let H_derivative = (v.cells.A[i] / v.consts.rho) * P_derivative;
        // let A_derivative = (v.cells.A[i + 1] - v.cells.A[i]) / v.consts.dx;
        // v.cells.F[i] = A_derivative * ((v.consts.cs2) - (H_derivative));
        // let i = v.x_last;
        // let P_derivative = (v.cells.beta[i]) / (2.0 * (v.cells.A[i] * v.cells.A0[i]).sqrt());
        // let H_derivative = (v.cells.A[i] / v.consts.rho) * P_derivative;
        // let A_derivative = (v.cells.A[i] - v.cells.A[i - 1]) / v.consts.dx;
        // v.cells.F[i] = A_derivative * ((v.consts.cs2) - (H_derivative));
        // v.cells.F = sim.compute.forcing_term(v);
        // v.cells.u = sim.compute.velocity_init(v);
        zip_for_each!(v.cells, |(A, mut u, F)| {
            *u = *u + ((dt / 2.0) * F) / A;
        });

        zip_for_each!(v.cells, |(A, u, mut f)| {
            let uc = u / owo;
            let uc2 = uc * uc;

            f.0 = (1.0 / 3.0) * A * (2.0 - 3.0 * uc2);
            f.1 = (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2);
            f.2 = (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2);
        });

        zip_for_each!(v.cells, |(mut A, f)| {
            *A = f.0 + f.1 + f.2;
        });
        //
        // izip2!(&v.cells.A, &mut v.cells.u, &v.cells.F).for_each(|(A, u, F)| {
        //     *u = *u + ((dt / 2.0) * *F) / *A;
        // });
        //
        // izip2!(&v.cells.A, &v.cells.u, &mut v.cells.f).for_each(|(A, u, f)| {
        //     let uc = *u / owo;
        //     let uc2 = uc * uc;
        //
        //     f.0 = (1.0 / 3.0) * A * (2.0 - 3.0 * uc2);
        //     f.1 = (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2);
        //     f.2 = (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2);
        // });
        //
        // izip2!(&mut v.cells.A, &v.cells.f).for_each(|(A, f)| {
        //     *A = f.0 + f.1 + f.2;
        // });

        // zip_for_each!(v.cells, |(A, u, F)| {
        //     *u = *u + ((dt / 2.0) * F) / A;
        // });
        // zip_for_each!(v.cells, |(f, A, u)| {
        //     let uc = u / owo;
        //     let uc2 = uc * uc;
        //
        //     f.0 = (1.0 / 3.0) * A * (2.0 - 3.0 * uc2);
        //     f.1 = (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2);
        //     f.2 = (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2);
        // });
        // v.cells.f0 = zip_map!(v.cells, |(A, u)| {
        //     let uc = u / v.consts.c;
        //     let uc2 = uc * uc;
        //     (1.0 / 3.0) * A * (2.0 - 3.0 * uc2)
        // })
        // .collect();
        // v.cells.f1 = zip_map!(v.cells, |(A, u)| {
        //     let uc = u / v.consts.c;
        //     let uc2 = uc * uc;
        //     (1.0 / 3.0) * A * (2.0 - 3.0 * uc2)
        // })
        // .collect();
        // v.cells.f2 = zip_map!(v.cells, |(A, u)| {
        //     let uc = u / v.consts.c;
        //     let uc2 = uc * uc;
        //     (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2)
        // })
        // .collect();
        // v.cells.f = sim.compute.FEQ(v);
        // v.cells.A = zip_map!(v.cells, |fff| fff.0 + fff.1 + fff.2); //sim.compute.area(v);
    }
}
// fn iterative(sim: &mut Simulation) {
//     for v in &mut sim.vessels {
//         for i in 0..v.x_dim {
//             v.compute_F(i); // In-place
//         }
//
//         for i in 0..v.x_dim {
//             v.populations_init(i);
//             // v.cells.u[i] += ((v.consts.dt / 2.0) * v.cells.F[i]) / v.cells.A[i];
//             // let uc = v.cells.u[i] / v.consts.c;
//             // let uc2 = uc * uc;
//             // v.cells.f0[i] = (1.0 / 3.0) * v.cells.A[i] * (2.0 - 3.0 * uc2);
//             // v.cells.f1[i] = (1.0 / 6.0) * v.cells.A[i] * (1.0 - 3.0 * uc + 3.0 * uc2);
//             // v.cells.f2[i] = (1.0 / 6.0) * v.cells.A[i] * (1.0 + 3.0 * uc + 3.0 * uc2);
//             // v.populations_init(i);
//             // v.cells.f[i] = v.computeEQ(i);
//             //
//             // v.cells.f0[i] = feq.0;
//             // v.cells.f1[i] = feq.1;
//             // v.cells.f2[i] = feq.2;
//             // v.cells.A[i] = v.cells.f[i].0 + v.cells.f[i].1 + v.cells.f[i].2;
//         }
//     }
// }
//
//

pub fn rotate(v: &mut Vec<f64>) {
    v.rotate_left(1);
    *v.last_mut().unwrap() = 10.0;
}

pub unsafe fn copy(v: &mut Vec<f64>) {
    ptr::copy(v.as_ptr().add(1), v.as_mut_ptr(), v.len() - 1);
    *v.last_mut().unwrap() = 10.0;
}

pub fn remove_push(v: &mut Vec<f64>) {
    v.remove(0);
    v.push(10.0);
}

pub fn vec_decque(v: &mut VecDeque<f64>) {
    v.pop_front();
    v.push_back(10.0);
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Vector shifting");

    // let mut parse = SimulationParsing::read_json(&"vascularNetworks/adan56.json".to_string());

    group.bench_function("Functional VecDeque", |b| {
        let mut sim = Simulation::<AlgoBase>::new("vascularNetworks/adan56.json");
        b.iter(|| black_box(functional(&mut sim)))
    });

    group.bench_function("Functional Vec ", |b| {
        let mut sim = Simulation::<AlgoBase>::new("vascularNetworks/adan56.json");
        b.iter(|| black_box(functional2(&mut sim)))
    });

    group.bench_function("Functional VecDeque2", |b| {
        let mut sim = Simulation::<AlgoBase>::new("vascularNetworks/adan56.json");
        b.iter(|| black_box(functional(&mut sim)))
    });

    group.bench_function("Functional Vec2 ", |b| {
        let mut sim = Simulation::<AlgoBase>::new("vascularNetworks/adan56.json");
        b.iter(|| black_box(functional2(&mut sim)))
    });
    // let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    // group.plot_config(plot_config);
    //
    // parse.dx = 0.00000001;
    // let consts = Constants::new(&parse);
    // let v = Vessel::new(&parse.vessels[0], consts.clone());
    // println!("owo {}", v.cells.A.len());
    //
    // for i in [10u32, 100u32, 1000u32, 10000u32, 100000u32, 1000000u32].iter() {
    //     // base.pow(*i);
    //
    //     group.bench_with_input(BenchmarkId::new("Vec::rotate", i), &i, |b, i| unsafe {
    //         let mut v: Vec<f64> = v.cells.A.clone();
    //         v.drain((**i as usize)..);
    //         b.iter(|| black_box(rotate(&mut v)))
    //     });
    //
    //     group.bench_with_input(BenchmarkId::new("ptr::copy", i), &i, |b, i| unsafe {
    //         let mut v: Vec<f64> = v.cells.A.clone();
    //         v.drain((**i as usize)..);
    //         b.iter(|| black_box(copy(&mut v)))
    //     });
    //     group.bench_with_input(BenchmarkId::new("Vec::remove_push", i), &i, |b, i| unsafe {
    //         let mut v: Vec<f64> = v.cells.A.clone();
    //         v.drain((**i as usize)..);
    //         b.iter(|| black_box(remove_push(&mut v)))
    //     });
    //
    //     group.bench_with_input(BenchmarkId::new("VecDecque::pop_push", i), &i, |b, i| unsafe {
    //         let mut v: Vec<f64> = v.cells.A.clone();
    //         v.drain((**i as usize)..);
    //         let mut v = VecDeque::from(v);
    //         b.iter(|| black_box(vec_decque(&mut v)))
    //     });

    // group.bench_with_input("ptr::copy", |b| unsafe {
    //     let mut v = v.cells.A.clone();
    //     b.iter(|| black_box(memmove(&mut v)))
    // });
    //
    // group.bench_with_input("Vec::remove_push", |b| unsafe {
    //     let mut v = v.cells.A.clone();
    //     b.iter(|| black_box(remove_push(&mut v)))
    // });
    //
    // group.bench_with_input("VecDecque::pop_push", |b| unsafe {
    //     b.iter(|| black_box(vec_decque(&mut v)))
    // });
    // }

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

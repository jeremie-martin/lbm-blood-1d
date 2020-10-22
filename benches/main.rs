// https://bheisler.github.io/criterion.rs/book/getting_started.html

#![allow(non_snake_case)]

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use lbm_blood_1d::lbm_algorithm::*;
use lbm_blood_1d::simulation::*;

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

pub fn functional<T: Algorithm>(sim: &mut Simulation<T>) {
    sim.vessels.iter_mut().for_each(|v| {
        v.cells.F = zip_map_enumerate!(v.cells, |(i, A0, beta)| {
            let P_derivative = (beta) / (2.0 * (v.cells.A[i] * A0).sqrt());
            let H_derivative = (v.cells.A[i] / v.consts.rho) * P_derivative;

            let A_derivative = match i {
                0 => (v.cells.A[i + 1] - v.cells.A[i]) / v.consts.dx,
                i if (i == v.x_last) => (v.cells.A[i] - v.cells.A[i - 1]) / v.consts.dx,
                _ => (v.cells.A[i + 1] - v.cells.A[i - 1]) / (2.0 * v.consts.dx),
            };

            A_derivative * ((v.consts.cs2) - (H_derivative))
        });

        v.cells.u = zip_map!(v.cells, |(A, u, F)| u + ((v.consts.dt / 2.0) * F) / A);

        v.cells.f = zip_map!(v.cells, |(A, u)| {
            let uc = u / v.consts.c;
            let uc2 = uc * uc;
            (
                (1.0 / 3.0) * A * (2.0 - 3.0 * uc2),
                (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2),
                (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2),
            )
        });

        v.cells.A = zip_map!(v.cells, |f| f.0 + f.1 + f.2);
    });
}

pub fn functional2<T: Algorithm>(sim: &mut Simulation<T>) {
    for v in &mut sim.vessels {
        v.cells.F = sim.algo.compute_forcing_term(v);
        v.cells.u = sim.algo.compute_velocity(v);
        v.cells.f = sim.algo.compute_FEQ(v);
        v.cells.A = sim.algo.compute_area(v);
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

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Vessels initialization");

    group.bench_function("Functional", |b| {
        let mut sim = Simulation::<AlgoBase>::new("vascularNetworks/adan56.json");
        b.iter(|| black_box(functional(&mut sim)))
    });
    group.bench_function("Function from lib", |b| {
        let mut sim = Simulation::<AlgoBase>::new("vascularNetworks/adan56.json");
        b.iter(|| black_box(functional2(&mut sim)))
    });

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

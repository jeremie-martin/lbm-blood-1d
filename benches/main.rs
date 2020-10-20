// https://bheisler.github.io/criterion.rs/book/getting_started.html

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lbm_blood_1d::simulation::*;

pub fn functional(sim: &Simulation) {
    let mut vessels = sim.vessels.clone();

    vessels.iter_mut().for_each(|v| {
        v.cells.F = (0..v.x_last).map(|i| v.compute_F_ret(sim.cs2, i)).collect();

        (0..v.x_last).for_each(|i| {
            v.cells.u[i] += ((sim.dt / 2.0) * v.cells.F[i]) / v.cells.A[i];
            let feq = sim.computeFEQ(v.cells.A[i], v.cells.u[i]);

            v.cells.f0[i] = feq.0;
            v.cells.f1[i] = feq.1;
            v.cells.f2[i] = feq.2;
            v.cells.A[i] = v.cells.f0[i] + v.cells.f1[i] + v.cells.f2[i];
        })
    });

    // for v in &mut vessels {
    //     for i in 0..v.x_last {
    //         v.compute_F(sim.cs2, i);
    //     }
    //
    //     for i in 0..v.x_last {
    //         v.cells.u[i] += ((sim.dt / 2.0) * v.cells.F[i]) / v.cells.A[i];
    //         let feq = sim.computeFEQ(v.cells.A[i], v.cells.u[i]);
    //
    //         v.cells.f0[i] = feq.0;
    //         v.cells.f1[i] = feq.1;
    //         v.cells.f2[i] = feq.2;
    //         v.cells.A[i] = v.cells.f0[i] + v.cells.f1[i] + v.cells.f2[i];
    //     }
    // }
}

fn iterative(sim: &Simulation) {
    let mut vessels = sim.vessels.clone();

    for v in &mut vessels {
        for i in 0..v.x_last {
            v.compute_F(sim.cs2, i);
        }

        for i in 0..v.x_last {
            v.cells.u[i] += ((sim.dt / 2.0) * v.cells.F[i]) / v.cells.A[i];
            let feq = sim.computeFEQ(v.cells.A[i], v.cells.u[i]);

            v.cells.f0[i] = feq.0;
            v.cells.f1[i] = feq.1;
            v.cells.f2[i] = feq.2;
            v.cells.A[i] = v.cells.f0[i] + v.cells.f1[i] + v.cells.f2[i];
        }
    }
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut sim = Simulation::new("vascularNetworks/bifur.json");
    let mut sim2 = Simulation::new("vascularNetworks/bifur.json");

    c.bench_function("Iterative", |b| b.iter(|| black_box(iterative(&sim2))));
    c.bench_function("Functional", |b| b.iter(|| black_box(functional(&sim))));
    c.bench_function("Iterative2", |b| b.iter(|| black_box(iterative(&sim2))));
    c.bench_function("Functional2", |b| b.iter(|| black_box(functional(&sim))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

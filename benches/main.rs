// https://bheisler.github.io/criterion.rs/book/getting_started.html

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lbm_blood_1d::simulation::*;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Marshalling single bifurcation", |b| {
        b.iter(|| black_box(read_json("vascularNetworks/bifur.json")))
    });
    c.bench_function("Marshalling adan56 vascular network", |b| {
        b.iter(|| black_box(read_json("vascularNetworks/adan56.json")))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

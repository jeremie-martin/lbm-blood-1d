//! This module

use crate::simulation_parsing::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Constants {
    /// Space step [m]
    #[serde(default)]
    pub dx: f64,
    /// Time step, depends on dx [m]
    pub dt: f64,
    /// Blood dynamic viscosity [Pa.s = kg.m-1.s-1]
    pub mu: f64,
    /// Blood density [kg.m-3]
    pub rho: f64,
    /// Lattice velocity [m.s-1]
    pub c: f64,
    /// Lattice velocity squared [m2.s-2]
    pub c2: f64,
    /// Lattice speed of sound [m.s-1]
    pub cs: f64,
    /// Lattice speed of sound squared [m2.s-1]
    pub cs2: f64,
    /// Kinematic viscosity (ratio of a fluid's dynamic viscosity to the fluid's density) [m2.s-1]
    pub nu: f64,
    /// Relaxation time [s]
    pub tau: f64,
    /// Adjusted Relaxation time [s]
    pub tau_bar: f64,
    /// Relaxation rate [dimensionless]
    pub omega: f64,
    /// Relaxation rate of the Force [dimensionless]
    pub omega_force: f64,
    pub Cl: f64,
    pub Ct: f64,
    pub CA: f64,
    pub Cu: f64,
    pub CQ: f64,
    pub CF: f64,
    pub cs2_lat: f64,

    pub l0_p: f64,
    pub u0_p: f64,
    pub t0_p: f64,
    pub nu_p: f64,

    pub l_d: f64,
    pub u_d: f64,
    pub t_d: f64,
    pub nu_d: f64,
    pub dx_d: f64,
    pub dt_d: f64,

    pub u0_lb: f64,
    pub nu_lb: f64,

    pub rho_p2d: f64,
    pub l_p2d: f64,
    pub t_p2d: f64,
    pub u_p2d: f64,

    pub rho_d2lb: f64,
    pub l_d2lb: f64,
    pub t_d2lb: f64,
    pub u_d2lb: f64,
}

impl Constants {
    pub fn new(parse: &SimulationParsing, dx: f64, fact: f64) -> Constants {
        let mu = parse.mu;
        let rho = parse.rho;

        // Kinematic viscosity (ratio of a fluid's dynamic viscosity to the fluid's density) [m2.s-1]

        // let beta = 85166.83575965902;
        // let beta = 88738.13415206889;
        // let beta = 94070.22632255145;
        // let beta = 74064.06939080187;
        let beta = 80803.92925723408;
        // let beta = 82661.79266214612;
        // let beta = 79340.36394324744;
        // let beta = 78118.48451568425;
        let l0_p = 1.0; // m
        let u0_p = (beta / (2.0 * rho)); // m.s
        let t0_p = l0_p / u0_p; // s
        let nu_p = mu / rho;

        let Re = u0_p * l0_p / nu_p;
        println!("Re: {} {} {}", Re, u0_p, t0_p);

        let dx_p = dx;
        let dt_p = dx / (u0_p);

        let t_d = 1.0;
        let l_d = 1.0;
        let u_d = u0_p * t0_p;
        let nu_d = 1.0 / Re;

        let l_p2d = 1.0 / l0_p;
        let t_p2d = 1.0 / t0_p;
        let u_p2d = t0_p / l0_p;
        let rho_p2d = 1.0 / rho;

        let dx_d = dx_p * l_p2d;
        let dt_d = dt_p * t_p2d;
        let dx_d = dx;
        let dt_d = dx / (u0_p);
        // let dt_d = dx / (u0_p * 3.0f64.sqrt());

        let l_d2lb = 1.0 / dx_d;
        let t_d2lb = 1.0 / dt_d;
        let u_d2lb = dt_d / dx_d;
        let rho_d2lb = 1.0;

        let u0_lb = dt_d / dx_d;
        let nu_lb = (dt_d / (dx_d * dx_d)) * nu_d;

        let nu = mu / rho;
        let Cl = 1.0 / (l_p2d * l_d2lb);
        let Ct = 1.0 / (t0_p * t_p2d * t_d2lb);
        let Cu = 1.0 / (u0_p * u_p2d * u_d2lb);
        let CA = 0.00048780487804878054;
        let CQ = CA * Cu;
        let CF = CA * dx_d / (dt_d * dt_d);
        // let CF = CA * Cl / (Ct * Ct);
        // let cs2_lat = (Cl / Ct) * (Cl / Ct) / 3.0;
        let cs2_lat = (dx_d / dt_d).powi(2) / 3.0;

        println!(
            "{} {} {} {} {} {} {} {} {} {} {}",
            dx_d, Cl, dt_d, Ct, t0_p, t_d, Cu, u_d2lb, CA, u0_p, nu_lb
        );

        let dx = 1.0;
        let dt = 1.0;

        let c = dx / dt;
        let c2 = c * c;
        let cs = c / 3.0f64.sqrt();
        let cs2 = c2 / 3.0;

        // Relaxation time [s]
        // let tau = nu / cs2;
        let tau = (nu / cs2) * (Cl / (Ct * Ct));
        // let tau_bar = 0.51;

        // Adjusted tau
        let tau_bar = tau + (dt / 2.0);
        let tau_bar = (nu_lb / cs2) + 0.5;

        let omega = dt / tau_bar;
        let omega_force = dt * (1.0 - (dt / (2.0 * tau_bar)));

        Constants {
            dx,
            dt,
            mu,
            rho,
            c,
            cs,
            c2,
            cs2,
            nu,
            tau,
            tau_bar,
            omega,
            omega_force,
            Cl,
            Ct,
            CA,
            Cu,
            CQ,
            CF,
            cs2_lat,
            l0_p,
            nu_p,
            t0_p,
            u0_p,
            l_d,
            t_d,
            nu_d,
            u_d,
            dx_d,
            dt_d,
            u0_lb,
            nu_lb,
            l_p2d,
            t_p2d,
            u_p2d,
            rho_p2d,
            l_d2lb,
            t_d2lb,
            u_d2lb,
            rho_d2lb,
        }
    }
}

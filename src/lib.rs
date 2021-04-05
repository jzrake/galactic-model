//! 10-parameter axisymmetric model of a disk galaxy, using
//! the Plummer model for the central bulge, two Miyamoto-
//! Nagai models for the thin and thick disks, and a logarithmic
//! model for the dark matter halo.

use std::f64::consts::PI;
use derive_more::{Add,Sub, Mul, Div};

// Define struct for galactic model parameters,
// including the gravitational constant. Here, slr
// stands for solar masses.
#[derive(Debug, Add, Sub, Mul, Div)]
struct GalacticModel {
    g:   f64, // gravitational constant (kpc*kpc*kpc/Myr/Myr/slr)
    m_b: f64, // mass of central bulge (slr)
    a_b: f64, // radial scale length of central bulge (kpc)
    v_h: f64, // radial velocities at large distances (kpc/Myr)
    a_h: f64, // radial scale length of dark matter halo (kpc)
    m_s: f64, // mass of thin disk (slr)
    a_s: f64, // radial scale length of thin disk (kpc)
    b_s: f64, // vertical scale length of thin disk (kpc)
    m_g: f64, // mass of thick disk (slr)
    a_g: f64, // radial scale length of thick disk (kpc)
    b_g: f64, // vertical scale length of thick disk (kpc)
}

impl GalacticModel {

    // Gravitational potential
    fn potential(&self, r: f64, z: f64) -> f64 {

        let phi_bulge      = -self.g * self.m_b / (r*r + z*z + self.a_b*self.a_b).sqrt();
        let phi_thin_disk  = -self.g * self.m_s / (r*r + (self.a_s + (z*z + self.b_s*self.b_s).sqrt()).powi(2)).sqrt();
        let phi_thick_disk = -self.g * self.m_g / (r*r + (self.a_g + (z*z + self.b_g*self.b_g).sqrt()).powi(2)).sqrt();
        let phi_halo       = 0.5*(self.v_h).powi(2)*((r*r + z*z + self.a_h*self.a_h).ln());

        return phi_bulge + phi_thin_disk + phi_thick_disk + phi_halo

    }

    // Mass density profile obtained via Poisson's equation using the above potential
    fn density(&self, r: f64, z: f64) -> f64 {

        let r_z_a            = r*r + z*z + self.a_b*self.a_b;
        let rho_bulge_z      = (1.0/(4.0*PI))*self.m_b * ((r_z_a).sqrt() - z*z*(r_z_a).powf(-1.0/2.0)) / (r_z_a);
        let rho_bulge_r      = (1.0/(4.0*PI))*self.m_b * (2.0*(r_z_a).sqrt() - r*r*(r_z_a).powf(-1.0/2.0)) / (r_z_a);
        let rho_bulge        = rho_bulge_z + rho_bulge_r;

        let z_b_s            = (z*z + self.b_s*self.b_s).sqrt();
        let z_b_a_s          = self.a_s + z_b_s;
        let r_z_b_a_s        = r*r + (z_b_a_s).powi(2);
        let rho_thin_disk_z  = self.m_s/4.0/PI * (z_b_s*r_z_b_a_s.powf(3.0/2.0) + z*z*r_z_b_a_s.powf(3.0/2.0)
                               - z*z*z_b_a_s*r_z_b_a_s.powf(3.0/2.0)/z_b_s - 3.0*z*z*z_b_a_s*z_b_a_s*r_z_b_a_s.sqrt())
                               / (z_b_s * r_z_b_a_s.powi(3));
        let rho_thin_disk_r  = self.m_s/4.0/PI * (2.0 * r_z_b_a_s.sqrt() - r*r*r_z_b_a_s.powf(-1.0/2.0)) / r_z_b_a_s;
        let rho_thin_disk    = rho_thin_disk_z + rho_thin_disk_r;

        let z_b_g            = (z*z + self.b_g*self.b_g).sqrt();
        let z_b_a_g          = self.a_g + z_b_g;
        let r_z_b_a_g        = r*r + (z_b_a_g).powi(2);
        let rho_thick_disk_z = self.m_g/4.0/PI * (z_b_g*r_z_b_a_g.powf(3.0/2.0) + z*z*r_z_b_a_g.powf(3.0/2.0)
                               - z*z*z_b_a_g*r_z_b_a_g.powf(3.0/2.0)/z_b_g - 3.0*z*z*z_b_a_g*z_b_a_g*r_z_b_a_g.sqrt())
                               / (z_b_g * r_z_b_a_g.powi(3));
        let rho_thick_disk_r = self.m_g/4.0/PI * (2.0 * r_z_b_a_g.sqrt() - r*r*r_z_b_a_g.powf(-1.0/2.0)) / r_z_b_a_g;
        let rho_thick_disk   = rho_thick_disk_z + rho_thick_disk_r;

        let rho_halo_z       = self.v_h*self.v_h*(r*r - z*z + self.a_h*self.a_h)/(r*r + z*z + self.a_h*self.a_h).powi(2);
        let rho_halo_r       = 2.0*self.v_h*self.v_h*(z*z + self.a_h*self.a_h)/(r*r + z*z + self.a_h*self.a_h).powi(2);
        let rho_halo         = rho_halo_z + rho_halo_r;

        return rho_bulge + rho_thin_disk + rho_thick_disk + rho_halo

    }

    // The z-component of the gravitational field obtained via the negative gradient of the above potential
    fn g_field_z(&self, r: f64, z: f64) -> f64 {

        let gfz_bulge      = -self.g*self.m_b/(r*r + z*z + self.a_b*self.a_b).powf(3.0/2.0);
        let gfz_thin_disk  = -self.g*self.m_b*z*(self.a_s + (z*z + self.b_s*self.b_s).sqrt())
                             / (r*r + (self.a_s + (z*z + self.b_s*self.b_s).sqrt()).powi(2)).powf(3.0/2.0);
        let gfz_thick_disk = -self.g*self.m_b*z*(self.a_g + (z*z + self.b_g*self.b_g).sqrt())
                             / (r*r + (self.a_g + (z*z + self.b_g*self.b_g).sqrt()).powi(2)).powf(3.0/2.0);
        let gfz_halo       = -self.v_h*self.v_h*z/(r*r + z*z + self.a_h*self.a_h);
        return gfz_bulge + gfz_thin_disk + gfz_thick_disk + gfz_halo

    }

}
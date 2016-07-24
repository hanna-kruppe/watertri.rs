//! An implementation of the ray-triangle intersection algorithm described in:
//!
//! > Sven Woop, Carsten Benthin, and Ingo Wald. "Watertight ray/triangle intersection."
//! > Journal of Computer Graphics Techniques (JCGT) 2.1 (2013): 65-82.
//!
//! Does not perform backface culling.
#![allow(non_snake_case)]
// Variable names from the paper (appendix A) are not snake_case
extern crate cgmath;

use cgmath::Vector3;

/// Precomputed data depending only on the ray.
#[derive(Clone, Debug)]
pub struct RayData {
    kx: usize,
    ky: usize,
    kz: usize,
    Sx: f32,
    Sy: f32,
    Sz: f32,
    org: Vector3<f32>,
}

impl RayData {
    /// Pre-compute the transformation that is applied to all triangles.
    pub fn new(org: Vector3<f32>, dir: Vector3<f32>) -> RayData {
        // The paper swaps kx and ky if dir[kz] is negative, to preserve winding order.
        // But winding order is only relevant for backface culling, which we don't perform.
        let kz = max_dim(dir);
        let kx = (kz + 1) % 3;
        let ky = (kz + 2) % 3;
        RayData {
            kx: kx,
            ky: ky,
            kz: kz,
            Sx: dir[kx] / dir[kz],
            Sy: dir[ky] / dir[kz],
            Sz: 1.0 / dir[kz],
            org: org,
        }
    }

    /// Perform the intersection calculation.
    pub fn intersect(&self, A: Vector3<f32>, B: Vector3<f32>, C: Vector3<f32>) -> Option<Intersection> {
        let (Sx, Sy, Sz, org) = (self.Sx, self.Sy, self.Sz, self.org);
        let (kx, ky, kz) = (self.kx, self.ky, self.kz);
        let (A, B, C) = (A - org, B - org, C - org);
        let Ax = A[kx] - Sx * A[kz];
        let Ay = A[ky] - Sy * A[kz];
        let Bx = B[kx] - Sx * B[kz];
        let By = B[ky] - Sy * B[kz];
        let Cx = C[kx] - Sx * C[kz];
        let Cy = C[ky] - Sy * C[kz];

        let mut U = Cx * By - Cy * Bx;
        let mut V = Ax * Cy - Ay * Cx;
        let mut W = Bx * Ay - By * Ax;

        if U == 0. || V == 0. || W == 0. {
            let CxBy = (Cx as f64) * (By as f64);
            let CyBx = (Cy as f64) * (Bx as f64);
            U = (CxBy - CyBx) as f32;
            let AxCy = (Ax as f64) * (Cy as f64);
            let AyCx = (Ay as f64) * (Cx as f64);
            V = (AxCy - AyCx) as f32;
            let BxAy = (Bx as f64) * (Ay as f64);
            let ByAx = (By as f64) * (Ax as f64);
            W = (BxAy - ByAx) as f32;
        }

        if (U < 0. || V < 0. || W < 0.) && (U > 0. || V > 0. || W > 0.) {
            return None;
        }

        let det = U + V + W;
        if det == 0. {
            return None;
        }

        let Az = Sz * A[kz];
        let Bz = Sz * B[kz];
        let Cz = Sz * C[kz];
        let T = U * Az + V * Bz + W * Cz;

        let rcpDet = 1.0 / det;
        Some(Intersection {
            t: T * rcpDet,
            u: U * rcpDet,
            v: V * rcpDet,
            w: W * rcpDet,
        })
    }
}

/// Geometric information about a ray-triangle intersection.
pub struct Intersection {
    /// Parametric distance from the ray origin to the intersection.
    pub t: f32,
    /// Barycentric coordinate.
    pub u: f32,
    /// Barycentric coordinate.
    pub v: f32,
    /// Barycentric coordinate.
    pub w: f32,
}

fn max_dim(v: Vector3<f32>) -> usize {
    let (x, y, z) = (v.x.abs(), v.y.abs(), v.z.abs());
    if x > y {
        // y isn't the maximum, so it's either x or z
        if x > z {
            0
        } else {
            2
        }
    } else if y > z {
        // x isn't the maximum, so it's either y or z
        1
    } else {
        2
    }
}

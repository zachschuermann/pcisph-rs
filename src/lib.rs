extern crate glam;
use glam::{vec2, vec3, Vec2, Vec3};
use kdbush::KDBush;
use rayon::prelude::*;
use std::f32::consts::PI;

#[derive(Clone, Copy, Default)]
pub struct Particle {
    pub x: Vec2,
    pub xlast: Vec2,
    pub v: Vec2,
    pub m: f32,
    pub p: f32,
    pub pv: f32,
    pub d: f32,
    pub dv: f32,
}

impl Particle {
    pub fn new() -> Self {
        Self {
            x: Vec2::ZERO,
            v: Vec2::ZERO,
            m: 1.0,
            ..Default::default()
        }
    }

    pub fn new_with_pos(x: f32, y: f32) -> Self {
        Self {
            x: Vec2::new(x, y),
            v: Vec2::ZERO,
            m: 1.0,
            ..Default::default()
        }
    }
}

const SOLVER_STEPS: usize = 10;
const REST_DENS: f32 = 45.0;
const STIFFNESS: f32 = 0.08;
const STIFF_APPROX: f32 = 0.1; // TODO: better naming
const SURFACE_TENSION: f32 = 0.0001;
const LINEAR_VISC: f32 = 0.25;
const QUAD_VISC: f32 = 0.5;
const PARTICLE_RADIUS: f32 = 0.03;
const H: f32 = 6.0 * PARTICLE_RADIUS;
const H2: f32 = H * H;
const DT: f32 = 0.000003; //(1.0 / 60.0) / SOLVER_STEPS as f32;
const DT2: f32 = DT * DT;
const KERN: f32 = 20.0 / (2.0 * PI * H * H);
const KERN_NORM: f32 = 30.0 / (2.0 * PI * H * H);
const EPS: f32 = 1e-6;
const EPS2: f32 = EPS * EPS;
pub const G: Vec2 = glam::const_vec2!([0.0, -9.8]);
pub const WINDOW_WIDTH: u32 = 1280;
pub const WINDOW_HEIGHT: u32 = 1024;
pub const VIEW_WIDTH: f32 = 20.0;
pub const VIEW_HEIGHT: f32 = WINDOW_HEIGHT as f32 * VIEW_WIDTH / WINDOW_WIDTH as f32;

pub struct State {
    pub particles: Vec<Particle>,
    boundaries: Vec<Vec3>,
    spatial_index: KDBush,
}

impl State {
    pub fn new() -> Self {
        Self {
            particles: Vec::new(),
            boundaries: Vec::new(),
            spatial_index: KDBush::new(0, kdbush::DEFAULT_NODE_SIZE),
        }
    }

    pub fn init_dam_break(&mut self, dam_max_particles: usize) {
        self.boundaries.clear();
        self.boundaries.push(vec3(1.0, 0.0, 0.0)); // left
        self.boundaries.push(vec3(0.0, 1.0, 0.0)); // bottom
        self.boundaries.push(vec3(-1.0, 0.0, -VIEW_WIDTH)); // right
        self.boundaries.push(vec3(0.0, -1.0, -VIEW_HEIGHT)); // top

        let mut start = vec2(0.25 * VIEW_WIDTH, 0.95 * VIEW_HEIGHT);
        let x0 = start.x;
        let num = f32::sqrt(dam_max_particles as f32) as usize;
        for i in 0..num {
            for j in 0..num {
                self.particles
                    .push(Particle::new_with_pos(start.x, start.y));
                start.x += 2.0 * PARTICLE_RADIUS + PARTICLE_RADIUS;
            }
            start.x = x0;
            start.y -= 2.0 * PARTICLE_RADIUS + PARTICLE_RADIUS;
        }
        println!(
            "Initialized dam break with {} particles",
            self.particles.len()
        );
    }

    fn integrate(&mut self) {
        self.particles.par_iter_mut().for_each(|p| {
            p.v += G * DT;
            p.xlast = p.x;
            p.x += DT * p.v;
        });
    }

    fn grid_insert(&mut self) {
        /*
        self.particles
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, p)| {
                let xind = (p.x.x / CELL_SIZE).floor() as usize;
                let yind = (p.x.y / CELL_SIZE).floor() as usize;
                let xind = usize::max(1, usize::min(GRID_WIDTH - 2, xind));
                let yind = usize::max(1, usize::min(GRID_HEIGHT - 2, yind));
                //self.grid[xind + yind * GRID_WIDTH].append(p);
                //self.grid_indices[i] = Vec2::new(xind as f32, yind as f32);
            });
             */
    }

    fn compute_forces(&mut self) {
        let particles_initial = self.particles.clone();
        self.particles
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, pi)| {
                let mut dens = 0.0;
                let mut dens_proj = 0.0;
                for (j, pj) in particles_initial.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    let dx = pj.x - pi.x;
                    let r2 = dx.length_squared();
                    if r2 < EPS2 || r2 > H2 {
                        continue;
                    }
                    let r = f32::sqrt(r2);
                    let a = 1.0 - r / H;
                    dens += pj.m * a * a * a * KERN;
                    dens_proj += pj.m * a * a * a * a * KERN_NORM;
                }
                pi.d = dens;
                pi.dv = dens_proj;
                pi.p = STIFFNESS * (dens - pi.m * REST_DENS);
                pi.pv = STIFF_APPROX * dens_proj;
            })
    }

    fn project_correct(&mut self) {
        let particles_initial = self.particles.clone();
        let bounds = self.boundaries.clone();
        self.particles
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, pi)| {
                // project
                let mut xproj = pi.x.clone();
                for (j, pj) in particles_initial.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    let dx = pj.x - pi.x;
                    let r2 = dx.length_squared();
                    if r2 < EPS2 || r2 > H2 {
                        continue;
                    }
                    let r = f32::sqrt(r2);
                    let a = 1.0 - r / H;
                    let d = DT2
                        * ((pi.pv + pj.pv) * a * a * a * KERN_NORM + (pi.p + pj.p) * a * a * KERN)
                        / 2.0;

                    // relaxation
                    xproj -= d * dx / (r * pi.m);

                    // surface tension
                    xproj += (SURFACE_TENSION / pi.m) * pj.m * a * a * KERN * dx;

                    // linear and quadratic visc
                    let dv = pi.v - pj.v;
                    let mut u = dv.dot(dx);
                    if u > 0.0 {
                        u /= r;
                        let I = 0.5 * DT * a * (LINEAR_VISC * u + QUAD_VISC * u * u);
                        xproj -= I * dx * DT;
                    }
                }

                // correct
                pi.x = xproj;
                pi.v = (xproj - pi.xlast) / DT;

                // boundary
                for b in &bounds {
                    let d = f32::max(pi.x.x * b.x + pi.x.x * b.x - b.z, 0.0);
                    if d < PARTICLE_RADIUS {
                        pi.v += (PARTICLE_RADIUS - d) * Vec2::new(b.x, b.y) / DT;
                    }
                }
            })
    }

    pub fn update(&mut self) {
        for _ in 0..SOLVER_STEPS {
            self.integrate();
            //self.grid_insert();
            self.compute_forces();
            self.project_correct();
        }
    }
}

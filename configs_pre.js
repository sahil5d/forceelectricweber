// Objects are all built from elementary neutral segments
//   At runtime, all segments are flattened into the global `particles` array,
//   so `calc_forces()` and `plotallparticles()` operate directly on particles.
// Segment (`add_segment`): 1 proton (+q) + 1 electron (-q).
//   Why `prop_p` (physics) and `prop_draw_x/y/z` (visual) are separated:
//   A discrete wire or coil uses ~60-288 macro-charged segments to represent ~10^20 real electrons.
//   Shifting physical electron positions (`prop_p`) away from their protons would create unphysical
//   discrete multipoles and open-wire end dipoles. Keeping `prop_p` at the neutral segment center
//   preserves exact macroscopic electrostatic neutrality while `prop_draw_x/y/z` shows visual drift.
// Line (straight wire)
//   `add_line`, `pre_line_ac`, `build_lines_parallel`, `pre_line_tension`
//   n segments along a line
//   electrons advance along x and wrap around at the ends in `calc_each_pos()`
// Coil (circular wire) 
//   `add_coil`, `add_secondary_coil`)
//   n segments around a ring of radius R
//   `set_coil_state()` advances each electron's angle `prop_phi` around the circle
// Bar magnet
//   (`add_bar_magnet`, `pre_magnet_coil`)
//   5x5x9 rectangular prism of 225 small circular coils (`CFG.magnets_attract.NU * CFG.magnets_attract.NW * CFG.magnets_attract.NS`),
//   with the 2x2x9 = 36 corner coils active in `calc_forces()` and all 225 coils updated
//   as a single rigid body via `sync_bar_magnet()`

const CFG = {
    line: {
        N: 4800,
        L: 12, // m (1200 cm)
        DX: 12 / 4800,
        D0: 0.045, // m (4.5 cm), initial distance of lone electron from wire
        VD: 2.25e3, // m/s, wire electron drift v (+x, visible drift)
        V0: 1.2e5, // m/s, lone test electron velocity (+x, parallel to wire)
        Q: 2e11, // charge per wire segment
    },
    coil: {
        N: 144,
        R: 0.08, // m (8 cm), radius of circular coil in xy-plane (z = 0)
        VD: 9e3, // m/s, tangential electron speed around the coil (visible rotation)
        V0: 3e4, // m/s, lone test electron velocity (+x, through center of coil)
        Q: 5e11 / 3, // charge per coil segment
    },
    line_ac: {
        N: 480,
        L: 1.2, // m (120 cm)
        DX: 1.2 / 480,
        D0: 0.025, // m (2.5 cm), initial distance of lone electron from AC wire
        AMP: 0.0035, // m (0.35 cm), oscillation amplitude of AC wire electrons
        CYCLES: 1.5, // number of AC cycles over stepsmax
        Q: 2e12 / 3, // charge per AC wire segment
    },
    coil_turnon: {
        N: 108,
        R: 0.08, // m (8 cm), radius of the coil
        VD: 9e3, // m/s, steady-state electron speed after voltage is turned ON
        Q: 5e12 / 3, // charge per coil segment
        STEP: 600, // step at which voltage is turned ON and current ramps up
        TAU: 100 * 1e-9, // ramp time constant (s) after voltage is turned on
    },
    coils_ac: {
        N1: 3, // 3-loop primary coil
        N2: 1, // 1-loop secondary coil (3:1 step-down transformer)
        N_LOOP: 72,
        R1: 0.08, // m (8 cm), radius of 3-loop primary coil
        R2: 0.065, // m (6.5 cm), radius of 1-loop secondary coil
        DZ: 0.008, // m (0.8 cm), axial spacing between primary loops
        CYCLES: 1, // number of AC cycles over stepsmax
        AMP_PHI1: (2 * Math.PI) / 72, // angular AC amplitude (rad) of primary coil electrons
        Q1: 6.85e11, // charge per primary wire segment (chosen so 1-loop secondary electron speed is 3x the 3-loop primary speed)
    },
    lines_parallel: {
        N: 600,
        L: 1.2, // m (120 cm, identical length for both wires)
        DX: 1.2 / 600,
        D0: 0.04, // m (4 cm), initial distance between parallel Wire 1 (y = -2 cm) and Wire 2 (y = +2 cm)
        VD: 9e3, // m/s, equal electron drift speed magnitude in both wires
        Q: 1.1e16 / 3, // equal charge per segment in both wires (scaled for proton mass M_PROT)
    },
    magnets_attract: {
        NS: 9, // 9 coil layers along the bar magnet's length (dipole axis)
        NU: 5, // 5 coils across width of cross-section
        NW: 5, // 5 coils across height of cross-section
        ACTIVE_COILS: 36, // 2x2x9 = 36 corner coils active in force calculations
        N_LOOP: 8, // segments per small circular coil
        R: 0.0009, // m (0.09 cm), radius of each small circular coil in the 5x5x9 prism
        DS: 0.021 / 8, // m, axial spacing between slices (bar magnet length = 2.1 cm)
        DU: 0.021 / 8, // m, transverse grid spacing across the 5x5 cross-section
        D0: 0.06, // m (6.0 cm), initial center-to-center distance for attracting bar magnets (Config 8)
        VD: 1.35e4, // m/s, tangential electron speed around each bar magnet loop
        Q: (1.08e18 * (0.007 / 0.0009)) / (36 * 8), // charge per active segment (gives ~3.36 N bar-magnet force at 6 cm)
        PMASS: (5.4e20 * M_ELEC) / (36 * 8), // mass per active particle
    },
    magnets_repel: {
        D0: 0.046, // m (4.6 cm), closer initial distance for repelling bar magnets (Config 9)
        TILT: (6 * Math.PI) / 180, // 6 deg slight initial misalignment for Config 9
    },
    magnet_coil: {
        VZ: 3.6e4, // m/s, axial speed (+z) of bar magnet moving through coil in Config 10
        Z0: -0.06, // m (-6.0 cm), initial center z-position of bar magnet in Config 10
        Q: ((1.08e18 * (0.007 / 0.0009)) / (36 * 8)) * 1e-3, // charge per segment for moving magnet in Config 10
    },
    line_magnet: {
        Z0: -0.028, // m (-2.8 cm), z-position (< 0) of bar magnet on z-axis in Config 11
        L: 0.45, // m (45 cm), length of DC wire in Config 11
        N: 120, // segments along DC wire in Config 11
        VD: 1.35e4, // m/s, DC wire electron drift speed (+x) in Config 11
        Q: 2.3e15, // charge per DC wire segment in Config 11
    },
    line_tension: {
        N: 60, // segments along the wire
        L: 0.08, // m (8 cm), initial unstretched wire length
        DX: 0.08 / 60,
        VD: 1.35e4, // m/s, electron drift speed (+x) along wire
        Q: 1.28e15, // charge per wire segment (scaled for proton mass M_PROT)
        PMASS: 1.28e15 * M_PROT, // mass per particle in segment
    },
};

var par_vd2_signed = CFG.lines_parallel.VD;
var mag_states = [];

// Helper: create a neutral wire/coil segment (1 positive proton/ion +q and 1 negative electron -q)
function add_segment(id_base, q, shared_opts, elec_opts) {
    for (var sgn of [1, -1]) {
        add_particle(Object.assign(
            {
                prop_id: id_base + (sgn > 0 ? '_prot' : '_elec'),
                prop_charge: sgn * q,
            },
            shared_opts,
            sgn < 0 ? elec_opts : null
        ));
    }
}

// Helper: create a straight DC wire of n segments along the x-axis from -len/2 to +len/2
function add_line(n, len, vd, q, id_prefix) {
    var dx = len / n;
    for (var i = 1; i <= n; i++) {
        var x = -len / 2 + (i - 0.5) * dx;
        add_segment(
            id_prefix + i,
            q,
            {
                prop_kind: 'line',
                prop_iswire: true,
                prop_vd: vd,
                prop_wire_l: len,
                prop_p: new p(x, 0, 0),
            },
            {
                prop_draw_x: x,
                prop_v: new vec(vd, 0, 0),
            }
        );
    }
}

// Helper: create a circular coil of n segments in the xy-plane at height z
function add_coil(n, R, z, vd, q, id_prefix, extra_shared) {
    var xc = (extra_shared && extra_shared.prop_xc) || 0;
    var yc = (extra_shared && extra_shared.prop_yc) || 0;
    var vz = (extra_shared && extra_shared.prop_vz) || 0;
    for (var i = 0; i < n; i++) {
        var phi = (2 * Math.PI * i) / n;
        var px = xc + R * Math.cos(phi), py = yc + R * Math.sin(phi);
        add_segment(
            id_prefix + (i + 1),
            q,
            Object.assign(
                {
                    prop_kind: 'coil',
                    prop_iswire: true,
                    prop_r: R,
                    prop_vd: vd,
                    prop_phi0: phi,
                    prop_phi: phi,
                    prop_p: new p(px, py, z),
                    prop_v: new vec(0, 0, vz),
                },
                extra_shared
            ),
            {
                prop_v: new vec(-vd * Math.sin(phi), vd * Math.cos(phi), vz),
            }
        );
    }
}

function pre_line() {
    var thisconfig = CFG.line;
    // Lone test electron moving parallel to the wire at y = thisconfig.D0, starting at x_start to the right of the left HUD
    var travel = thisconfig.V0 * (stepsmax * timestep);
    var x_start = -0.02;
    add_particle({
        prop_id: 'elec_test',
        prop_p: new p(x_start, thisconfig.D0, 0),
        prop_v: new vec(thisconfig.V0, 0, 0),
        prop_charge: -1,
    });

    // Wire along x-axis from -thisconfig.L/2 to +thisconfig.L/2 made of protons and drifting electrons
    add_line(thisconfig.N, thisconfig.L, thisconfig.VD, thisconfig.Q, 'wire');

    // Reference line along +x at y = thisconfig.D0 showing the initial straight direction
    plotlinesegment(new vec(1, 0, 0), new p(x_start, thisconfig.D0, 0), 0.99 * travel, 0x686260);
    set_camera(0.07, 0.03, 0.38, 0.07, 0.03, 0);
}

function pre_coil() {
    var thisconfig = CFG.coil;
    // Lone test electron moving along +x right through the center (0, 0, 0) of the coil
    var x_start = -0.5 * thisconfig.V0 * (stepsmax * timestep);
    add_particle({
        prop_id: 'elec_test_coil',
        prop_p: new p(x_start, 0, 0),
        prop_v: new vec(thisconfig.V0, 0, 0),
        prop_charge: -1,
    });

    // Single circular coil in the xy-plane (z = 0)
    add_coil(thisconfig.N, thisconfig.R, 0, thisconfig.VD, thisconfig.Q, 'coil1');

    // Reference line along +x (y = 0) showing the undeflected straight path
    plotlinesegment(new vec(1, 0, 0), new p(x_start, 0, 0), -2 * x_start, 0x686260);
    set_camera(0, 0, 0.35, 0, 0, 0);
}

function pre_line_ac() {
    var thisconfig = CFG.line_ac;
    var omega = (2 * Math.PI * thisconfig.CYCLES) / (stepsmax * timestep);

    // Lone test electron initially at rest at (0, thisconfig.D0, 0) near the AC wire
    add_particle({
        prop_id: 'elec_test_ac',
        prop_p: new p(0, thisconfig.D0, 0),
        prop_charge: -1,
    });

    // Straight wire along x-axis with stationary protons and harmonically oscillating AC electrons:
    // x_e(t) = x0 - A * cos(omega * t), v_e(t) = A * omega * sin(omega * t), a_e(t) = A * omega^2 * cos(omega * t)
    var ax0 = thisconfig.AMP * (omega ** 2);
    for (var i = 1; i <= thisconfig.N; i++) {
        var x0 = -thisconfig.L / 2 + (i - 0.5) * thisconfig.DX;
        var xe = wrap_coord(x0 - thisconfig.AMP, thisconfig.L);
        add_segment(
            'ac' + i,
            thisconfig.Q,
            {
                prop_kind: 'line_ac',
                prop_iswire: true,
                prop_x0: x0,
                prop_p: new p(x0, 0, 0),
            },
            {
                prop_draw_x: xe,
                prop_a: new vec(ax0, 0, 0),
            }
        );
    }

    // Vertical reference line at x = 0 showing the test electron's starting x-position
    plotlinesegment(new vec(0, 1, 0), new p(0, thisconfig.D0, 0), 0.03, 0x686260);
    set_camera(0, 0.05, 0.32, 0, 0.05, 0);
}

function pre_coil_turnon() {
    var thisconfig = CFG.coil_turnon;
    // Lone test electron initially at rest at the same starting position as Config 2 inside the 1-loop coil
    var x_start = -0.5 * CFG.coil.V0 * (stepsmax * timestep);
    add_particle({
        prop_id: 'elec_test_coil_turnon',
        prop_p: new p(x_start, 0, 0),
        prop_charge: -1,
    });

    // Circular loop centered at z = 0, initially at rest before voltage turns ON at step thisconfig.STEP
    add_coil(thisconfig.N, thisconfig.R, 0, 0, thisconfig.Q, 'coil_turnon_', { prop_kind: 'coil_turnon' });

    // Reference line along +y from (x_start, 0, 0) showing the tangential (-phi, Lenz) induction direction
    plotlinesegment(new vec(0, 1, 0), new p(x_start, 0, 0), 0.015, 0x686260);
    set_camera(0, 0, 0.35, 0, 0, 0);
}

function add_secondary_coil() {
    var thisconfig = CFG.coils_ac;
    // Push the lead secondary electron FIRST (prop_iswire: false, index 0)
    // so calc_forces() computes Weber forces onto the secondary loop.
    // (by circular symmetry, all electrons on the secondary loop experience the same tangential Weber force as the lead electron)
    for (var i = 0; i < thisconfig.N_LOOP; i++) {
        var phi = (2 * Math.PI * i) / thisconfig.N_LOOP;
        var pos = new p(thisconfig.R2 * Math.cos(phi), thisconfig.R2 * Math.sin(phi), 0);
        var shared = { prop_group: 'sec', prop_r: thisconfig.R2, prop_p: pos };
        add_particle(Object.assign({
            prop_id: 'sec_l1_elec' + (i + 1),
            prop_kind: i === 0 ? 'coil_secondary_lead' : 'coil_secondary',
            prop_iswire: i > 0,
            prop_phi0: phi,
            prop_phi: phi,
            prop_vphi: 0,
            prop_charge: -1,
        }, shared));
        add_particle(Object.assign({
            prop_id: 'sec_l1_prot' + (i + 1),
            prop_kind: 'coil_secondary',
            prop_iswire: true,
            prop_charge: 1,
        }, shared));
    }
}

function pre_coils_ac() {
    var thisconfig = CFG.coils_ac;
    add_secondary_coil();

    // 3) Push the 3-loop primary coil (radius thisconfig.R1) driven with alternating current:
    //    phi_1(t) = phi0 - thisconfig.AMP_PHI1 * cos(omega * t)
    //    At t = 0, v_phi1(0) = 0 and tangential acceleration is at peak a_phi1(0) = thisconfig.R1 * thisconfig.AMP_PHI1 * omega^2
    var omega_tfr = (2 * Math.PI * thisconfig.CYCLES) / (stepsmax * timestep);
    var aphi1_0 = thisconfig.R1 * thisconfig.AMP_PHI1 * (omega_tfr ** 2);

    for (var k = 0; k < thisconfig.N1; k++) {
        var zk = (k - (thisconfig.N1 - 1) / 2) * thisconfig.DZ;
        for (var i = 0; i < thisconfig.N_LOOP; i++) {
            var phi0 = (2 * Math.PI * i) / thisconfig.N_LOOP;
            add_segment(
                'pri_l' + (k + 1) + '_' + (i + 1),
                thisconfig.Q1,
                {
                    prop_kind: 'coil_ac',
                    prop_iswire: true,
                    prop_r: thisconfig.R1,
                    prop_phi0: phi0,
                    prop_p: new p(thisconfig.R1 * Math.cos(phi0), thisconfig.R1 * Math.sin(phi0), zk),
                }
            );
            set_coil_state(particles[particles.length - 1], phi0 - thisconfig.AMP_PHI1, 0, aphi1_0);
        }
    }

    set_camera(0, 0, 0.35, 0, 0, 0);
}

function build_lines_parallel(vd2_dir) {
    var thisconfig = CFG.lines_parallel;
    par_vd2_signed = vd2_dir * thisconfig.VD;
    var y1_0 = -thisconfig.D0 / 2;
    var y2_0 = thisconfig.D0 / 2;

    // 1) Push central neutral lead atoms of BOTH Wire 1 (indices 0, 1 at y = -D0/2) and Wire 2 (indices 2, 3 at y = +D0/2):
    //    Both wires have identical length, segment spacing, charge per segment, and drift speed magnitude,
    //    so both wires carry the exact same current magnitude |I1| = |I2| and move the exact same distance!
    var wire_specs = [
        { id: 1, group: 'par1', lead_idx: 0, y0: y1_0, vd: thisconfig.VD },
        { id: 2, group: 'par2', lead_idx: 2, y0: y2_0, vd: par_vd2_signed },
    ];

    for (var ws of wire_specs) {
        add_segment(
            'par_w' + ws.id + '_lead',
            1,
            {
                prop_kind: 'line_plll_lead',
                prop_group: ws.group,
                prop_lead_idx: ws.lead_idx,
                prop_vd: ws.vd,
                prop_ystart: ws.y0,
                prop_p: new p(0, ws.y0, 0),
                prop_mass: M_PROT,
            },
            {
                prop_v: new vec(ws.vd, 0, 0),
            }
        );
    }

    // 2) Push the protons and drifting electrons of both wires
    for (var ws2 of wire_specs) {
        for (var j = 1; j <= thisconfig.N; j++) {
            var xj = -thisconfig.L / 2 + (j - 0.5) * thisconfig.DX;
            add_segment(
                'par_w' + ws2.id + '_' + j,
                thisconfig.Q,
                {
                    prop_kind: 'line_plll',
                    prop_iswire: true,
                    prop_group: ws2.group,
                    prop_lead_idx: ws2.lead_idx,
                    prop_vd: ws2.vd,
                    prop_p: new p(xj, ws2.y0, 0),
                },
                {
                    prop_draw_x: xj,
                    prop_v: new vec(ws2.vd, 0, 0),
                }
            );
        }
    }

    // Short grey reference lines centered at x = 0 showing both wires' initial vertical positions y = -D0/2 and y = +D0/2
    var par_ref_len = 0.08;
    plotlinesegment(new vec(1, 0, 0), new p(-par_ref_len / 2, y1_0, 0), par_ref_len, 0x686260);
    plotlinesegment(new vec(1, 0, 0), new p(-par_ref_len / 2, y2_0, 0), par_ref_len, 0x686260);
    set_camera(0, 0, 0.28, 0, 0, 0);
}

function pre_lines_parallel() {
    build_lines_parallel(1);
}

function pre_lines_antiparallel() {
    build_lines_parallel(-1);
}

function for_each_mag_seg(fn, include_visual) {
    var thisconfig = CFG.magnets_attract;
    var passes = include_visual ? [true, false] : [true];
    var coil_idx = 0;
    for (var pass_active of passes) {
        for (var ks = 0; ks < thisconfig.NS; ks++) {
            var s = (ks - (thisconfig.NS - 1) / 2) * thisconfig.DS;
            for (var ku = 0; ku < thisconfig.NU; ku++) {
                var is_corner_u = (ku === 0 || ku === thisconfig.NU - 1);
                var uc = (ku - (thisconfig.NU - 1) / 2) * thisconfig.DU;
                for (var kw = 0; kw < thisconfig.NW; kw++) {
                    var is_corner = is_corner_u && (kw === 0 || kw === thisconfig.NW - 1);
                    if (is_corner !== pass_active) continue;
                    var wc = (kw - (thisconfig.NW - 1) / 2) * thisconfig.DU;
                    coil_idx++;
                    for (var i = 0; i < thisconfig.N_LOOP; i++) {
                        fn(s, uc, wc, (2 * Math.PI * i) / thisconfig.N_LOOP, is_corner, coil_idx, i);
                    }
                }
            }
        }
    }
}

function add_bar_magnet(id, x0, theta0, vd, z0) {
    var thisconfig = CFG.magnets_attract;
    z0 = z0 || 0;
    var pmass = thisconfig.PMASS;
    var mag = {
        id: id, x0: x0, theta0: theta0,
        x: x0, y: 0, z: z0, vx: 0, vy: 0, ax: 0, ay: 0, ax_old: 0, ay_old: 0,
        theta: theta0, omega: 0, alpha: 0, alpha_old: 0, phi_shift: 0,
        vd: vd, fx: 0, fy: 0, tau: 0, fx0: 0, tau0: 0,
        mass: 0, inertia: 0, particles: [], all_particles: [],
    };
    for_each_mag_seg(function (s, uc, wc, phi, is_active, coil_idx, i) {
        if (is_active) {
            mag.mass += 2 * pmass;
            mag.inertia += 2 * pmass * (s * s + (uc + thisconfig.R * Math.cos(phi)) ** 2);
        }
        for (var sgn of [1, -1]) {
            var pt = add_particle({
                prop_id: 'mag' + id + '_c' + coil_idx + (sgn > 0 ? '_prot' : '_elec') + (i + 1),
                prop_kind: 'magnet',
                prop_group: 'mag' + id,
                prop_iswire: !is_active,
                prop_visual_only: !is_active,
                prop_s: s,
                prop_uc: uc,
                prop_wc: wc,
                prop_phi: phi,
                prop_charge: sgn * thisconfig.Q,
                prop_mass: pmass,
            });
            if (is_active) {
                mag.particles.push(pt);
            }
            mag.all_particles.push(pt);
        }
    }, true);
    sync_bar_magnet(mag, true);
    mag_states.push(mag);

    // Grey reference lines showing initial axis and perpendicular axis at each bar magnet center
    var ref_len = 0.04, ux = Math.cos(theta0), uy = Math.sin(theta0);
    plotlinesegment(new vec(ux, uy, 0), new p(x0 - 0.5 * ref_len * ux, -0.5 * ref_len * uy, z0), ref_len, 0x686260);
}

function build_bar_magnets(vd2_dir, tilt_rad, d0) {
    var thisconfig = CFG.magnets_attract;
    add_bar_magnet(1, -d0 / 2, tilt_rad, thisconfig.VD);
    add_bar_magnet(2, d0 / 2, tilt_rad, vd2_dir * thisconfig.VD);
    set_camera(0, -0.1, 0.24, 0, 0, 0);
}

function pre_magnets_attract() {
    var thisconfig = CFG.magnets_attract;
    build_bar_magnets(1, 0, thisconfig.D0);
}

function pre_magnets_repel() {
    var thisconfig = CFG.magnets_repel;
    // Initially N-to-N (S->N and N<-S): both magnets repel and torque counterclockwise toward theta = 90 deg
    // (side-by-side anti-parallel equilibrium). Over longer runs, angular momentum carries them past 90 deg to 180 deg
    // (S-to-S repulsion), where the 1/r^4 force and 1/r^3 torque are much weaker as the magnets have already moved apart.
    build_bar_magnets(-1, thisconfig.TILT, thisconfig.D0);
}

function pre_magnet_coil() {
    var thisconfig = CFG.magnet_coil;
    // 1) 1-loop conducting secondary coil in the xy-plane (z = 0, radius CFG.coils_ac.R2)
    add_secondary_coil();

    // 2) 5x5x9 bar magnet moving along +z through the coil (36 corner coils active first, 189 visual coils second)
    for_each_mag_seg(function (s, uc, wc, phi, is_active, coil_idx, i) {
        if (i === 0) {
            add_coil(CFG.magnets_attract.N_LOOP, CFG.magnets_attract.R, thisconfig.Z0 + s, CFG.magnets_attract.VD, thisconfig.Q, 'mcoil_c' + coil_idx + '_', {
                prop_visual_only: !is_active,
                prop_xc: uc,
                prop_yc: wc,
                prop_vz: thisconfig.VZ,
            });
        }
    }, true);

    // Reference line along the z-axis through the center of the coil showing the bar magnet's path
    plotlinesegment(new vec(0, 0, 1), new p(0, 0, thisconfig.Z0 - 0.015), 0.08, 0x686260);
    set_camera(0.2, 0, 0.2, 0, 0, 0);
}

function pre_line_magnet() {
    var thisconfig = CFG.line_magnet;
    // 1) Bar magnet on the z-axis at (0, 0, thisconfig.Z0) with z < 0, initially parallel to the wire (theta0 = 0 along x-axis)
    add_bar_magnet(1, 0, 0, CFG.magnets_attract.VD, thisconfig.Z0);

    // 2) Straight DC wire along the x-axis (y = 0, z = 0)
    add_line(thisconfig.N, thisconfig.L, thisconfig.VD, thisconfig.Q, 'oersted');

    set_camera(0.05, 0.15, 0.15, 0, 0, 0);
}

function pre_line_tension() {
    var thisconfig = CFG.line_tension;
    // Straight wire along x-axis where all segments interact pairwise via calc_forces():
    // Each segment k has its own prop_group ('tseg' + k) and prop_iswire: false,
    // so same-segment proton & electron skip each other while all distinct segments repel longitudinally.
    // Note: each segment's proton and electron receive equal half-shares (+kq^2/r^2 * vd^2/2c^2) of the longitudinal
    // repulsion and are coupled together at x_new with proton mass (thisconfig.PMASS), since any charge separation
    // in a real wire is immediately locked together by the much stronger (c^2/vd^2 ~ 500x) Coulomb restoring force.
    for (var k = 1; k <= thisconfig.N; k++) {
        var x0 = -thisconfig.L / 2 + (k - 0.5) * thisconfig.DX;
        var lead_idx = 2 * (k - 1);
        add_segment(
            'tseg' + k,
            thisconfig.Q,
            {
                prop_kind: 'line_tension',
                prop_group: 'tseg' + k,
                prop_iswire: false,
                prop_lead_idx: lead_idx,
                prop_x0: x0,
                prop_vd: thisconfig.VD,
                prop_p: new p(x0, 0, 0),
                prop_v: new vec(0, 0, 0),
                prop_mass: thisconfig.PMASS,
            },
            {
                prop_v: new vec(thisconfig.VD, 0, 0),
            }
        );
    }

    // Grey reference line and end ticks at x = +-thisconfig.L/2 showing initial unstretched wire length
    plotlinesegment(new vec(1, 0, 0), new p(-thisconfig.L / 2, 0, 0), thisconfig.L, 0x686260);
    plotlinesegment(new vec(0, 1, 0), new p(-thisconfig.L / 2, -0.006, 0), 0.012, 0x686260);
    plotlinesegment(new vec(0, 1, 0), new p(thisconfig.L / 2, -0.006, 0), 0.012, 0x686260);
    set_camera(0, 0, 0.32, 0, 0, 0);
}

/*
const SHELL_N = 500;
const SHELL_R = 7; // m, radius of hollow charged spherical shell
const SHELL_Q_PER_PT = 7e12; // positive charge per shell vertex (500 pts) at meter scale for Weber F_shell
const SHELL_EXT_X = 25; // m, symmetric external dipole/field charges at (+-SHELL_EXT_X, 0, 0)
const SHELL_EXT_Q = 4e13; // external charge pulling test electron along +x at meter scale

function pre_spherical_shell() {
    // 1) Lone test electron initially at rest at the exact center (0, 0, 0) of the hollow sphere
    add_particle({
        prop_id: 'elec_test_shell',
        prop_history: [new p(0, 0, 0)],
        prop_p: new p(0, 0, 0),
        prop_charge: -1,
    });

    // 2) Build a 500-point spherical shell of radius SHELL_R around (0, 0, 0) using getapproxequidistantpointsonsphere:
    var shell_pts = getapproxequidistantpointsonsphere(SHELL_N, SHELL_R);
    for (var k = 0; k < shell_pts.length; k++) {
        add_particle({
            prop_id: 'shell_prot' + (k + 1),
            prop_iswire: true,
            prop_isshell: true,
            prop_p: shell_pts[k],
            prop_charge: SHELL_Q_PER_PT,
        });
    }

    // 3) Symmetric external field charges at (+SHELL_EXT_X, 0, 0) with +SHELL_EXT_Q
    //    and (-SHELL_EXT_X, 0, 0) with -SHELL_EXT_Q:
    //    Their Coulomb forces add along +x (F_ext = 2 * KEE * SHELL_EXT_Q / SHELL_EXT_X^2),
    //    while their 1/r Weber F_shell terms (+Q/R and -Q/R) cancel to zero at the origin!
    add_particle({
        prop_id: 'ext_attractor_pos',
        prop_iswire: true,
        prop_p: new p(SHELL_EXT_X, 0, 0),
        prop_charge: SHELL_EXT_Q,
    });
    add_particle({
        prop_id: 'ext_repulsor_neg',
        prop_iswire: true,
        prop_p: new p(-SHELL_EXT_X, 0, 0),
        prop_charge: -SHELL_EXT_Q,
        prop_mass: M_PROT,
    });

    plotlinesegment(new vec(1, 0, 0), new p(0, 0, 0), SHELL_R, 0x686260);
    set_camera(0, -18, 14, 0, 0, 0);
}
*/

var CONFIG_PRE = {
    line: pre_line,
    coil: pre_coil,
    line_ac: pre_line_ac,
    coil_turnon: pre_coil_turnon,
    coils_ac: pre_coils_ac,
    lines_parallel: pre_lines_parallel,
    lines_antiparallel: pre_lines_antiparallel,
    magnets_attract: pre_magnets_attract,
    magnets_repel: pre_magnets_repel,
    magnet_coil: pre_magnet_coil,
    line_magnet: pre_line_magnet,
    line_tension: pre_line_tension,
    // spherical_shell: pre_spherical_shell
};

var config_descriptions = {
    line: '• Lone electron (v = 1.2e5 m/s) parallel to DC wire at y = 4.5 cm\n• Direct Weber velocity force (vᵣ²/2c², v⊥²/c² terms)\n• Replaces magnetic field B and F = qv x B',
    coil: '• Lone electron (v = 3e4 m/s) through center of coil\n• Direct Weber pairwise force from rotating coil electrons\n• Replaces coil B-field and F = qv x B',
    line_ac: '• Stationary electron at y = 2.5 cm near 120 cm AC wire\n• Direct Weber acc force (r aᵣ/c² term)\n• Replaces Faraday inductive E = −∂A/∂t',
    coil_turnon: '• Stationary electron in 1-loop coil, turn-on at t = 600 ns\n• Direct Weber acceleration kick from ramping coil electrons\n• Replaces Faraday/Lenz changing-flux B-field induction',
    coils_ac: '• 3-loop AC primary and 1-loop secondary\n• Direct primary-to-secondary Weber acceleration force\n• Replaces mag flux (V2/V1 = N2/N1 = 1/3)',
    lines_parallel: '• Two 120 cm wires at y = ±2 cm, same-direction currents (|I1| = |I2|)\n• Direct pairwise Weber attraction (both move equally toward y = 0)\n• Replaces wire B-field and Ampère force',
    lines_antiparallel: '• Two 120 cm wires at y = ±2 cm, opp-direction currents (|I1| = |I2|)\n• Direct pairwise Weber repulsion (both move equally apart)\n• Replaces wire B-field and Ampère force',
    magnets_attract: '• Two 5×5×9-coil bar magnets at x = ±3.0 cm, aligned (N–S)\n• Direct pairwise Weber attraction between circulating Ampèrian loops\n• Replaces magnetic pole / Biot-Savart dipole attraction',
    magnets_repel: '• Two 5×5×9-coil bar magnets at x = ±2.3 cm, opp currents (N–N) tilted 6°\n• Direct pairwise Weber torque and repulsion rotates and pushes apart\n• Replaces mag dipole torque τ = μ × B',
    magnet_coil: '• 5×5×9-coil bar magnet (vz = 3.6e4 m/s) moving through coil (R = 6.5 cm)\n• Direct Weber velocity force (vᵣ²/2c², v⊥²/c²) drives coil current\n• Replaces Faraday induction E = −dΦ_B/dt',
    line_magnet: '• 5×5×9-coil bar magnet at (0, 0, −2.8 cm) below DC wire (Ørsted exp)\n• Direct pairwise Weber torque twists magnet to be perp\n• Replaces wire B-field and compass torque τ = μ × B',
    line_tension: '• 60-segment DC wire (L = 8 cm) with inter-segment Weber forces\n• Collinear segments repel via −vᵣ²/2c², stretching & exploding wire\n• Ampère–Weber longitudinal tension (Lorentz I dl × B predicts 0)',
    // spherical_shell: '• Stationary electron inside 500-point charged sphere (R = 7 cm)\n• Direct Weber acceleration reaction from shell charges\n• Produces F_shell = -(Qq / 3c²R) a'
};

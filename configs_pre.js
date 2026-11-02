// note: what is p
//     point, pos, or particle.....
// particle is the {prop_id, prop_p, prop_v, ...}.
// p is a point is a pos is a {x: , y: , z:}

// Helper: create a neutral wire/coil segment (1 positive proton/ion +q and 1 negative electron -q)
function add_wire_segment(id_base, q, shared_opts, elec_opts) {
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
function add_dc_wire(n, len, vd, q, id_prefix) {
    var dx = len / n;
    for (var i = 1; i <= n; i++) {
        var x = -len / 2 + (i - 0.5) * dx;
        add_wire_segment(
            id_prefix + i,
            q,
            {
                prop_iswire: true,
                prop_isdcwire: true,
                prop_vd: vd,
                prop_wire_l: len,
                prop_p: new p(x, 0, 0),
            },
            {
                prop_v: new vec(vd, 0, 0),
            }
        );
    }
}

// Helper: create a circular coil loop of n segments in the xy-plane at height z
function add_coil_loop(n, R, z, vd, q, id_prefix, extra_shared) {
    for (var i = 0; i < n; i++) {
        var phi = (2 * Math.PI * i) / n;
        var px = R * Math.cos(phi), py = R * Math.sin(phi);
        var vz = (extra_shared && extra_shared.prop_vz) || 0;
        add_wire_segment(
            id_prefix + (i + 1),
            q,
            Object.assign(
                {
                    prop_iswire: true,
                    prop_r: R,
                    prop_vd: vd,
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

const N_WIRE = 3600;
const WIRE_L = 600; // m
const WIRE_DX = WIRE_L / N_WIRE;
const WIRE_D0 = 4.5; // m, initial distance of lone electron from wire
const WIRE_VD = 1.5e-3 * CSPEED; // wire electron drift v (+x, vd << v0)
const TEST_V0 = 0.1 * CSPEED; // lone test electron velocity (+x, parallel to wire)
const WIRE_Q = 2e14 / 3; // charge per wire segment at meter scale (1/3 per particle with 3x particles)

function pre_wire() {
    // Lone test electron moving parallel to the wire at y = WIRE_D0, centered across x = 0
    var x_start = -0.5 * TEST_V0 * (stepsmax * timestep);
    add_particle({
        prop_id: 'elec_test',
        // prop_history: [new p(x_start, WIRE_D0, 0)],
        prop_p: new p(x_start, WIRE_D0, 0),
        prop_v: new vec(TEST_V0, 0, 0),
        prop_charge: -1,
    });

    // Wire along x-axis from -WIRE_L/2 to +WIRE_L/2 made of protons and drifting electrons
    add_dc_wire(N_WIRE, WIRE_L, WIRE_VD, WIRE_Q, 'wire');

    // Reference line along +x at y = WIRE_D0 showing the initial straight direction
    plotlinesegment(new vec(1, 0, 0), new p(x_start, WIRE_D0, 0), -1.99 * x_start, 0x686260);
    set_camera(0, 2, 18, 0, 2, 0);
}

const N_PER_COIL = 144;
const COIL_R = 8; // m, radius of circular coil in xy-plane (z = 0)
const COIL_VD = 0.015 * CSPEED; // tangential electron speed around the coil (visible rotation)
const COIL_V0 = 0.05 * CSPEED; // lone test electron velocity (+x, through center of coil)
const COIL_Q = 5e13 / 3; // charge per coil segment at meter scale (1/3 per particle with 3x particles)

function pre_coil() {
    // Lone test electron moving along +x right through the center (0, 0, 0) of the coil
    var x_start = -0.5 * COIL_V0 * (stepsmax * timestep);
    add_particle({
        prop_id: 'elec_test_coil',
        prop_xstart: x_start,
        // prop_history: [new p(x_start, 0, 0)],
        prop_p: new p(x_start, 0, 0),
        prop_v: new vec(COIL_V0, 0, 0),
        prop_charge: -1,
    });

    // Single circular coil in the xy-plane (z = 0)
    add_coil_loop(N_PER_COIL, COIL_R, 0, COIL_VD, COIL_Q, 'coil1', { prop_iscoil: true });

    // Reference line along +x (y = 0) showing the undeflected straight path
    plotlinesegment(new vec(1, 0, 0), new p(x_start, 0, 0), -2 * x_start, 0x686260);
    set_camera(0, 0, 35, 0, 0, 0);
}

const N_AC_WIRE = 480;
const AC_WIRE_L = 120; // m
const AC_WIRE_DX = AC_WIRE_L / N_AC_WIRE;
const AC_WIRE_D0 = 2.5; // m, initial distance of lone electron from AC wire
const AC_AMP = 0.8; // m, oscillation amplitude of AC wire electrons
const AC_CYCLES = 2.5; // number of AC cycles over stepsmax
const AC_Q = 7.5e13 / 3; // charge per AC wire segment at meter scale (1/3 per particle with 3x particles)

function pre_ac_wire() {
    var omega = (2 * Math.PI * AC_CYCLES) / (stepsmax * timestep);

    // Lone test electron initially at rest at (0, AC_WIRE_D0, 0) near the AC wire
    add_particle({
        prop_id: 'elec_test_ac',
        // prop_history: [new p(0, AC_WIRE_D0, 0)],
        prop_p: new p(0, AC_WIRE_D0, 0),
        prop_charge: -1,
    });

    // Straight wire along x-axis with stationary protons and harmonically oscillating AC electrons:
    // x_e(t) = x0 - A * cos(omega * t), v_e(t) = A * omega * sin(omega * t), a_e(t) = A * omega^2 * cos(omega * t)
    var ax0 = AC_AMP * (omega ** 2);
    for (var i = 1; i <= N_AC_WIRE; i++) {
        var x0 = -AC_WIRE_L / 2 + (i - 0.5) * AC_WIRE_DX;
        var xe = wrap_coord(x0 - AC_AMP, AC_WIRE_L);
        add_wire_segment(
            'ac' + i,
            AC_Q,
            {
                prop_iswire: true,
                prop_isacwire: true,
                prop_x0: x0,
                prop_p: new p(x0, 0, 0),
            },
            {
                prop_p: new p(xe, 0, 0),
                prop_a: new vec(ax0, 0, 0),
                prop_anext: new vec(ax0, 0, 0),
            }
        );
    }

    // Vertical reference line at x = 0 showing the test electron's starting x-position
    plotlinesegment(new vec(0, 1, 0), new p(0, AC_WIRE_D0, 0), 3, 0x686260);
    set_camera(0, 5, 32, 0, 5, 0);
}

const COIL_TURNON_N = 108;
const COIL_TURNON_R = 8; // m, radius of the coil
const COIL_TURNON_VD = 0.015 * CSPEED; // steady-state electron speed after voltage is turned ON
const COIL_TURNON_Q = 5e14 / 3; // charge per coil segment at meter scale (1/3 per particle with 3x particles)
const COIL_TURNON_R0 = 3; // m, initial position (r0, 0, 0) of stationary test electron inside coil
const COIL_TURNON_STEP = 600; // step at which voltage is turned ON and current ramps up
const COIL_TURNON_TAU = 100 * 2e-10; // ramp time constant (s) after voltage is turned on

function pre_coil_turnon() {
    // Lone test electron initially at rest at (COIL_TURNON_R0, 0, 0) inside the 1-loop coil
    add_particle({
        prop_id: 'elec_test_coil_turnon',
        // prop_history: [new p(COIL_TURNON_R0, 0, 0)],
        prop_p: new p(COIL_TURNON_R0, 0, 0),
        prop_charge: -1,
    });

    // Circular loop centered at z = 0, initially at rest before voltage turns ON at step COIL_TURNON_STEP
    add_coil_loop(COIL_TURNON_N, COIL_TURNON_R, 0, 0, COIL_TURNON_Q, 'coil_turnon_', { prop_iscoil_turnon: true });

    // Reference line along -y from (COIL_TURNON_R0, 0, 0) showing the tangential (-phi, Lenz's law) induction direction
    plotlinesegment(new vec(0, -1, 0), new p(COIL_TURNON_R0, 0, 0), 1.5, 0x686260);
    set_camera(0, 0, 35, 0, 0, 0);
}

const TFR_N1 = 3; // 3-loop primary coil
const TFR_N2 = 1; // 1-loop secondary coil (3:1 step-down transformer)
const TFR_N_PER_LOOP = 72;
const TFR_R1 = 8; // m, radius of 3-loop primary coil
const TFR_R2 = 6.5; // m, radius of 1-loop secondary coil
const TFR_DZ = 0.8; // m, axial spacing between primary loops
const TFR_CYCLES = 2; // number of AC cycles over stepsmax
const TFR_AMP_PHI1 = (2 * Math.PI) / 24; // angular AC amplitude (rad) of primary coil electrons
const TFR_Q1 = 8e13 / 3; // charge per primary wire segment at meter scale (1/3 per particle with 3x particles)

function add_secondary_coil() {
    // Push the lead secondary electron FIRST (prop_iswire: false, index 0)
    // so calc_forces() computes Weber forces onto the secondary loop.
    // (by circular symmetry, all electrons on the secondary loop experience the same tangential Weber force as the lead electron)
    for (var i = 0; i < TFR_N_PER_LOOP; i++) {
        var phi = (2 * Math.PI * i) / TFR_N_PER_LOOP;
        var pos = new p(TFR_R2 * Math.cos(phi), TFR_R2 * Math.sin(phi), 0);
        var shared = { prop_group: 'sec', prop_issecondary: true, prop_r: TFR_R2, prop_p: pos };
        add_particle(Object.assign({
            prop_id: 'sec_l1_elec' + (i + 1),
            prop_iswire: i > 0,
            prop_is_sec_lead: i === 0,
            prop_phi0: phi,
            prop_phi: phi,
            prop_vphi: 0,
            // prop_history: i === 0 ? [pos.clone()] : undefined,
            prop_charge: -1,
        }, shared));
        add_particle(Object.assign({
            prop_id: 'sec_l1_prot' + (i + 1),
            prop_iswire: true,
            prop_charge: 1,
        }, shared));
    }
}

function pre_transformer() {
    add_secondary_coil();

    // 3) Push the 3-loop primary coil (radius TFR_R1) driven with alternating current:
    //    phi_1(t) = phi0 - TFR_AMP_PHI1 * cos(omega * t)
    //    At t = 0, v_phi1(0) = 0 and tangential acceleration is at peak a_phi1(0) = TFR_R1 * TFR_AMP_PHI1 * omega^2
    var omega_tfr = (2 * Math.PI * TFR_CYCLES) / (stepsmax * timestep);
    var aphi1_0 = TFR_R1 * TFR_AMP_PHI1 * (omega_tfr ** 2);

    for (var k = 0; k < TFR_N1; k++) {
        var zk = (k - (TFR_N1 - 1) / 2) * TFR_DZ;
        for (var i = 0; i < TFR_N_PER_LOOP; i++) {
            var phi0 = (2 * Math.PI * i) / TFR_N_PER_LOOP;
            add_wire_segment(
                'pri_l' + (k + 1) + '_' + (i + 1),
                TFR_Q1,
                {
                    prop_iswire: true,
                    prop_isprimary: true,
                    prop_r: TFR_R1,
                    prop_phi0: phi0,
                    prop_p: new p(TFR_R1 * Math.cos(phi0), TFR_R1 * Math.sin(phi0), zk),
                }
            );
            set_coil_state(particles[particles.length - 1], TFR_R1, phi0 - TFR_AMP_PHI1, 0, aphi1_0);
        }
    }

    set_camera(0, 0, 35, 0, 0, 0);
}

const PAR_N = 600;
const PAR_L = 120; // m (identical length for both wires)
const PAR_DX = PAR_L / PAR_N;
const PAR_D0 = 4; // m, initial distance between parallel Wire 1 (y = -2 m) and Wire 2 (y = +2 m)
const PAR_VD = 0.015 * CSPEED; // equal electron drift speed magnitude in both wires
const PAR_Q = 6e14 / 3; // equal charge per segment in both wires at meter scale (1/3 per particle with 3x particles)
var par_vd2_signed = PAR_VD;

function build_parallel_wires(vd2_dir) {
    par_vd2_signed = vd2_dir * PAR_VD;
    var y1_0 = -PAR_D0 / 2;
    var y2_0 = PAR_D0 / 2;

    // 1) Push central neutral lead atoms of BOTH Wire 1 (indices 0, 1 at y = -PAR_D0/2) and Wire 2 (indices 2, 3 at y = +PAR_D0/2):
    //    Both wires have identical length, segment spacing, charge per segment (PAR_Q), and drift speed magnitude (PAR_VD),
    //    so both wires carry the exact same current magnitude |I1| = |I2| and move the exact same distance!
    var wire_specs = [
        { id: 1, group: 'par1', lead_idx: 0, y0: y1_0, vd: PAR_VD },
        { id: 2, group: 'par2', lead_idx: 2, y0: y2_0, vd: par_vd2_signed },
    ];

    for (var ws of wire_specs) {
        add_wire_segment(
            'par_w' + ws.id + '_lead',
            1,
            {
                prop_group: ws.group,
                prop_ispar: true,
                prop_is_par_lead: true,
                prop_wire_id: ws.id,
                prop_lead_idx: ws.lead_idx,
                prop_vd: ws.vd,
                prop_ystart: ws.y0,
                prop_p: new p(0, ws.y0, 0),
                prop_mass: M_ELEC,
            },
            {
                prop_v: new vec(ws.vd, 0, 0),
            }
        );
    }

    // 2) Push the protons and drifting electrons of both wires
    for (var w = wire_specs.length - 1; w >= 0; w--) {
        var ws2 = wire_specs[w];
        for (var j = 1; j <= PAR_N; j++) {
            var xj = -PAR_L / 2 + (j - 0.5) * PAR_DX;
            add_wire_segment(
                'par_w' + ws2.id + '_' + j,
                PAR_Q,
                {
                    prop_iswire: true,
                    prop_group: ws2.group,
                    prop_ispar: true,
                    prop_wire_id: ws2.id,
                    prop_lead_idx: ws2.lead_idx,
                    prop_vd: ws2.vd,
                    prop_p: new p(xj, ws2.y0, 0),
                },
                {
                    prop_v: new vec(ws2.vd, 0, 0),
                }
            );
        }
    }

    // Short grey reference lines centered at x = 0 showing both wires' initial vertical positions y = -PAR_D0/2 and y = +PAR_D0/2
    var par_ref_len = 8;
    plotlinesegment(new vec(1, 0, 0), new p(-par_ref_len / 2, y1_0, 0), par_ref_len, 0x686260);
    plotlinesegment(new vec(1, 0, 0), new p(-par_ref_len / 2, y2_0, 0), par_ref_len, 0x686260);
    set_camera(0, 0, 28, 0, 0, 0);
}

function pre_parallel_wires() {
    build_parallel_wires(1);
}

function pre_antiparallel_wires() {
    build_parallel_wires(-1);
}

const MAG_LOOPS = 3;
const MAG_N_PER_LOOP = 10;
const MAG_R = 1.4; // m, radius of each bar magnet loop
const MAG_DZ = 1.05; // m, axial spacing between loops (bar magnet length = 2.1 m)
const MAG_D0 = 6.26; // m, initial center-to-center distance for attracting bar magnets (Config 8)
const MAG_D0_REPEL = 4.6; // m, closer initial distance for repelling bar magnets (Config 9)
const MAG_VD = 0.045 * CSPEED; // tangential electron speed around each bar magnet loop
const MAG_Q = 1.8e15; // charge per bar magnet segment
const MAG_TILT = (6 * Math.PI) / 180; // 6 deg slight initial misalignment for Config 9
const MCOIL_VZ = 0.06 * CSPEED; // axial speed (+z) of bar magnet moving through coil in Config 10
const MCOIL_Z0 = -6.0; // initial center z-position of bar magnet in Config 10
const OERSTED_Z0 = -2.8; // m, z-position (< 0) of bar magnet on z-axis in Config 11
const OERSTED_L = 45; // m, length of DC wire in Config 11
const OERSTED_N = 120; // segments along DC wire in Config 11
const OERSTED_DX = OERSTED_L / OERSTED_N;
const OERSTED_VD = 0.045 * CSPEED; // DC wire electron drift speed (+x) in Config 11
const OERSTED_Q = 1.15e14; // charge per DC wire segment in Config 11

var mag_states = [];

function for_each_mag_seg(fn) {
    for (var k = 0; k < MAG_LOOPS; k++) {
        var s = (k - (MAG_LOOPS - 1) / 2) * MAG_DZ;
        for (var i = 0; i < MAG_N_PER_LOOP; i++) {
            fn(k, i, s, (2 * Math.PI * i) / MAG_N_PER_LOOP);
        }
    }
}

function sync_bar_magnet(mag) {
    var cosT = Math.cos(mag.theta), sinT = Math.sin(mag.theta);
    for (var pt of mag.particles) {
        var cosP = Math.cos(pt.prop_phi), sinP = Math.sin(pt.prop_phi);
        var dx = pt.prop_s * cosT - MAG_R * cosP * sinT;
        var dy = pt.prop_s * sinT + MAG_R * cosP * cosT;
        pt.prop_dx = dx;
        pt.prop_dy = dy;
        pt.prop_p.set(mag.x + dx, mag.y + dy, mag.z + MAG_R * sinP);
        var vrx = mag.vx - mag.omega * dy, vry = mag.vy + mag.omega * dx;
        if (pt.prop_charge > 0) {
            pt.prop_v.set(vrx, vry, 0);
        } else {
            pt.prop_v.set(vrx + mag.vd * sinP * sinT, vry - mag.vd * sinP * cosT, mag.vd * cosP);
        }
    }
}

function add_bar_magnet(id, x0, theta0, vd, z0) {
    z0 = z0 || 0;
    var pmass = MAG_Q * M_ELEC;
    var mag = {
        id: id, x0: x0, theta0: theta0,
        x: x0, y: 0, z: z0, vx: 0, vy: 0, ax: 0, ay: 0, ax_old: 0, ay_old: 0,
        theta: theta0, omega: 0, alpha: 0, alpha_old: 0,
        vd: vd, fx: 0, fy: 0, tau: 0, fx0: 0, fy0: 0, tau0: 0,
        mass: 0, inertia: 0, particles: [],
    };
    for_each_mag_seg(function (k, i, s, phi) {
        mag.mass += 2 * pmass;
        mag.inertia += 2 * pmass * (s * s + (MAG_R * Math.cos(phi)) ** 2);
        for (var sgn of [1, -1]) {
            mag.particles.push(add_particle({
                prop_id: 'mag' + id + '_l' + (k + 1) + (sgn > 0 ? '_prot' : '_elec') + (i + 1),
                prop_group: 'mag' + id,
                prop_ismagnet: true,
                prop_mag_id: id,
                prop_s: s,
                prop_phi: phi,
                prop_charge: sgn * MAG_Q,
                prop_mass: pmass,
            }));
        }
    });
    sync_bar_magnet(mag);
    mag_states.push(mag);

    // Grey reference lines showing initial axis and perpendicular axis at each bar magnet center
    var ref_len = 4, ux = Math.cos(theta0), uy = Math.sin(theta0);
    plotlinesegment(new vec(ux, uy, 0), new p(x0 - 0.5 * ref_len * ux, -0.5 * ref_len * uy, z0), ref_len, 0x686260);
    // plotlinesegment(new vec(0, 1, 0), new p(x0, -MAG_R * 1.3, z0), MAG_R * 2.6, 0x686260);
}

function build_bar_magnets(vd2_dir, tilt_rad, d0) {
    mag_states = [];
    add_bar_magnet(1, -d0 / 2, tilt_rad, MAG_VD);
    add_bar_magnet(2, d0 / 2, tilt_rad, vd2_dir * MAG_VD);
    set_camera(0, -10, 24, 0, 0, 0);
}

function pre_magnets_attract() {
    build_bar_magnets(1, 0, MAG_D0);
}

function pre_magnets_repel() {
    build_bar_magnets(-1, MAG_TILT, MAG_D0_REPEL);
}

function pre_magnet_coil() {
    mag_states = [];
    // 1) 1-loop conducting secondary coil in the xy-plane (z = 0, radius TFR_R2)
    add_secondary_coil();

    // 2) Bar magnet (MAG_LOOPS coaxial circular loops of radius MAG_R) moving along +z through the coil
    for (var k = 0; k < MAG_LOOPS; k++) {
        var s = (k - (MAG_LOOPS - 1) / 2) * MAG_DZ;
        add_coil_loop(MAG_N_PER_LOOP, MAG_R, MCOIL_Z0 + s, MAG_VD, MAG_Q, 'mcoil_l' + (k + 1) + '_', {
            prop_iscoil: true,
            prop_vz: MCOIL_VZ,
        });
    }

    // Reference line along the z-axis through the center of the coil showing the bar magnet's path
    plotlinesegment(new vec(0, 0, 1), new p(0, 0, MCOIL_Z0 - 1.5), 8, 0x686260);
    set_camera(0, -45, 5, 0, 0, -1.5);
}

function pre_wire_magnet() {
    mag_states = [];
    // 1) Bar magnet on the z-axis at (0, 0, OERSTED_Z0) with z < 0, initially parallel to the wire (theta0 = 0 along x-axis)
    add_bar_magnet(1, 0, 0, MAG_VD, OERSTED_Z0);

    // 2) Straight DC wire along the x-axis (y = 0, z = 0)
    add_dc_wire(OERSTED_N, OERSTED_L, OERSTED_VD, OERSTED_Q, 'oersted');

    set_camera(0, 15, 15, 0, 0, 0);
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
    wire: pre_wire,
    coil: pre_coil,
    ac_wire: pre_ac_wire,
    coil_turnon: pre_coil_turnon,
    transformer: pre_transformer,
    parallel_wires: pre_parallel_wires,
    antiparallel_wires: pre_antiparallel_wires,
    magnets_attract: pre_magnets_attract,
    magnets_repel: pre_magnets_repel,
    magnet_coil: pre_magnet_coil,
    wire_magnet: pre_wire_magnet,
    // spherical_shell: pre_spherical_shell
};

var config_descriptions = {
    wire: '• Lone electron (v = 0.1c) parallel to DC wire at y = 4.5 m\n• Direct Weber velocity force (vᵣ²/2c², v⊥²/c² terms)\n• Replaces magnetic field B and F = qv x B',
    coil: '• Lone electron (v = 0.05c) through center of circular coil (R = 8 m)\n• Direct Weber pairwise force from rotating coil electrons\n• Replaces circular coil B-field and F = qv x B',
    ac_wire: '• Stationary electron at y = 2.5 m near 120 m AC wire\n• Direct Weber acceleration force (r aᵣ/c² term)\n• Replaces Faraday inductive field E = −∂A/∂t',
    coil_turnon: '• Stationary electron in 1-loop coil (R = 8 m), turn-on at step 600\n• Direct Weber acceleration kick from ramping coil electrons\n• Replaces Faraday/Lenz changing-flux B-field induction',
    transformer: '• 3-loop AC primary (R1 = 8 m) and 1-loop secondary (R2 = 6.5 m)\n• Direct primary-to-secondary Weber acceleration force\n• Replaces magnetic flux linkage (V2/V1 = N2/N1 = 1/3)',
    parallel_wires: '• Two 120 m wires at y = ±2 m, same-direction currents (|I1| = |I2|)\n• Direct pairwise Weber attraction (both move equally toward y = 0)\n• Replaces wire-to-wire B-field and Ampère force law',
    antiparallel_wires: '• Two 120 m wires at y = ±2 m, opp-direction currents (|I1| = |I2|)\n• Direct pairwise Weber repulsion (both move equally apart)\n• Replaces wire-to-wire B-field and Ampère force law',
    magnets_attract: '• Two 3-loop bar magnets (R = 1.4 m) at x = ±3.13 m, aligned (N–S)\n• Direct pairwise Weber attraction between circulating Ampèrian loops\n• Replaces magnetic pole / Biot-Savart dipole attraction',
    magnets_repel: '• Two 3-loop bar magnets at x = ±2.3 m, opp currents (N–N) tilted 6°\n• Direct pairwise Weber torque and repulsion rotates and pushes apart\n• Replaces magnetic dipole torque τ = μ × B and repulsion',
    magnet_coil: '• 3-loop bar magnet (vz = 0.06c) moving through coil (R = 6.5 m)\n• Direct Weber velocity force (vᵣ²/2c², v⊥²/c²) drives coil current\n• Replaces Faraday moving-magnet flux induction E = −dΦ_B/dt',
    wire_magnet: '• 3-loop bar magnet at (0, 0, −2.8 m) below DC wire (Ørsted exp)\n• Direct pairwise Weber torque twists magnet to be perp\n• Replaces wire B-field and compass torque τ = μ × B',
    // spherical_shell: '• Stationary electron inside 500-point charged sphere (R = 7 m)\n• Direct Weber acceleration reaction from shell charges\n• Produces F_shell = -(Qq / 3c²R) a'
};

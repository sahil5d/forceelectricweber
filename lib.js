// note: what is p
//     point, pos, or particle.....
// particle is the {prop_id, prop_p, prop_v, ...}.
// p is a point is a pos is a {x: , y: , z:}

const p = THREE.Vector3;
const vec = THREE.Vector3;

const M_PROT = 1.672622 * 10 ** -27;  // kg
const M_ELEC = 9.109384 * 10 ** -31;  // kg

const KEE = 2.307078 * 10 ** -28;   // N*m^2/C^2, Coulomb constant * (1 electron charge)^2

// const BOHR_X = 5.29177 * 10 ** -11; // m
// const PLANCKH = 6.62607015 * 10 ** -34; // m^2 kg / s
// const ALPHA = 0.00729735256;
const CSPEED = 2.9979 * 10 ** 8; // m/s
const CSQUARED = CSPEED ** 2;
const INV_CSQUARED = 1 / CSQUARED;

// const KEEDIVCALPHAALLSQ = (KEE / (ALPHA * CSPEED)) ** 2; // 1.112*10^-68 kg^2 m^4 / s^2. for massmassfx

const JOULETOEV = 6.241509 * 10 ** 18;

const MAX_PARTICLES = 10000;
const ELEC_DRAW_OFFSET = 0.0002; // m

var points_nuke = null;
var points_elec = null;
const pos_nuke = new Float32Array(MAX_PARTICLES * 3);
const pos_elec = new Float32Array(MAX_PARTICLES * 3);
var circle_tex = null;

// Particle kinds (`prop_kind`) and force-target flag (`prop_iswire`):
//   Target kinds (`prop_iswire: false` — accumulate Weber forces in `calc_forces()`):
//     'free'                : unconstrained test particle (Configs 1–4)
//     'coil_secondary_lead' : lead electron on secondary coil receiving induced tangential force (Configs 5, 10)
//     'line_plll_lead'      : central lead atom of a parallel wire receiving transverse force (Configs 6, 7)
//     'line_tension'        : wire segment receiving longitudinal inter-segment force (Config 12)
//     'magnet'              : bar-magnet segment (active corner coils have `prop_iswire: false`; interior coils have `prop_visual_only: true`)
//   Source-only kinds (`prop_iswire: true` — exert Weber forces as source j, skipped as target i and in F = ma):
//     'line'                : straight DC wire segment with constant electron drift (Configs 1, 11)
//     'line_ac'             : straight AC wire segment with oscillating electrons (Config 3)
//     'line_plll'           : non-lead parallel wire segment that follows 'line_plll_lead' (Configs 6, 7)
//     'coil'                : circular coil segment with steady circulating electrons (Configs 2, 10)
//     'coil_turnon'         : circular coil segment whose electron speed ramps up after turn-on (Config 4)
//     'coil_ac'             : AC primary coil segment in transformer (Config 5)
//     'coil_secondary'      : non-lead secondary coil segment that follows 'coil_secondary_lead' (Configs 5, 10)
function add_particle(opts) {
    var charge = opts.prop_charge !== undefined ? opts.prop_charge : -1;
    var p_obj = Object.assign({
        prop_idx: particles.length, // particles.length grows with each push below, assigning sequential index 0, 1, 2, ...
        prop_kind: 'free',
        prop_iswire: false, // true = source-only particle (exerts force as j, skipped as target i and in F = ma)
        prop_p: new p(0, 0, 0),
        prop_v: new vec(0, 0, 0),
        prop_a: new vec(0, 0, 0),
        prop_a_old: new vec(),
        prop_charge: charge,
        prop_mass: charge > 0 ? M_PROT : M_ELEC,
    }, opts);
    if (opts.prop_p) p_obj.prop_p = opts.prop_p.clone();
    if (opts.prop_v) p_obj.prop_v = opts.prop_v.clone();
    if (opts.prop_a) p_obj.prop_a = opts.prop_a.clone();

    if (p_obj.prop_charge < 0 && p_obj.prop_kind !== 'free') {
        var theta_off = Math.random() * 2 * Math.PI;
        var phi_off = Math.acos(2 * Math.random() - 1);
        var r_off = (0.4 + 0.6 * Math.random()) * ELEC_DRAW_OFFSET;
        p_obj.prop_draw_offset = new vec(
            r_off * Math.sin(phi_off) * Math.cos(theta_off),
            r_off * Math.sin(phi_off) * Math.sin(theta_off),
            r_off * Math.cos(phi_off)
        );
    }

    particles.push(p_obj);
    return p_obj;
}

function wrap_coord(x, L) {
    if (x > L / 2) return x - L;
    if (x < -L / 2) return x + L;
    return x;
}

function get_circle_texture() {
    if (circle_tex) return circle_tex;
    if (typeof document !== 'undefined' && document.createElement) {
        var canvas = document.createElement('canvas');
        canvas.width = 64;
        canvas.height = 64;
        var ctx = canvas.getContext && canvas.getContext('2d');
        if (ctx) {
            ctx.beginPath();
            ctx.arc(32, 32, 28, 0, 2 * Math.PI);
            ctx.fillStyle = '#ffffff';
            ctx.fill();
            circle_tex = new THREE.CanvasTexture(canvas);
            return circle_tex;
        }
    }
    return null;
}

function make_point_cloud(buf, color, tex) {
    var attr = new THREE.InstancedBufferAttribute(buf, 3);
    var mat = new THREE.PointsNodeMaterial({
        color: color,
        size: 5,
        sizeAttenuation: false,
        map: tex,
        alphaTest: 0.5,
        positionNode: new THREE.BufferAttributeNode(attr),
    });
    var pts = new THREE.Sprite(mat);
    pts.attr = attr;
    pts.frustumCulled = false;
    scene.add(pts);
    return pts;
}

// in reality estimated electron r is ~ 3 fm, .0003 * 10^-11m
// in reality estimated proton r is ~ 1 fm, .0001 * 10^-11m
function init_plot() {
    var tex = get_circle_texture();
    points_nuke = make_point_cloud(pos_nuke, '#7f2e0c', tex);
    points_elec = make_point_cloud(pos_elec, '#fd4f0b', tex);
}

function plotallparticles() {
    if (!points_nuke) init_plot();

    for (var mag of mag_states) {
        sync_bar_magnet(mag, true);
    }

    var n_nuke = 0;
    var n_elec_plot = 0;

    for (var i = 0; i < particles.length; i++) {
        var particle = particles[i];
        if (particle.prop_kind === 'line_plll_lead') {
            continue;
        }

        if (particle.prop_charge >= 1) {
            pos_nuke[n_nuke * 3] = particle.prop_p.x;
            pos_nuke[n_nuke * 3 + 1] = particle.prop_p.y;
            pos_nuke[n_nuke * 3 + 2] = particle.prop_p.z;
            n_nuke++;
        } else {
            if (particle.prop_visual_only && particle.prop_kind === 'coil') {
                particle.prop_draw_x = (particle.prop_xc || 0) + particle.prop_r * Math.cos(particle.prop_phi);
                particle.prop_draw_y = (particle.prop_yc || 0) + particle.prop_r * Math.sin(particle.prop_phi);
            }
            // prop_p holds the neutral segment physics position; prop_draw_x/y/z holds the drifted visual electron position
            var doff = particle.prop_draw_offset;
            var draw_x = particle.prop_draw_x !== undefined ? particle.prop_draw_x : particle.prop_p.x;
            var draw_y = particle.prop_draw_y !== undefined ? particle.prop_draw_y : particle.prop_p.y;
            var draw_z = particle.prop_draw_z !== undefined ? particle.prop_draw_z : particle.prop_p.z;
            pos_elec[n_elec_plot * 3] = draw_x + (doff ? doff.x : 0);
            pos_elec[n_elec_plot * 3 + 1] = draw_y + (doff ? doff.y : 0);
            pos_elec[n_elec_plot * 3 + 2] = draw_z + (doff ? doff.z : 0);
            n_elec_plot++;
        }
    }

    points_nuke.count = n_nuke;
    points_nuke.attr.needsUpdate = true;

    points_elec.count = n_elec_plot;
    points_elec.attr.needsUpdate = true;
}

function plotlinesegment(argdir, argorigin, arglength, arghex) {
    const points = [argorigin, argorigin.clone().add(argdir.clone().multiplyScalar(arglength))];
    const geometry = new THREE.BufferGeometry();
    geometry.setFromPoints(points);
    const material = new THREE.LineBasicMaterial({ color: arghex, opacity: 0.5, transparent: true });
    scene.add(new THREE.Line(geometry, material));
}

function trim(arg) {
    if (typeof arg === 'number') {
        return parseFloat(arg.toPrecision(3));
    }
    return { x: trim(arg.x), y: trim(arg.y), z: trim(arg.z) };
}

function pct_diff(approx, sim) {
    return trim((approx - sim) / Math.abs(sim) * 100) + '%';
}

function log_force_cmp(label_sim, f_sim, label_approx, f_approx) {
    hud_log(label_sim, trim(f_sim), 'N');
    hud_log(label_approx, trim(f_approx), 'N');
    hud_log('Diff:', pct_diff(f_approx, f_sim));
}

// from https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
function getapproxequidistantpointsonsphere(n, R) {
    var out = [];
    var phiradians = Math.PI * (Math.sqrt(5) - 1); // golden angle in radians
    for (var i = 0; i < n; i++) {
        var y = 1 - 2 * i / (n - 1); // y goes from 1 to -1
        var flatradius = Math.sqrt(1 - y * y); // horizontal radius at y. the radius swinging around yaxis.
        var theta = phiradians * i; // golden angle increment
        out.push(new p(Math.cos(theta) * flatradius * R, y * R, Math.sin(theta) * flatradius * R));
    }
    return out;
}

function set_camera(px, py, pz, tx, ty, tz, fov) {
    if (typeof camera !== 'undefined' && camera) {
        camera.fov = fov || 35;
        if (camera.updateProjectionMatrix) camera.updateProjectionMatrix();
        camera.position.set(px, py, pz);
        camera.lookAt(new THREE.Vector3(tx, ty, tz));
    }
    if (typeof controls !== 'undefined' && controls) {
        controls.target.set(tx, ty, tz);
        controls.update();
    }
}

function hud_clear() {
    if (typeof document !== 'undefined') {
        var el = document.getElementById('hud-log');
        if (el) el.innerHTML = '';
    }
}

function hud_log() {
    var args = Array.prototype.slice.call(arguments);
    var msg = args.map(function (a) {
        return typeof a === 'object' ? JSON.stringify(a) : String(a);
    }).join(' ');
    if (typeof document !== 'undefined') {
        var el = document.getElementById('hud-log');
        if (el) {
            var ph = el.querySelector('.hud-placeholder');
            if (ph) ph.remove();
            var line = document.createElement('div');
            if (args.length >= 2 && typeof args[0] === 'string' && args[0].endsWith(':')) {
                line.className = 'hud-row';
                var key_span = document.createElement('span');
                key_span.className = 'hud-key';
                key_span.textContent = args[0].slice(0, -1);
                var val_span = document.createElement('span');
                val_span.className = 'hud-val';
                val_span.textContent = args.slice(1).map(function (a) {
                    return typeof a === 'object' ? JSON.stringify(a) : String(a);
                }).join(' ');
                line.appendChild(key_span);
                line.appendChild(val_span);
            } else {
                line.className = 'hud-note';
                line.textContent = msg;
            }
            el.appendChild(line);
            return;
        }
    }
    console.log.apply(console, arguments);
}

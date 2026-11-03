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

// const KEEDIVCALPHAALLSQ = (KEE / (ALPHA * CSPEED)) ** 2; // 1.112*10^-68 kg^2 m^4 / s^2. for massmassfx

const JOULETOEV = 6.241509 * 10 ** 18;

const ZEROALL = new p(0, 0, 0);

const MAX_PARTICLES = 10000;
// const TRAIL_LEN = 100;
const ELEC_DRAW_OFFSET = 0.02;

var points_nuke = null;
var points_elec = null;
const pos_nuke = new Float32Array(MAX_PARTICLES * 3);
const pos_elec = new Float32Array(MAX_PARTICLES * 3);
var circle_tex = null;

function add_particle(opts) {
    var charge = opts.prop_charge !== undefined ? opts.prop_charge : -1;
    var p_obj = Object.assign({
        prop_idx: particles.length, // particles.length grows with each push below, assigning sequential index 0, 1, 2, ...
        prop_iswire: false,
        prop_p: new p(0, 0, 0),
        prop_v: new vec(0, 0, 0),
        prop_vtemp: new vec(),
        prop_a: ZEROALL.clone(),
        prop_atemp: new vec(),
        prop_anext: ZEROALL.clone(),
        prop_charge: charge,
        prop_mass: charge > 0 ? M_PROT : M_ELEC,
    }, opts);
    if (opts.prop_p) p_obj.prop_p = opts.prop_p.clone();
    if (opts.prop_v) p_obj.prop_v = opts.prop_v.clone();
    if (opts.prop_a) p_obj.prop_a = opts.prop_a.clone();
    if (opts.prop_anext) p_obj.prop_anext = opts.prop_anext.clone();

    if (p_obj.prop_charge < 0 && (p_obj.prop_iswire || p_obj.prop_is_sec_lead || p_obj.prop_ismagnet || p_obj.prop_istension)) {
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

/*
function record_trail(particle) {
    if (!particle.prop_history) return;
    if (particle.prop_history.length < TRAIL_LEN) {
        particle.prop_history.push(particle.prop_p.clone());
    } else {
        var oldest = particle.prop_history.shift();
        oldest.copy(particle.prop_p);
        particle.prop_history.push(oldest);
    }
}
*/

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
    var geom = new THREE.BufferGeometry();
    geom.setAttribute('position', new THREE.BufferAttribute(buf, 3));
    var pts = new THREE.Points(
        geom,
        new THREE.PointsMaterial({ color: color, size: 6, sizeAttenuation: false, map: tex, alphaTest: 0.5 })
    );
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

    if (typeof mag_states !== 'undefined') {
        for (var m = 0; m < mag_states.length; m++) {
            sync_bar_magnet(mag_states[m], true);
        }
    }

    var n_nuke = 0;
    var n_elec_plot = 0;

    for (var i = 0; i < particles.length; i++) {
        var particle = particles[i];
        if (particle.prop_is_par_lead) {
            continue;
        }

        if (particle.prop_charge >= 1) {
            pos_nuke[n_nuke * 3] = particle.prop_p.x;
            pos_nuke[n_nuke * 3 + 1] = particle.prop_p.y;
            pos_nuke[n_nuke * 3 + 2] = particle.prop_p.z;
            n_nuke++;
        } else {
            var doff = particle.prop_draw_offset;
            pos_elec[n_elec_plot * 3] = particle.prop_p.x + (doff ? doff.x : 0);
            pos_elec[n_elec_plot * 3 + 1] = particle.prop_p.y + (doff ? doff.y : 0);
            pos_elec[n_elec_plot * 3 + 2] = particle.prop_p.z + (doff ? doff.z : 0);
            n_elec_plot++;
        }

        /*
        if (particle.prop_history && particle.prop_history.length > 1) {
            if (!particle.prop_trail_line) {
                var geom_trail = new THREE.BufferGeometry();
                particle.prop_trail_buf = new Float32Array(TRAIL_LEN * 3);
                geom_trail.setAttribute('position', new THREE.BufferAttribute(particle.prop_trail_buf, 3));
                var mat_trail = new THREE.LineBasicMaterial({ color: 0xbcb8b4 });
                particle.prop_trail_line = new THREE.Line(geom_trail, mat_trail);
                particle.prop_trail_line.frustumCulled = false;
                scene.add(particle.prop_trail_line);
            }
            var hist = particle.prop_history;
            var n_pts = hist.length;
            var buf = particle.prop_trail_buf;
            var t_off = particle.prop_draw_offset;
            for (var h = 0; h < n_pts; h++) {
                buf[h * 3] = hist[h].x + (t_off ? t_off.x : 0);
                buf[h * 3 + 1] = hist[h].y + (t_off ? t_off.y : 0);
                buf[h * 3 + 2] = hist[h].z + (t_off ? t_off.z : 0);
            }
            particle.prop_trail_line.geometry.setDrawRange(0, n_pts);
            particle.prop_trail_line.geometry.attributes.position.needsUpdate = true;
        }
        */
    }

    points_nuke.geometry.setDrawRange(0, n_nuke);
    points_nuke.geometry.attributes.position.needsUpdate = true;

    points_elec.geometry.setDrawRange(0, n_elec_plot);
    points_elec.geometry.attributes.position.needsUpdate = true;
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
    var msg = Array.prototype.slice.call(arguments).map(function (a) {
        return typeof a === 'object' ? JSON.stringify(a) : String(a);
    }).join(' ');
    if (typeof document !== 'undefined') {
        var el = document.getElementById('hud-log');
        if (el) {
            var ph = el.querySelector('.hud-placeholder');
            if (ph) ph.remove();
            var line = document.createElement('div');
            line.className = 'hud-line';
            line.textContent = msg;
            el.appendChild(line);
            return;
        }
    }
    console.log.apply(console, arguments);
}

// Accumulate pairwise Weber forces directly onto each particle's prop_f
var dir = new vec();
var pdiff = new vec();
var vdiff = new vec();
var adiff = new vec();
var force = new vec(); // the initial -1 is so -1*oppchargeproduct will be 1, -1*samechargeproduct will be -1.
var minr; // save smallest pairwise distance

// Flat Structure-of-Arrays (SoA) buffers for fast unboxed pairwise force evaluation
const _px = new Float64Array(MAX_PARTICLES);
const _py = new Float64Array(MAX_PARTICLES);
const _pz = new Float64Array(MAX_PARTICLES);
const _vx = new Float64Array(MAX_PARTICLES);
const _vy = new Float64Array(MAX_PARTICLES);
const _vz = new Float64Array(MAX_PARTICLES);
const _ax = new Float64Array(MAX_PARTICLES);
const _ay = new Float64Array(MAX_PARTICLES);
const _az = new Float64Array(MAX_PARTICLES);
const _fx = new Float64Array(MAX_PARTICLES);
const _fy = new Float64Array(MAX_PARTICLES);
const _fz = new Float64Array(MAX_PARTICLES);
const _q = new Float64Array(MAX_PARTICLES);
const _iswire = new Uint8Array(MAX_PARTICLES);
const _jstart = new Int32Array(MAX_PARTICLES);
const _jend = new Int32Array(MAX_PARTICLES);

// call this before loop, and every loop
function calc_forces() {
    var local_minr = 1000; // 1000 meters
    var pts = particles;
    var n = pts.length;
    var px = _px, py = _py, pz = _pz;
    var vx = _vx, vy = _vy, vz = _vz;
    var ax = _ax, ay = _ay, az = _az;
    var fx = _fx, fy = _fy, fz = _fz;
    var q = _q, iswire = _iswire, jstart = _jstart, jend = _jend;

    for (var i = 0; i < n; i++) {
        var pt = pts[i];
        var p_i = pt.prop_p, v_i = pt.prop_v, a_i = pt.prop_a;
        px[i] = p_i.x; py[i] = p_i.y; pz[i] = p_i.z;
        vx[i] = v_i.x; vy[i] = v_i.y; vz[i] = v_i.z;
        ax[i] = a_i.x; ay[i] = a_i.y; az[i] = a_i.z;
        fx[i] = 0; fy[i] = 0; fz[i] = 0;
        q[i] = pt.prop_charge;
        iswire[i] = pt.prop_iswire ? 1 : 0;
    }

    // Precompute contiguous interaction range [jstart[i], jend[i]) at step 0 to skip same-group / wire particles
    if (steps === 0) {
        for (var i = 0; i < n; i++) {
            if (iswire[i]) {
                jstart[i] = n;
                jend[i] = n;
                continue;
            }
            var pi = pts[i];
            var pi_group = pi.prop_group;
            var pi_par = pi.prop_is_par_lead;
            var js = i + 1;
            // Ignore internal forces within the same wire, coil, or bar magnet
            while (js < n && ((pi_group && pts[js].prop_group === pi_group) || (pi_par && pts[js].prop_is_par_lead))) {
                js++;
            }
            var je = js;
            while (je < n && !(pi_group && pts[je].prop_group === pi_group) && !(pi_par && pts[je].prop_is_par_lead)) {
                je++;
            }
            jstart[i] = js;
            jend[i] = je;
        }
    }

    var inv_c2 = 1 / CSQUARED;
    var half_inv_c2 = 0.5 * inv_c2;
    var kee = KEE;

    // note: could change iloop max from '< .length-1' to '< .length' if need to do a check on the last element in array (like max r from origin). causes no issue on pairdata (inner jloop won't run when i=.length-1)

    for (var i = 0; i < n - 1; i++) {
        var js = jstart[i];
        var je = jend[i];
        if (js >= je) {
            continue;
        }

        var qi_kee = -kee * q[i];
        var pix = px[i], piy = py[i], piz = pz[i];
        var vix = vx[i], viy = vy[i], viz = vz[i];
        var aix = ax[i], aiy = ay[i], aiz = az[i];
        var fix = 0, fiy = 0, fiz = 0;

        for (var j = js; j < je; j++) {
            var dx = px[j] - pix;
            var dy = py[j] - piy;
            var dz = pz[j] - piz;
            var distsq = dx * dx + dy * dy + dz * dz; // better than calling .distanceToSquared() which repeats .sub operation
            var dist = Math.sqrt(distsq);
            var inv_dist = 1 / dist;

            // kqq/r^2 Coulomb
            // Weber force:
            //   F = (kq1q2 / r^2) * (1 - rdot^2 / (2 * c^2) + r * rddot / c^2)
            //     = (kq1q2 / r^2) * (1 - v_rad^2 / (2 * c^2) + v_perp^2 / c^2 + r * a_r / c^2)
            // easy definitions:
            //
            //      v is velocity vector
            //      v_radial is radial velocity, vector, https://en.wikipedia.org/wiki/Radial_velocity
            //      rdot is runit dot v, scalar
            //
            //      a is acceleration vector
            //      a_radial is radial acceleration, vector
            //      runit dot a is scalar radial acceleration
            //      rdoubledot is change in rdot, is a scalar, it is not the same as radial acceleration
            //          = v_tangent^2 / r    + (r_unit dot a)
            //          = (v^2 - rdot^2) / r + (r_unit dot a)
            //          case of object flyby constant speed by stationary object
            //              acceleration is 0, so radial acceleration runit dot a is 0
            //              but during the flyby the rdot increases (- to 0 to +) so the rdoubledot is positive
            //          case of object going in circle around stationary object
            //              acceleration is -v_tangent^2/r
            //              but rdot=0, change in rdot=0, so rdoubledot=0
            //              why? because differentiating the dist twice is not same as differentiating vector r twice

            var dvx = vx[j] - vix;
            var dvy = vy[j] - viy;
            var dvz = vz[j] - viz;
            var dax = ax[j] - aix;
            var day = ay[j] - aiy;
            var daz = az[j] - aiz;

            var rdot = (dx * dvx + dy * dvy + dz * dvz) * inv_dist;
            var vdiff_sq = dvx * dvx + dvy * dvy + dvz * dvz;
            var rddot = (vdiff_sq - rdot * rdot + dx * dax + dy * day + dz * daz) * inv_dist;
            var weberfactor = 1 - rdot * rdot * half_inv_c2 + dist * rddot * inv_c2;
            var fscale = qi_kee * q[j] * (inv_dist / distsq) * weberfactor;

            var dfx = dx * fscale;
            var dfy = dy * fscale;
            var dfz = dz * fscale;

            fix += dfx;
            fiy += dfy;
            fiz += dfz;
            fx[j] -= dfx;
            fy[j] -= dfy;
            fz[j] -= dfz;

            // optional
            if (dist < local_minr) {
                local_minr = dist;
            }
        }

        fx[i] += fix;
        fy[i] += fiy;
        fz[i] += fiz;
    }

    minr = local_minr;
    for (var i = 0; i < n; i++) {
        if (!iswire[i]) {
            pts[i].prop_f.set(fx[i], fy[i], fz[i]);
        }
    }
}

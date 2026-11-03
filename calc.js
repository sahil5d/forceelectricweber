// Accumulate pairwise Weber forces into flat Structure-of-Arrays (SoA) buffers (_fx, _fy, _fz)

// Flat Structure-of-Arrays (SoA) buffers for fast unboxed pairwise force evaluation
const _px = new Float64Array(MAX_PARTICLES), _py = new Float64Array(MAX_PARTICLES), _pz = new Float64Array(MAX_PARTICLES);
const _vx = new Float64Array(MAX_PARTICLES), _vy = new Float64Array(MAX_PARTICLES), _vz = new Float64Array(MAX_PARTICLES);
const _ax = new Float64Array(MAX_PARTICLES), _ay = new Float64Array(MAX_PARTICLES), _az = new Float64Array(MAX_PARTICLES);
const _fx = new Float64Array(MAX_PARTICLES), _fy = new Float64Array(MAX_PARTICLES), _fz = new Float64Array(MAX_PARTICLES);
const _q = new Float64Array(MAX_PARTICLES), _iswire = new Uint8Array(MAX_PARTICLES);
const _jstart = new Int32Array(MAX_PARTICLES), _jend = new Int32Array(MAX_PARTICLES);

// call this before loop, and every loop
function calc_forces() {
    // Local aliases keep SoA typed arrays and constants in function stack slots for the O(N^2) inner loop
    var pts = particles;
    var n = pts.length;
    var px = _px, py = _py, pz = _pz;
    var vx = _vx, vy = _vy, vz = _vz;
    var ax = _ax, ay = _ay, az = _az;
    var fx = _fx, fy = _fy, fz = _fz;
    var q = _q, iswire = _iswire, jstart = _jstart, jend = _jend;

    for (var i = 0; i < n; i++) {
        var pt = pts[i];
        if (pt.prop_visual_only) {
            iswire[i] = 1;
            continue;
        }
        var p_i = pt.prop_p, v_i = pt.prop_v, a_i = pt.prop_a;
        px[i] = p_i.x; py[i] = p_i.y; pz[i] = p_i.z;
        vx[i] = v_i.x; vy[i] = v_i.y; vz[i] = v_i.z;
        ax[i] = a_i.x; ay[i] = a_i.y; az[i] = a_i.z;
        fx[i] = 0; fy[i] = 0; fz[i] = 0;
        q[i] = pt.prop_charge;
        iswire[i] = pt.prop_iswire ? 1 : 0;
    }

    // Precompute contiguous interaction range [jstart[i], jend[i]) at step 0 to skip same-group / wire / visual-only particles
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
            // Ignore internal forces within the same wire (line), coil, or bar magnet, and skip visual-only coils
            while (js < n && (pts[js].prop_visual_only || (pi_group && pts[js].prop_group === pi_group) || (pi_par && pts[js].prop_is_par_lead))) {
                js++;
            }
            var je = js;
            while (je < n && !pts[je].prop_visual_only && !(pi_group && pts[je].prop_group === pi_group) && !(pi_par && pts[je].prop_is_par_lead)) {
                je++;
            }
            jstart[i] = js;
            jend[i] = je;
        }
    }

    var inv_c2 = INV_CSQUARED;
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

            // kqq/r^2 Coulomb
            // Weber force:
            //   F = (kq1q2 / r^2) * (1 - rdot^2 / (2 * c^2) + r rddot / c^2)
            //     = (kq1q2 / r^2) * (1 - v_rad^2 / (2 * c^2) + v_perp^2 / c^2 + r a_r / c^2)
            //     = (kq1q2 / r^2) * (1 + (v^2 - 1.5 v_rad^2 + r dot a) / c^2)
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

            var r_dot_v = dx * dvx + dy * dvy + dz * dvz;
            var vr_sq = (r_dot_v * r_dot_v) / distsq;
            var v_sq = dvx * dvx + dvy * dvy + dvz * dvz;
            var r_ar = dx * dax + dy * day + dz * daz;
            var weberfactor = 1 + (v_sq - 1.5 * vr_sq + r_ar) * inv_c2;
            var fmag = (qi_kee * q[j] / distsq) * weberfactor;

            var dfx = fmag * (dx / dist);
            var dfy = fmag * (dy / dist);
            var dfz = fmag * (dz / dist);

            fix += dfx;
            fiy += dfy;
            fiz += dfz;
            fx[j] -= dfx;
            fy[j] -= dfy;
            fz[j] -= dfz;
        }

        fx[i] += fix;
        fy[i] += fiy;
        fz[i] += fiz;
    }
}

// calc accel for each particle
function calc_each_accel() {
    // add +or- 1/r^2 vectors
    for (var i = 0; i < particles.length; i++) {
        var p_i = particles[i];
        if (p_i.prop_iswire || p_i.prop_ismagnet) {
            continue;
        }

        if (p_i.prop_istension || p_i.prop_is_par_lead) {
            if (i === p_i.prop_lead_idx) {
                var ax = p_i.prop_istension ? 0.5 * (_fx[i] + _fx[i + 1]) / p_i.prop_mass : 0;
                var ay = p_i.prop_is_par_lead ? 0.5 * (_fy[i] + _fy[i + 1]) / p_i.prop_mass : 0;
                p_i.prop_a_old.copy(p_i.prop_a);
                particles[i + 1].prop_a_old.copy(particles[i + 1].prop_a);
                p_i.prop_a.set(ax, ay, 0);
                particles[i + 1].prop_a.set(ax, ay, 0);
            }
            continue;
        }

        p_i.prop_a_old.copy(p_i.prop_a);
        var inv_m = 1 / p_i.prop_mass;
        p_i.prop_a.set(_fx[i] * inv_m, _fy[i] * inv_m, _fz[i] * inv_m);
    }

    for (var mag of mag_states) {
        var fx = 0, fy = 0, tau = 0;
        for (var pt of mag.particles) {
            var idx = pt.prop_idx;
            var fpx = _fx[idx], fpy = _fy[idx];
            fx += fpx;
            fy += fpy;
            tau += pt.prop_dx * fpy - pt.prop_dy * fpx;
        }
        mag.fx = fx;
        mag.fy = fy;
        mag.tau = tau;
        mag.ax_old = mag.ax;
        mag.ay_old = mag.ay;
        mag.alpha_old = mag.alpha;
        mag.ax = fx / mag.mass;
        mag.ay = fy / mag.mass;
        mag.alpha = tau / mag.inertia;
    }
}

function tangential_comp(v, phi) {
    return -v.x * Math.sin(phi) + v.y * Math.cos(phi);
}

function set_coil_state(p_obj, phi, vd, ad) {
    var R = p_obj.prop_r, cosP = Math.cos(phi), sinP = Math.sin(phi);
    p_obj.prop_phi = phi;
    p_obj.prop_p.x = (p_obj.prop_xc || 0) + R * cosP;
    p_obj.prop_p.y = (p_obj.prop_yc || 0) + R * sinP;
    p_obj.prop_v.x = -vd * sinP;
    p_obj.prop_v.y = vd * cosP;
    if (ad !== undefined) {
        p_obj.prop_a.x = -ad * sinP;
        p_obj.prop_a.y = ad * cosP;
    }
}

function step_coil_electron(p_obj, vd, ad) {
    set_coil_state(p_obj, p_obj.prop_phi + (vd / p_obj.prop_r) * timestep, vd, ad);
}

function sync_bar_magnet(mag, sync_visual) {
    var cosT = Math.cos(mag.theta), sinT = Math.sin(mag.theta);
    for (var pt of mag.particles) {
        var cosP = Math.cos(pt.prop_phi), sinP = Math.sin(pt.prop_phi);
        var u = pt.prop_uc + MAG_R * cosP;
        var w = pt.prop_wc + MAG_R * sinP;
        var dx = pt.prop_s * cosT - u * sinT;
        var dy = pt.prop_s * sinT + u * cosT;
        pt.prop_dx = dx;
        pt.prop_dy = dy;
        pt.prop_p.set(mag.x + dx, mag.y + dy, mag.z + w);
        var vrx = mag.vx - mag.omega * dy, vry = mag.vy + mag.omega * dx;
        if (pt.prop_charge > 0) {
            pt.prop_v.set(vrx, vry, 0);
        } else {
            pt.prop_v.set(vrx + mag.vd * sinP * sinT, vry - mag.vd * sinP * cosT, mag.vd * cosP);
        }
    }
    if (sync_visual) {
        var phi_shift = mag.phi_shift;
        for (var vpt of mag.visual_particles) {
            var phi = vpt.prop_charge < 0 ? vpt.prop_phi + phi_shift : vpt.prop_phi;
            var u_v = vpt.prop_uc + MAG_R * Math.cos(phi);
            var w_v = vpt.prop_wc + MAG_R * Math.sin(phi);
            vpt.prop_p.set(
                mag.x + vpt.prop_s * cosT - u_v * sinT,
                mag.y + vpt.prop_s * sinT + u_v * cosT,
                mag.z + w_v
            );
        }
    }
}

function calc_each_vel_and_accel() {
    calc_each_accel();

    for (var mag of mag_states) {
        mag.vx += 0.5 * (mag.ax_old + mag.ax) * timestep;
        mag.vy += 0.5 * (mag.ay_old + mag.ay) * timestep;
        mag.omega += 0.5 * (mag.alpha_old + mag.alpha) * timestep;
        sync_bar_magnet(mag);
    }

    for (var i = 0; i < particles.length; i++) {
        var p_i = particles[i];
        if (p_i.prop_iswire || p_i.prop_ismagnet) {
            continue;
        }
        if (p_i.prop_is_sec_lead) {
            var phi_m = p_i.prop_phi;
            var aphi_prev = tangential_comp(p_i.prop_a_old, phi_m);
            var aphi_next = tangential_comp(p_i.prop_a, phi_m);
            p_i.prop_vphi += 0.5 * (aphi_prev + aphi_next) * timestep;
            p_i.prop_v.set(-p_i.prop_vphi * Math.sin(phi_m), p_i.prop_vphi * Math.cos(phi_m), 0);
            continue;
        }
        p_i.prop_v.x += 0.5 * (p_i.prop_a_old.x + p_i.prop_a.x) * timestep;
        p_i.prop_v.y += 0.5 * (p_i.prop_a_old.y + p_i.prop_a.y) * timestep;
        p_i.prop_v.z += 0.5 * (p_i.prop_a_old.z + p_i.prop_a.z) * timestep;
    }
}

function calc_each_pos() {
    var omega_ac = (2 * Math.PI * AC_CYCLES) / (stepsmax * timestep);
    var omega_tfr = (2 * Math.PI * TFR_CYCLES) / (stepsmax * timestep);
    var t_next = totaltime + timestep;

    var coil_vd_curr = 0;
    var coil_ad_curr = 0;
    if (steps >= COIL_TURNON_STEP) {
        var dt_on = t_next - COIL_TURNON_STEP * timestep;
        coil_vd_curr = COIL_TURNON_VD * (1 - Math.exp(-dt_on / COIL_TURNON_TAU));
        coil_ad_curr = (COIL_TURNON_VD / COIL_TURNON_TAU) * Math.exp(-dt_on / COIL_TURNON_TAU);
    }

    for (var mag of mag_states) {
        mag.x += mag.vx * timestep + 0.5 * mag.ax * (timestep ** 2);
        mag.y += mag.vy * timestep + 0.5 * mag.ay * (timestep ** 2);
        mag.theta += mag.omega * timestep + 0.5 * mag.alpha * (timestep ** 2);
        var dphi = (mag.vd / MAG_R) * timestep;
        mag.phi_shift += dphi;
        for (var pt of mag.particles) {
            if (pt.prop_charge < 0) {
                pt.prop_phi += dphi;
            }
        }
        sync_bar_magnet(mag);
    }

    for (var i = 0; i < particles.length; i++) {
        var p_i = particles[i];
        if (p_i.prop_ismagnet) {
            continue;
        }

        if (p_i.prop_istension || p_i.prop_is_par_lead) {
            if (i === p_i.prop_lead_idx) {
                p_i.prop_p.x += p_i.prop_v.x * timestep + 0.5 * p_i.prop_a.x * (timestep ** 2);
                p_i.prop_p.y += p_i.prop_v.y * timestep + 0.5 * p_i.prop_a.y * (timestep ** 2);
                particles[i + 1].prop_p.copy(p_i.prop_p);
            }
            continue;
        }

        if (p_i.prop_ispar) {
            var lead = particles[p_i.prop_lead_idx];
            p_i.prop_p.y = lead.prop_p.y;
            p_i.prop_v.y = lead.prop_v.y;
            if (p_i.prop_charge < 0) {
                p_i.prop_p.x = wrap_coord(p_i.prop_p.x + p_i.prop_vd * timestep, PAR_L);
            }
            continue;
        }

        if (p_i.prop_isprimary) {
            if (p_i.prop_charge < 0) {
                var phi_p = p_i.prop_phi0 - TFR_AMP_PHI1 * Math.cos(omega_tfr * t_next);
                var vphi_p = p_i.prop_r * TFR_AMP_PHI1 * omega_tfr * Math.sin(omega_tfr * t_next);
                var aphi_p = p_i.prop_r * TFR_AMP_PHI1 * (omega_tfr ** 2) * Math.cos(omega_tfr * t_next);
                set_coil_state(p_i, phi_p, vphi_p, aphi_p);
            }
            continue;
        }

        if (p_i.prop_issecondary) {
            if (p_i.prop_is_sec_lead) {
                var phi_old = p_i.prop_phi;
                var aphi_curr = tangential_comp(p_i.prop_a, phi_old);
                var phi_new = phi_old + (p_i.prop_vphi * timestep + 0.5 * aphi_curr * (timestep ** 2)) / p_i.prop_r;
                p_i.prop_phi = phi_new;
                p_i.prop_p.x = p_i.prop_r * Math.cos(phi_new);
                p_i.prop_p.y = p_i.prop_r * Math.sin(phi_new);
                // record_trail(p_i);
            } else if (p_i.prop_charge < 0) {
                set_coil_state(p_i, p_i.prop_phi0 + particles[0].prop_phi, particles[0].prop_vphi);
            }
            continue;
        }

        if (p_i.prop_iscoil_turnon) {
            if (p_i.prop_charge < 0) {
                step_coil_electron(p_i, coil_vd_curr, coil_ad_curr);
            }
            continue;
        }

        if (p_i.prop_isacwire) {
            if (p_i.prop_charge < 0) {
                p_i.prop_p.x = wrap_coord(p_i.prop_x0 - AC_AMP * Math.cos(omega_ac * t_next), AC_WIRE_L);
                p_i.prop_v.x = AC_AMP * omega_ac * Math.sin(omega_ac * t_next);
                p_i.prop_a.x = AC_AMP * (omega_ac ** 2) * Math.cos(omega_ac * t_next);
            }
            continue;
        }

        if (p_i.prop_iscoil) {
            if (p_i.prop_vz) {
                p_i.prop_p.z += p_i.prop_vz * timestep;
            }
            if (p_i.prop_charge < 0) {
                if (p_i.prop_visual_only) {
                    p_i.prop_phi += (p_i.prop_vd / p_i.prop_r) * timestep;
                } else {
                    step_coil_electron(p_i, p_i.prop_vd);
                }
            }
            continue;
        }

        if (p_i.prop_isdcwire) {
            if (p_i.prop_charge < 0) {
                p_i.prop_p.x = wrap_coord(p_i.prop_p.x + p_i.prop_vd * timestep, p_i.prop_wire_l);
            }
            continue;
        }

        if (p_i.prop_iswire) {
            continue;
        }

        p_i.prop_p.x += p_i.prop_v.x * timestep + 0.5 * p_i.prop_a.x * (timestep ** 2);
        p_i.prop_p.y += p_i.prop_v.y * timestep + 0.5 * p_i.prop_a.y * (timestep ** 2);
        p_i.prop_p.z += p_i.prop_v.z * timestep + 0.5 * p_i.prop_a.z * (timestep ** 2);

        // record_trail(p_i);
    }
}

// Accumulate pairwise Weber forces into flat Structure-of-Arrays (SoA) buffers (_fx, _fy, _fz)

// Flat Structure-of-Arrays (SoA) buffers for fast unboxed pairwise force evaluation
const _px = new Float64Array(MAX_PARTICLES), _py = new Float64Array(MAX_PARTICLES), _pz = new Float64Array(MAX_PARTICLES);
const _vx = new Float64Array(MAX_PARTICLES), _vy = new Float64Array(MAX_PARTICLES), _vz = new Float64Array(MAX_PARTICLES);
const _ax = new Float64Array(MAX_PARTICLES), _ay = new Float64Array(MAX_PARTICLES), _az = new Float64Array(MAX_PARTICLES);
const _fx = new Float64Array(MAX_PARTICLES), _fy = new Float64Array(MAX_PARTICLES), _fz = new Float64Array(MAX_PARTICLES);
const _q = new Float64Array(MAX_PARTICLES), _iswire = new Uint8Array(MAX_PARTICLES);
const _group = new Int32Array(MAX_PARTICLES), _jstart = new Int32Array(MAX_PARTICLES), _jend = new Int32Array(MAX_PARTICLES);

// call this before loop, and every loop
function calc_forces() {
    // Local aliases keep SoA typed arrays and constants in function stack slots for the O(N^2) inner loop
    var pts = particles;
    var n = pts.length;
    var px = _px, py = _py, pz = _pz;
    var vx = _vx, vy = _vy, vz = _vz;
    var ax = _ax, ay = _ay, az = _az;
    var fx = _fx, fy = _fy, fz = _fz;
    var q = _q, iswire = _iswire, grp = _group, jstart = _jstart, jend = _jend;

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

    // Precompute group IDs and outer interaction bounds [jstart[i], jend[i]) at step 0.
    // Note: Do NOT assume valid target particles j > i form a single contiguous range in `particles`.
    // [jstart[i], jend[i]) only trims leading/trailing excluded particles; any same-group or visual-only
    // particles in the middle of [jstart[i], jend[i]) are skipped in the inner j-loop via `grp[j]`,
    // so configs can push wires, coils, and magnets into `particles` in any order.
    if (steps === 0) {
        var group_map = {};
        var next_gid = 1;
        for (var i = 0; i < n; i++) {
            var pt = pts[i];
            if (pt.prop_visual_only || pt.prop_kind === 'line_plll_lead') {
                grp[i] = -1;
            } else if (pt.prop_group) {
                var gname = pt.prop_group;
                if (!group_map[gname]) group_map[gname] = next_gid++;
                grp[i] = group_map[gname];
            } else {
                grp[i] = 0;
            }
        }
        for (var i = 0; i < n; i++) {
            if (iswire[i]) {
                jstart[i] = n;
                jend[i] = n;
                continue;
            }
            var gi = grp[i] > 0 ? grp[i] : (pts[i].prop_group ? (group_map[pts[i].prop_group] || 0) : 0);
            var js = i + 1;
            while (js < n && (grp[js] < 0 || (gi !== 0 && grp[js] === gi))) {
                js++;
            }
            var je = n;
            while (je > js && (grp[je - 1] < 0 || (gi !== 0 && grp[je - 1] === gi))) {
                je--;
            }
            // Store active particle's own positive group ID so the inner loop can skip any non-contiguous same-group blocks
            grp[i] = gi;
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

        var gi = grp[i];
        var qi_kee = -kee * q[i];
        var pix = px[i], piy = py[i], piz = pz[i];
        var vix = vx[i], viy = vy[i], viz = vz[i];
        var aix = ax[i], aiy = ay[i], aiz = az[i];
        var fix = 0, fiy = 0, fiz = 0;

        for (var j = js; j < je; j++) {
            var gj = grp[j];
            if (gj < 0 || (gi !== 0 && gj === gi)) {
                continue;
            }
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
        if (p_i.prop_iswire || p_i.prop_kind === 'magnet') {
            continue;
        }

        if (p_i.prop_kind === 'line_tension' || p_i.prop_kind === 'line_plll_lead') {
            if (i === p_i.prop_lead_idx) {
                var ax = p_i.prop_kind === 'line_tension' ? 0.5 * (_fx[i] + _fx[i + 1]) / p_i.prop_mass : 0;
                var ay = p_i.prop_kind === 'line_plll_lead' ? 0.5 * (_fy[i] + _fy[i + 1]) / p_i.prop_mass : 0;
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

function set_coil_state(p_obj, phi, vd, ad) {
    var R = p_obj.prop_r;
    var phi_phys = p_obj.prop_phi0 !== undefined ? p_obj.prop_phi0 : phi;
    var cosP = Math.cos(phi_phys), sinP = Math.sin(phi_phys);
    p_obj.prop_phi = phi;
    p_obj.prop_draw_x = (p_obj.prop_xc || 0) + R * Math.cos(phi);
    p_obj.prop_draw_y = (p_obj.prop_yc || 0) + R * Math.sin(phi);
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
    var R = CFG.magnets_attract.R;
    var cosT = Math.cos(mag.theta), sinT = Math.sin(mag.theta);
    for (var pt of mag.particles) {
        var cosP = Math.cos(pt.prop_phi), sinP = Math.sin(pt.prop_phi);
        var u = pt.prop_uc + R * cosP;
        var w = pt.prop_wc + R * sinP;
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
        for (var vpt of mag.all_particles) {
            if (vpt.prop_charge > 0 && !vpt.prop_visual_only) continue;
            var phi = vpt.prop_charge < 0 ? vpt.prop_phi + mag.phi_shift : vpt.prop_phi;
            var u_v = vpt.prop_uc + R * Math.cos(phi);
            var w_v = vpt.prop_wc + R * Math.sin(phi);
            var x = mag.x + vpt.prop_s * cosT - u_v * sinT;
            var y = mag.y + vpt.prop_s * sinT + u_v * cosT;
            var z = mag.z + w_v;
            if (vpt.prop_charge < 0) {
                vpt.prop_draw_x = x;
                vpt.prop_draw_y = y;
                vpt.prop_draw_z = z;
            } else {
                vpt.prop_p.set(x, y, z);
            }
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
        if (p_i.prop_iswire || p_i.prop_kind === 'magnet') {
            continue;
        }
        if (p_i.prop_kind === 'coil_secondary_lead') {
            p_i.prop_vphi += 0.5 * (p_i.prop_a_old.y + p_i.prop_a.y) * timestep;
            p_i.prop_v.set(0, p_i.prop_vphi, 0);
            continue;
        }
        p_i.prop_v.x += 0.5 * (p_i.prop_a_old.x + p_i.prop_a.x) * timestep;
        p_i.prop_v.y += 0.5 * (p_i.prop_a_old.y + p_i.prop_a.y) * timestep;
        p_i.prop_v.z += 0.5 * (p_i.prop_a_old.z + p_i.prop_a.z) * timestep;
    }
}

function calc_each_pos() {
    var omega_ac = (2 * Math.PI * CFG.line_ac.CYCLES) / (stepsmax * timestep);
    var omega_tfr = (2 * Math.PI * CFG.coils_ac.CYCLES) / (stepsmax * timestep);
    var t_next = totaltime + timestep;

    var coil_vd_curr = 0;
    var coil_ad_curr = 0;
    if (steps >= CFG.coil_turnon.STEP) {
        var dt_on = t_next - CFG.coil_turnon.STEP * timestep;
        coil_vd_curr = CFG.coil_turnon.VD * (1 - Math.exp(-dt_on / CFG.coil_turnon.TAU));
        coil_ad_curr = (CFG.coil_turnon.VD / CFG.coil_turnon.TAU) * Math.exp(-dt_on / CFG.coil_turnon.TAU);
    }

    for (var mag of mag_states) {
        mag.x += mag.vx * timestep + 0.5 * mag.ax * (timestep ** 2);
        mag.y += mag.vy * timestep + 0.5 * mag.ay * (timestep ** 2);
        mag.theta += mag.omega * timestep + 0.5 * mag.alpha * (timestep ** 2);
        mag.phi_shift += (mag.vd / CFG.magnets_attract.R) * timestep;
        sync_bar_magnet(mag);
    }

    for (var i = 0; i < particles.length; i++) {
        var p_i = particles[i];
        switch (p_i.prop_kind) {
            case 'magnet':
                break;

            case 'line_tension':
            case 'line_plll_lead':
                if (i === p_i.prop_lead_idx) {
                    p_i.prop_p.x += p_i.prop_v.x * timestep + 0.5 * p_i.prop_a.x * (timestep ** 2);
                    p_i.prop_p.y += p_i.prop_v.y * timestep + 0.5 * p_i.prop_a.y * (timestep ** 2);
                    particles[i + 1].prop_p.copy(p_i.prop_p);
                }
                break;

            case 'line_plll': {
                var lead = particles[p_i.prop_lead_idx];
                p_i.prop_p.y = lead.prop_p.y;
                p_i.prop_v.y = lead.prop_v.y;
                if (p_i.prop_charge < 0) {
                    p_i.prop_draw_x = wrap_coord(p_i.prop_draw_x + p_i.prop_vd * timestep, CFG.lines_parallel.L);
                }
                break;
            }

            case 'coil_ac':
                if (p_i.prop_charge < 0) {
                    var phi_p = p_i.prop_phi0 - CFG.coils_ac.AMP_PHI1 * Math.cos(omega_tfr * t_next);
                    var vphi_p = p_i.prop_r * CFG.coils_ac.AMP_PHI1 * omega_tfr * Math.sin(omega_tfr * t_next);
                    var aphi_p = p_i.prop_r * CFG.coils_ac.AMP_PHI1 * (omega_tfr ** 2) * Math.cos(omega_tfr * t_next);
                    set_coil_state(p_i, phi_p, vphi_p, aphi_p);
                }
                break;

            case 'coil_secondary_lead': {
                var aphi_curr = p_i.prop_a.y;
                var phi_new = p_i.prop_phi + (p_i.prop_vphi * timestep + 0.5 * aphi_curr * (timestep ** 2)) / p_i.prop_r;
                p_i.prop_phi = phi_new;
                p_i.prop_draw_x = p_i.prop_r * Math.cos(phi_new);
                p_i.prop_draw_y = p_i.prop_r * Math.sin(phi_new);
                break;
            }

            case 'coil_secondary':
                if (p_i.prop_charge < 0) {
                    set_coil_state(p_i, p_i.prop_phi0 + particles[0].prop_phi, particles[0].prop_vphi);
                }
                break;

            case 'coil_turnon':
                if (p_i.prop_charge < 0) {
                    step_coil_electron(p_i, coil_vd_curr, coil_ad_curr);
                }
                break;

            case 'line_ac':
                if (p_i.prop_charge < 0) {
                    p_i.prop_draw_x = wrap_coord(p_i.prop_x0 - CFG.line_ac.AMP * Math.cos(omega_ac * t_next), CFG.line_ac.L);
                    p_i.prop_v.x = CFG.line_ac.AMP * omega_ac * Math.sin(omega_ac * t_next);
                    p_i.prop_a.x = CFG.line_ac.AMP * (omega_ac ** 2) * Math.cos(omega_ac * t_next);
                }
                break;

            case 'coil':
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
                break;

            case 'line':
                if (p_i.prop_charge < 0) {
                    p_i.prop_draw_x = wrap_coord(p_i.prop_draw_x + p_i.prop_vd * timestep, p_i.prop_wire_l);
                }
                break;

            case 'free':
                if (!p_i.prop_iswire) {
                    p_i.prop_p.x += p_i.prop_v.x * timestep + 0.5 * p_i.prop_a.x * (timestep ** 2);
                    p_i.prop_p.y += p_i.prop_v.y * timestep + 0.5 * p_i.prop_a.y * (timestep ** 2);
                    p_i.prop_p.z += p_i.prop_v.z * timestep + 0.5 * p_i.prop_a.z * (timestep ** 2);
                }
                break;
        }
    }
}

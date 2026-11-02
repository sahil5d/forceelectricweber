// calc accel for each particle
// or calc accelnext (arg foranext)
function calc_each_accel(foranext) {
    // add +or- 1/r^2 vectors
    for (var i = 0; i < particles.length; i++) {
        var p_i = particles[i];
        if (p_i.prop_iswire || p_i.prop_ismagnet) {
            continue;
        }

        var target_a = foranext ? p_i.prop_anext : p_i.prop_a;
        target_a.copy(p_i.prop_f).multiplyScalar(1 / p_i.prop_mass);
    }

    for (var mag of mag_states) {
        var fx = 0, fy = 0, tau = 0;
        for (var pt of mag.particles) {
            fx += pt.prop_f.x;
            fy += pt.prop_f.y;
            tau += pt.prop_dx * pt.prop_f.y - pt.prop_dy * pt.prop_f.x;
        }
        mag.fx = fx;
        mag.fy = fy;
        mag.tau = tau;
        mag.ax = fx / mag.mass;
        mag.ay = fy / mag.mass;
        mag.alpha = tau / mag.inertia;
    }
}

function tangential_comp(v, phi) {
    return -v.x * Math.sin(phi) + v.y * Math.cos(phi);
}

function set_coil_state(p_obj, R, phi, vd, ad) {
    p_obj.prop_phi = phi;
    p_obj.prop_p.x = (p_obj.prop_xc || 0) + R * Math.cos(phi);
    p_obj.prop_p.y = (p_obj.prop_yc || 0) + R * Math.sin(phi);
    p_obj.prop_v.x = -vd * Math.sin(phi);
    p_obj.prop_v.y = vd * Math.cos(phi);
    if (ad !== undefined) {
        p_obj.prop_a.x = -ad * Math.sin(phi);
        p_obj.prop_a.y = ad * Math.cos(phi);
        p_obj.prop_anext.x = p_obj.prop_a.x;
        p_obj.prop_anext.y = p_obj.prop_a.y;
    }
}

function step_coil_electron(p_obj, R, vd, ad) {
    set_coil_state(p_obj, R, p_obj.prop_phi + (vd / R) * timestep, vd, ad);
}

function calc_each_vel_and_accelnext() {
    for (var i = 0; i < particles.length; i++) {
        if (!particles[i].prop_iswire && !particles[i].prop_ismagnet) {
            particles[i].prop_atemp.copy(particles[i].prop_a);
        }
    }

    for (var mag of mag_states) {
        mag.ax_old = mag.ax;
        mag.ay_old = mag.ay;
        mag.alpha_old = mag.alpha;
    }

    calc_each_accel(true); // arg foranext

    for (var mag2 of mag_states) {
        mag2.vx += 0.5 * (mag2.ax_old + mag2.ax) * timestep;
        mag2.vy += 0.5 * (mag2.ay_old + mag2.ay) * timestep;
        mag2.omega += 0.5 * (mag2.alpha_old + mag2.alpha) * timestep;
        sync_bar_magnet(mag2);
    }

    for (var i = 0; i < particles.length; i++) {
        var p_i = particles[i];
        if (p_i.prop_iswire || p_i.prop_ismagnet) {
            continue;
        }
        if (p_i.prop_istension) {
            var t_idx = p_i.prop_lead_idx;
            if (i === t_idx) {
                var ax_prev = 0.5 * (particles[t_idx].prop_atemp.x + particles[t_idx + 1].prop_atemp.x);
                var ax_next = 0.5 * (particles[t_idx].prop_anext.x + particles[t_idx + 1].prop_anext.x);
                var vx_new = particles[t_idx].prop_v.x + 0.5 * (ax_prev + ax_next) * timestep;
                particles[t_idx].prop_v.set(vx_new, 0, 0);
                particles[t_idx + 1].prop_v.set(vx_new + particles[t_idx].prop_vd, 0, 0);
                particles[t_idx].prop_anext.set(ax_next, 0, 0);
                particles[t_idx + 1].prop_anext.set(ax_next, 0, 0);
            }
            continue;
        }
        if (p_i.prop_is_par_lead) {
            var idx = p_i.prop_lead_idx;
            if (i === idx) {
                var ay_prev = 0.5 * (particles[idx].prop_atemp.y + particles[idx + 1].prop_atemp.y);
                var ay_next = 0.5 * (particles[idx].prop_anext.y + particles[idx + 1].prop_anext.y);
                var vy_new = particles[idx].prop_v.y + 0.5 * (ay_prev + ay_next) * timestep;
                particles[idx].prop_v.set(0, vy_new, 0);
                particles[idx + 1].prop_v.set(particles[idx].prop_vd, vy_new, 0);
            }
            continue;
        }
        if (p_i.prop_is_sec_lead) {
            var phi_m = p_i.prop_phi;
            var aphi_prev = tangential_comp(p_i.prop_atemp, phi_m);
            var aphi_next = tangential_comp(p_i.prop_anext, phi_m);
            p_i.prop_vphi += 0.5 * (aphi_prev + aphi_next) * timestep;
            p_i.prop_v.set(-p_i.prop_vphi * Math.sin(phi_m), p_i.prop_vphi * Math.cos(phi_m), 0);
            continue;
        }
        p_i.prop_v.add(p_i.prop_atemp.add(p_i.prop_anext).multiplyScalar(0.5 * timestep));
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

        if (p_i.prop_istension) {
            var t_idx = p_i.prop_lead_idx;
            if (i === t_idx) {
                var ax_curr = 0.5 * (particles[t_idx].prop_a.x + particles[t_idx + 1].prop_a.x);
                var x_new = particles[t_idx].prop_p.x + particles[t_idx].prop_v.x * timestep + 0.5 * ax_curr * (timestep ** 2);
                particles[t_idx].prop_p.x = x_new;
                particles[t_idx + 1].prop_p.x = x_new;
            }
            continue;
        }

        if (p_i.prop_ispar) {
            var idx = p_i.prop_lead_idx;
            if (p_i.prop_is_par_lead) {
                if (i === idx) {
                    var ay_curr = 0.5 * (particles[idx].prop_a.y + particles[idx + 1].prop_a.y);
                    var y_new = particles[idx].prop_p.y + particles[idx].prop_v.y * timestep + 0.5 * ay_curr * (timestep ** 2);
                    particles[idx].prop_p.y = y_new;
                    particles[idx + 1].prop_p.y = y_new;
                    // record_trail(particles[idx]);
                }
            } else {
                p_i.prop_p.y = particles[idx].prop_p.y;
                p_i.prop_v.y = particles[idx].prop_v.y;
                if (p_i.prop_charge < 0) {
                    p_i.prop_p.x = wrap_coord(p_i.prop_p.x + p_i.prop_vd * timestep, PAR_L);
                }
            }
            continue;
        }

        if (p_i.prop_isprimary) {
            if (p_i.prop_charge < 0) {
                var phi_p = p_i.prop_phi0 - TFR_AMP_PHI1 * Math.cos(omega_tfr * t_next);
                var vphi_p = p_i.prop_r * TFR_AMP_PHI1 * omega_tfr * Math.sin(omega_tfr * t_next);
                var aphi_p = p_i.prop_r * TFR_AMP_PHI1 * (omega_tfr ** 2) * Math.cos(omega_tfr * t_next);
                set_coil_state(p_i, p_i.prop_r, phi_p, vphi_p, aphi_p);
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
                set_coil_state(p_i, p_i.prop_r, p_i.prop_phi0 + particles[0].prop_phi, particles[0].prop_vphi);
            }
            continue;
        }

        if (p_i.prop_iscoil_turnon) {
            if (p_i.prop_charge < 0) {
                step_coil_electron(p_i, p_i.prop_r, coil_vd_curr, coil_ad_curr);
            }
            continue;
        }

        if (p_i.prop_isacwire) {
            if (p_i.prop_charge < 0) {
                p_i.prop_p.x = wrap_coord(p_i.prop_x0 - AC_AMP * Math.cos(omega_ac * t_next), AC_WIRE_L);
                p_i.prop_v.x = AC_AMP * omega_ac * Math.sin(omega_ac * t_next);
                var ax = AC_AMP * (omega_ac ** 2) * Math.cos(omega_ac * t_next);
                p_i.prop_a.x = ax;
                p_i.prop_anext.x = ax;
            }
            continue;
        }

        if (p_i.prop_iscoil) {
            if (p_i.prop_vz) {
                p_i.prop_p.z += p_i.prop_vz * timestep;
            }
            if (p_i.prop_charge < 0) {
                step_coil_electron(p_i, p_i.prop_r, p_i.prop_vd);
            }
            continue;
        }

        if (p_i.prop_isdcwire) {
            if (p_i.prop_charge < 0) {
                p_i.prop_p.x = wrap_coord(p_i.prop_p.x + (p_i.prop_vd || WIRE_VD) * timestep, p_i.prop_wire_l || WIRE_L);
            }
            continue;
        }

        if (p_i.prop_iswire) {
            continue;
        }

        p_i.prop_vtemp.copy(p_i.prop_v);
        p_i.prop_atemp.copy(p_i.prop_a);
        p_i.prop_p.add(p_i.prop_vtemp.multiplyScalar(timestep)).add(p_i.prop_atemp.multiplyScalar(0.5 * timestep ** 2));

        // record_trail(p_i);
    }
}

function post_line() {
    var thisconfig = CFG.line;
    var testelec = particles[0];
    // Analytical magnetic force: F_y = q * v0 * B, where B = (mu0 * I) / (2 * pi * d0)
    // Since mu0 / (4 * pi) = k_e / c^2 and I = -e * (thisconfig.Q / thisconfig.DX) * thisconfig.VD:
    // Note: Weber v_diff terms scale as v0^2 - (v0 - vd)^2 = 2*v0*vd * (1 - vd / (2*v0)).
    // F_Bvq keeps only 2*v0*vd and omits the -vd^2 drift term; in real wires vd ~ 1e-4 m/s so vd/(2*v0) ~ 0,
    // whereas here thisconfig.VD = 2.25e3 m/s is exaggerated for visible drift, giving a vd/(2*v0) = 0.94% diff.
    var lambda = thisconfig.Q / thisconfig.DX;
    var fy_lab_Bvq = -KEE * lambda * (2 * thisconfig.V0 * thisconfig.VD) / (CSQUARED * thisconfig.D0);

    log_force_cmp('F_Weber:', testelec.prop_fy0, 'F_qvB (approx):', fy_lab_Bvq);
    hud_log('Deflection Δy:', trim((testelec.prop_p.y - thisconfig.D0) * 100), 'cm');
    hud_log('% diff vd/(2·v0) from visible drift speed');
}

function post_coil() {
    var thisconfig = CFG.coil;
    var testelec_coil = particles[0];
    // Analytical magnetic force at center (0, 0, 0) of a circular coil of radius R in the xy-plane:
    // B_z = (mu0 * I) / (2 * R), I = -thisconfig.N * thisconfig.Q * e * thisconfig.VD / (2 * pi * R)
    // With mu0 / (4 * pi) = k_e / c^2, F_y = (-e) * (v0 x B)_y = +e * v0 * B_z:
    var fy_lab_Bvq = -KEE * thisconfig.N * thisconfig.Q * thisconfig.V0 * thisconfig.VD / (CSQUARED * thisconfig.R ** 2);

    log_force_cmp('F_Weber:', testelec_coil.prop_fy_mid, 'F_qvB (approx):', fy_lab_Bvq);
    hud_log('Deflection Δy:', trim(testelec_coil.prop_p.y * 100), 'cm');
}

function post_line_ac() {
    var thisconfig = CFG.line_ac;
    // Analytical Faraday/Weber inductive force on stationary test electron from AC wire acceleration a_e(0):
    // F_x = -(2 * k_e * e^2 * lambda / c^2) * [asinh(L / (2*d0)) - (L/2)/sqrt((L/2)^2 + d0^2)] * a_e(0)
    var omega_ac = (2 * Math.PI * thisconfig.CYCLES) / (stepsmax * timestep);
    var halfL = thisconfig.L / 2;
    var geom_factor = 2 * (Math.asinh(halfL / thisconfig.D0) - halfL / Math.sqrt(halfL ** 2 + thisconfig.D0 ** 2));
    var fx_lab_ind = -(KEE * (thisconfig.Q / thisconfig.DX) / CSQUARED) * geom_factor * (thisconfig.AMP * omega_ac ** 2);

    log_force_cmp('F_Weber:', particles[0].prop_fx0, 'F_Faraday (approx):', fx_lab_ind);
}

function post_coil_turnon() {
    var thisconfig = CFG.coil_turnon;
    var testelec_coil_turnon = particles[0];
    // Analytical Faraday/Weber tangential inductive force when voltage is turned on (a_phi = d(vd)/dt > 0):
    // F_y = -(KEE * Q * R * r0 * a_phi / c^2) * sum_i sin^2(phi_i) / r_i^3
    var ad_turnon = (thisconfig.VD / thisconfig.TAU) * Math.exp(-timestep / thisconfig.TAU);
    var r0 = testelec_coil_turnon.prop_r0_turnon;
    var sum_loop = 0;
    for (var i = 0; i < thisconfig.N; i++) {
        var phi = (2 * Math.PI * i) / thisconfig.N;
        var dist = Math.sqrt(thisconfig.R ** 2 + r0 ** 2 - 2 * thisconfig.R * r0 * Math.cos(phi));
        sum_loop += (Math.sin(phi) ** 2) / (dist ** 3);
    }
    var fy_lab_turnon = -(KEE * thisconfig.Q * thisconfig.R * r0 * ad_turnon / CSQUARED) * sum_loop;

    log_force_cmp('F_Weber (turn-on):', testelec_coil_turnon.prop_fy_turnon, 'F_Faraday (approx):', fy_lab_turnon);
    hud_log('Tangential kick Δy:', trim(testelec_coil_turnon.prop_p.y * 100), 'cm');
}

function post_coils_ac() {
    var thisconfig = CFG.coils_ac;
    // Analytical Faraday/Neumann mutual-inductance tangential force & EMF on 1-loop secondary from 3-loop AC primary:
    // F_phi = -(KEE * Q1 * R1 * R2 * a_phi1(0) / c^2) * sum_{k=0..N1-1, i=0..N-1} sin^2(phi_i) / r_{ki}^3
    // Induced EMF on secondary loop is E_2 = 2 * pi * R2 * (F_phi / E_CHARGE)
    var omega_tfr = (2 * Math.PI * thisconfig.CYCLES) / (stepsmax * timestep);
    var aphi1_0 = thisconfig.R1 * thisconfig.AMP_PHI1 * (omega_tfr ** 2);

    var sum_pri = 0;
    for (var k = 0; k < thisconfig.N1; k++) {
        var zk = (k - (thisconfig.N1 - 1) / 2) * thisconfig.DZ;
        for (var i = 0; i < thisconfig.N_LOOP; i++) {
            var phi = (2 * Math.PI * i) / thisconfig.N_LOOP;
            var dist = Math.sqrt(thisconfig.R1 ** 2 + thisconfig.R2 ** 2 - 2 * thisconfig.R1 * thisconfig.R2 * Math.cos(phi) + zk ** 2);
            sum_pri += (Math.sin(phi) ** 2) / (dist ** 3);
        }
    }
    var fphi_lab = -(KEE * thisconfig.Q1 * thisconfig.R1 * thisconfig.R2 * aphi1_0 / CSQUARED) * sum_pri;
    var fphi_sim = particles[0].prop_fy0;

    // Convert tangential force per electron into loop EMF (Volts): E_loop = 2 * pi * R2 * F_phi * JOULETOEV
    var v2_emf_lab = 2 * Math.PI * thisconfig.R2 * (-fphi_lab) * JOULETOEV;
    var v2_emf_sim = 2 * Math.PI * thisconfig.R2 * (-fphi_sim) * JOULETOEV;
    var v1_n1turn_equiv = thisconfig.N1 * (v2_emf_sim / thisconfig.N2);

    hud_log('F_Weber (Secondary):', trim(fphi_sim), 'N');
    hud_log('F_Faraday (approx):', trim(fphi_lab), 'N');
    hud_log('Secondary EMF Weber:', trim(v2_emf_sim), 'V');
    hud_log('Secondary EMF Faraday:', trim(v2_emf_lab), 'V');
    hud_log('Diff:', pct_diff(v2_emf_lab, v2_emf_sim));
    hud_log('Ratio V2 / V1:', trim(v2_emf_sim / v1_n1turn_equiv), '(N2/N1 = ' + thisconfig.N2 + '/' + thisconfig.N1 + ')');
}

function post_lines_parallel() {
    var thisconfig = CFG.lines_parallel;
    // Analytical Ampère force on a neutral atom of Wire 2 (at y = +D0/2) and Wire 1 (at y = -D0/2):
    // Both wires carry equal current magnitude |I1| = |I2| = (Q / DX) * e * VD.
    // For finite wire length L: geom_amp = (L/2)/sqrt((L/2)^2 + d0^2) * [1 + d0^2 / ((L/2)^2 + d0^2)] -> 1 as L -> inf
    // F_y2 = -2 * KEE * (Q / DX) * VD * par_vd2_signed / (c^2 * D0) * geom_amp, F_y1 = -F_y2
    var halfL = thisconfig.L / 2;
    var r_end = Math.sqrt(halfL ** 2 + thisconfig.D0 ** 2);
    var geom_amp = (halfL / r_end) * (1 + (thisconfig.D0 ** 2) / (r_end ** 2));
    var fy2_lab_ampere = -2 * KEE * (thisconfig.Q / thisconfig.DX) * thisconfig.VD * par_vd2_signed / (CSQUARED * thisconfig.D0) * geom_amp;
    var fy1_lab_ampere = -fy2_lab_ampere;
    var fy1_sim_weber = particles[0].prop_fy_atom0;
    var fy2_sim_weber = particles[2].prop_fy_atom0;
    var dy1 = particles[0].prop_p.y - particles[0].prop_ystart;
    var dy2 = particles[2].prop_p.y - particles[2].prop_ystart;
    var current_amp = (thisconfig.Q / thisconfig.DX) * (1 / JOULETOEV) * thisconfig.VD;

    hud_log('Current |I1| = |I2|:', trim(current_amp), 'A');
    hud_log('Wire1 (Weber/Ampère):', trim(fy1_sim_weber) + ' / ' + trim(fy1_lab_ampere), 'N');
    hud_log('Wire2 (Weber/Ampère):', trim(fy2_sim_weber) + ' / ' + trim(fy2_lab_ampere), 'N');
    hud_log('Diff:', pct_diff(fy2_lab_ampere, fy2_sim_weber));
    hud_log('Moved |Δy1| / |Δy2|:', trim(Math.abs(dy1) * 100) + ' / ' + trim(Math.abs(dy2) * 100), 'cm');
}

function calc_mag_biotsavart(d0, theta1, vd1, theta2, vd2) {
    var thisconfig = CFG.magnets_attract;
    // Analytical Biot-Savart / Ampère magnetic force F = sum (I2 dl2 x B1) and torque tau = sum r2 x (I2 dl2 x B1)
    var build_segs = function (x_center, theta, vd) {
        var out = [];
        var cosT = Math.cos(theta), sinT = Math.sin(theta);
        for_each_mag_seg(function (s, uc, wc, phi) {
            var cosP = Math.cos(phi), sinP = Math.sin(phi);
            var u = uc + thisconfig.R * cosP, w = wc + thisconfig.R * sinP;
            var dx = s * cosT - u * sinT;
            var dy = s * sinT + u * cosT;
            out.push({
                x: x_center + dx, y: dy, z: w, dx: dx, dy: dy,
                qvx: -thisconfig.Q * vd * sinP * sinT,
                qvy: thisconfig.Q * vd * sinP * cosT,
                qvz: -thisconfig.Q * vd * cosP,
            });
        });
        return out;
    };
    var segs1 = build_segs(-d0 / 2, theta1, vd1);
    var segs2 = build_segs(d0 / 2, theta2, vd2);
    var fx2 = 0, fy2 = 0, tau2 = 0;
    for (var s2 of segs2) {
        var bx = 0, by = 0, bz = 0;
        for (var s1 of segs1) {
            var rx = s2.x - s1.x, ry = s2.y - s1.y, rz = s2.z - s1.z;
            var r3 = (rx * rx + ry * ry + rz * rz) ** 1.5;
            bx += (KEE / CSQUARED) * (s1.qvy * rz - s1.qvz * ry) / r3;
            by += (KEE / CSQUARED) * (s1.qvz * rx - s1.qvx * rz) / r3;
            bz += (KEE / CSQUARED) * (s1.qvx * ry - s1.qvy * rx) / r3;
        }
        var dfx = s2.qvy * bz - s2.qvz * by;
        var dfy = s2.qvz * bx - s2.qvx * bz;
        fx2 += dfx;
        fy2 += dfy;
        tau2 += s2.dx * dfy - s2.dy * dfx;
    }
    return { fx2: fx2, fy2: fy2, tau2: tau2 };
}

function post_magnets_attract() {
    var thisconfig = CFG.magnets_attract;
    var m1 = mag_states[0], m2 = mag_states[1];
    var bs = calc_mag_biotsavart(thisconfig.D0, 0, thisconfig.VD, 0, thisconfig.VD);
    var dx1 = m1.x - m1.x0;
    var dx2 = m2.x - m2.x0;

    hud_log('Mag 1 F (Weber / B-S):', trim(m1.fx0) + ' / ' + trim(-bs.fx2), 'N');
    hud_log('Mag 2 F (Weber / B-S):', trim(m2.fx0) + ' / ' + trim(bs.fx2), 'N');
    hud_log('Diff:', pct_diff(bs.fx2, m2.fx0));
    hud_log('Attraction |Δx1| / |Δx2|:', trim(Math.abs(dx1) * 100) + ' / ' + trim(Math.abs(dx2) * 100), 'cm');
}

function post_magnets_repel() {
    var thisconfig = CFG.magnets_repel;
    var m1 = mag_states[0], m2 = mag_states[1];
    var bs = calc_mag_biotsavart(thisconfig.D0, thisconfig.TILT, CFG.magnets_attract.VD, thisconfig.TILT, -CFG.magnets_attract.VD);
    var dtheta1_deg = ((m1.theta - m1.theta0) * 180) / Math.PI;
    var dtheta2_deg = ((m2.theta - m2.theta0) * 180) / Math.PI;
    var dr1 = Math.hypot(m1.x - m1.x0, m1.y);
    var dr2 = Math.hypot(m2.x - m2.x0, m2.y);

    hud_log('F (Weber / B-S):', trim(m2.fx0) + ' / ' + trim(bs.fx2), 'N');
    hud_log('Torque τ_z (Weber / B-S):', trim(m2.tau0) + ' / ' + trim(bs.tau2), 'N·m');
    hud_log('Diff (F / τ_z):', pct_diff(bs.fx2, m2.fx0) + ' / ' + pct_diff(bs.tau2, m2.tau0));
    hud_log('Rotation Δθ1 / Δθ2:', trim(dtheta1_deg) + '° / ' + trim(dtheta2_deg) + '°');
    hud_log('Repulsion |Δr1| / |Δr2|:', trim(dr1 * 100) + ' / ' + trim(dr2 * 100), 'cm');
}

function post_magnet_coil() {
    var thisconfig = CFG.magnet_coil;
    // Analytical Faraday/Lenz tangential force on stationary coil (radius CFG.coils_ac.R2) from 5x5x9 bar magnet moving at thisconfig.VZ
    var sum_mag = 0;
    for_each_mag_seg(function (s, uc, wc, phi) {
        var xp = uc + CFG.magnets_attract.R * Math.cos(phi), yp = wc + CFG.magnets_attract.R * Math.sin(phi), zk = thisconfig.Z0 + s;
        var dx = CFG.coils_ac.R2 - xp, dy = -yp;
        var dist = Math.sqrt(dx * dx + dy * dy + zk * zk);
        sum_mag += zk * yp * (dx * Math.sin(phi) + yp * Math.cos(phi)) / (dist ** 5);
    });
    var fphi_lab = (3 * KEE * thisconfig.Q * thisconfig.VZ * CFG.magnets_attract.VD / CSQUARED) * sum_mag;
    var fphi_sim = particles[0].prop_fy0;
    var emf_sim = 2 * Math.PI * CFG.coils_ac.R2 * (-fphi_sim) * JOULETOEV;
    var emf_lab = 2 * Math.PI * CFG.coils_ac.R2 * (-fphi_lab) * JOULETOEV;
    var dphi_deg = (particles[0].prop_phi * 180) / Math.PI;

    log_force_cmp('F_Weber (Coil):', fphi_sim, 'F_Faraday (approx):', fphi_lab);
    hud_log('EMF (Weber / Faraday):', trim(emf_sim) + ' / ' + trim(emf_lab), 'V');
    hud_log('Coil Δφ (|v_φ|):', trim(dphi_deg) + '° (' + trim(Math.abs(particles[0].prop_vphi)) + ' m/s)');
}

function post_line_magnet() {
    var thisconfig = CFG.line_magnet;
    var m1 = mag_states[0];
    // Analytical Biot-Savart / Ampère torque from finite DC wire of length thisconfig.L on x-axis (z = 0) onto bar magnet at (0, 0, thisconfig.Z0):
    var b_fac = (2 * KEE * thisconfig.N * thisconfig.Q * thisconfig.VD) / (CSQUARED * thisconfig.L);
    var halfL = thisconfig.L / 2;
    var tau_bs = 0;
    for_each_mag_seg(function (s, uc, wc, phi) {
        var yp = uc + CFG.magnets_attract.R * Math.cos(phi), dz = thisconfig.Z0 + wc + CFG.magnets_attract.R * Math.sin(phi);
        var rho2 = yp * yp + dz * dz;
        var xp_pos = halfL - s, xp_neg = halfL + s;
        var r_pos = Math.sqrt(xp_pos * xp_pos + rho2), r_neg = Math.sqrt(xp_neg * xp_neg + rho2);
        var geom = 0.5 * ((xp_pos / r_pos) ** 3 + (xp_neg / r_neg) ** 3);
        var by = b_fac * geom * dz / rho2;
        var bz = -b_fac * geom * yp / rho2;
        var qvy = CFG.magnets_attract.Q * CFG.magnets_attract.VD * Math.sin(phi);
        var qvz = -CFG.magnets_attract.Q * CFG.magnets_attract.VD * Math.cos(phi);
        var dfx = qvy * bz - qvz * by;
        var dfy = 0.5 * b_fac * yp * (qvy * yp + qvz * dz) * (1 / (r_neg ** 3) - 1 / (r_pos ** 3));
        tau_bs += s * dfy - yp * dfx;
    });
    var dtheta_deg = ((m1.theta - m1.theta0) * 180) / Math.PI;

    hud_log('Torque τ_z Weber:', trim(m1.tau0), 'N·m');
    hud_log('Biot-Savart (approx):', trim(tau_bs), 'N·m');
    hud_log('Diff:', pct_diff(tau_bs, m1.tau0));
    hud_log('Magnet Rotation Δθ:', trim(dtheta_deg) + '° (toward 90°)');
}

function post_line_tension() {
    var thisconfig = CFG.line_tension;
    // Analytical Ampère–Weber longitudinal repulsion between collinear wire segments:
    // F_x(i <- j) = (KEE * q^2 * vd^2) / (c^2 * (x_i - x_j)^2) = (mu0 / 4pi) * (I * dx)^2 / (x_i - x_j)^2
    // Midpoint axial tension T_mid is the sum of longitudinal repulsion between the left half (x < 0) and right half (x > 0):
    var f_mid_ana = 0;
    var f_end_ana = 0;
    var x_end0 = -thisconfig.L / 2 + (thisconfig.N - 0.5) * thisconfig.DX;
    for (var i = 0; i < thisconfig.N; i++) {
        var xi = -thisconfig.L / 2 + (i + 0.5) * thisconfig.DX;
        if (i < thisconfig.N - 1) {
            f_end_ana += (KEE * (thisconfig.Q ** 2) * (thisconfig.VD ** 2)) / (CSQUARED * ((x_end0 - xi) ** 2));
        }
        if (i >= thisconfig.N / 2) {
            for (var j = 0; j < thisconfig.N / 2; j++) {
                var xj = -thisconfig.L / 2 + (j + 0.5) * thisconfig.DX;
                f_mid_ana += (KEE * (thisconfig.Q ** 2) * (thisconfig.VD ** 2)) / (CSQUARED * ((xi - xj) ** 2));
            }
        }
    }
    var f_mid_sim = particles[0].prop_f_mid0;
    var f_end_sim = particles[0].prop_f_end0;
    var dx_end = particles[2 * (thisconfig.N - 1)].prop_p.x - x_end0;
    var l_final = thisconfig.L + 2 * dx_end;

    hud_log('Midpoint Tension (W / A):', trim(f_mid_sim) + ' / ' + trim(f_mid_ana), 'N');
    hud_log('End Segment F (W / A):', trim(f_end_sim) + ' / ' + trim(f_end_ana), 'N');
    hud_log('Diff (Lorentz = 0 N):', pct_diff(f_mid_ana, f_mid_sim));
    hud_log('Wire Length (Δx):', trim(thisconfig.L * 100) + ' → ' + trim(l_final * 100) + ' cm (+' + trim(2 * dx_end * 100) + ' cm)');
}

/*
function post_spherical_shell() {
    var testelec_shell = particles[0];
    // Analytical Weber force exerted by a hollow spherical shell of SHELL_N = 500 charges (+SHELL_Q_PER_PT) at radius SHELL_R:
    // F_shell = -(KEE * SHELL_N * SHELL_Q_PER_PT / (3 * c^2 * SHELL_R)) * a_x
    var f_ext_x = 2 * KEE * SHELL_EXT_Q / (SHELL_EXT_X ** 2);
    var ax0 = testelec_shell.prop_ax0;
    var fshell_lab = -((KEE * SHELL_N * SHELL_Q_PER_PT) / (3 * CSQUARED * SHELL_R)) * ax0;
    var fshell_sim = M_ELEC * ax0 - f_ext_x;

    hud_log('F_ext:', trim(f_ext_x), 'N');
    hud_log('F_shell Weber:', trim(fshell_sim), 'N');
    hud_log('F_shell continuum approx:', trim(fshell_lab), 'N');
    hud_log('Diff:', pct_diff(fshell_lab, fshell_sim));
    hud_log('Electron Δx:', trim(testelec_shell.prop_p.x), 'm');
}
*/

var CONFIG_POST = {
    line: post_line,
    coil: post_coil,
    line_ac: post_line_ac,
    coil_turnon: post_coil_turnon,
    coils_ac: post_coils_ac,
    lines_parallel: post_lines_parallel,
    lines_antiparallel: post_lines_parallel,
    magnets_attract: post_magnets_attract,
    magnets_repel: post_magnets_repel,
    magnet_coil: post_magnet_coil,
    line_magnet: post_line_magnet,
    line_tension: post_line_tension,
    // spherical_shell: post_spherical_shell
};

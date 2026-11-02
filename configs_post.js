function pct_diff(approx, sim) {
    return trim((approx - sim) / Math.abs(sim) * 100) + '%';
}

function log_force_cmp(label_sim, f_sim, label_approx, f_approx) {
    hud_log(label_sim, trim(f_sim), 'N');
    hud_log(label_approx, trim(f_approx), 'N');
    hud_log('Diff:', pct_diff(f_approx, f_sim));
}

function post_wire() {
    var testelec = particles[0];
    // Analytical magnetic force: F_y = q * v0 * B, where B = (mu0 * I) / (2 * pi * d0)
    // Since mu0 / (4 * pi) = k_e / c^2 and I = -e * (WIRE_Q / WIRE_DX) * WIRE_VD:
    // Note: Weber relative velocity terms scale as v0^2 - (v0 - vd)^2 = 2*v0*vd * (1 - vd / (2*v0)).
    // F_Bvq keeps only 2*v0*vd and omits the -vd^2 drift term; in real wires vd ~ 1e-4 m/s so vd/(2*v0) ~ 0,
    // whereas here WIRE_VD = 0.0015c is exaggerated for visible drift, giving a vd/(2*v0) = 0.75% diff.
    var lambda = WIRE_Q / WIRE_DX;
    var fy_lab_Bvq = -KEE * lambda * (2 * TEST_V0 * WIRE_VD) / (CSQUARED * WIRE_D0);

    log_force_cmp('F_Weber:', testelec.prop_fy0, 'F_qvB (field approx):', fy_lab_Bvq);
    hud_log('Note: 0.75% diff is vd/(2·v0) from visible drift speed');
    hud_log('Δy deflection:', trim(testelec.prop_p.y - WIRE_D0), 'm');
}

function post_coil() {
    var testelec_coil = particles[0];
    // Analytical magnetic force at center (0, 0, 0) of a circular coil of radius R in the xy-plane:
    // B_z = (mu0 * I) / (2 * R), I = -N_PER_COIL * COIL_Q * e * COIL_VD / (2 * pi * R)
    // With mu0 / (4 * pi) = k_e / c^2, F_y = (-e) * (v0 x B)_y = +e * v0 * B_z:
    var fy_lab_Bvq = -KEE * N_PER_COIL * COIL_Q * COIL_V0 * COIL_VD / (CSQUARED * COIL_R ** 2);

    log_force_cmp('F_Weber:', testelec_coil.prop_fy_mid, 'F_qvB (field approx):', fy_lab_Bvq);
    hud_log('Δy deflection:', trim(testelec_coil.prop_p.y), 'm');
}

function post_ac_wire() {
    // Analytical Faraday/Weber inductive force on stationary test electron from AC wire acceleration a_e(0):
    // F_x = -(2 * k_e * e^2 * lambda / c^2) * [asinh(L / (2*d0)) - (L/2)/sqrt((L/2)^2 + d0^2)] * a_e(0)
    var omega_ac = (2 * Math.PI * AC_CYCLES) / (stepsmax * timestep);
    var halfL = AC_WIRE_L / 2;
    var geom_factor = 2 * (Math.asinh(halfL / AC_WIRE_D0) - halfL / Math.sqrt(halfL ** 2 + AC_WIRE_D0 ** 2));
    var fx_lab_ind = -(KEE * (AC_Q / AC_WIRE_DX) / CSQUARED) * geom_factor * (AC_AMP * omega_ac ** 2);

    log_force_cmp('F_Weber:', particles[0].prop_fx0, 'F_Faraday (field approx):', fx_lab_ind);
}

function post_coil_turnon() {
    var testelec_coil_turnon = particles[0];
    // Analytical Faraday/Weber tangential inductive force when voltage is turned on (a_phi = d(vd)/dt > 0):
    // F_y = -(KEE * COIL_TURNON_Q * COIL_TURNON_R * r0 * a_phi / c^2) * sum_i sin^2(phi_i) / r_i^3
    var ad_turnon = (COIL_TURNON_VD / COIL_TURNON_TAU) * Math.exp(-timestep / COIL_TURNON_TAU);
    var r0 = testelec_coil_turnon.prop_r0_turnon;
    var sum_loop = 0;
    for (var i = 0; i < COIL_TURNON_N; i++) {
        var phi = (2 * Math.PI * i) / COIL_TURNON_N;
        var dist = Math.sqrt(COIL_TURNON_R ** 2 + r0 ** 2 - 2 * COIL_TURNON_R * r0 * Math.cos(phi));
        sum_loop += (Math.sin(phi) ** 2) / (dist ** 3);
    }
    var fy_lab_turnon = -(KEE * COIL_TURNON_Q * COIL_TURNON_R * r0 * ad_turnon / CSQUARED) * sum_loop;

    log_force_cmp('F_Weber turn-on:', testelec_coil_turnon.prop_fy_turnon, 'F_Faraday (field approx):', fy_lab_turnon);
    hud_log('Tangential kick Δy:', trim(testelec_coil_turnon.prop_p.y), 'm');
}

function post_transformer() {
    // Analytical Faraday/Neumann mutual-inductance tangential force & EMF on 1-loop secondary from 3-loop AC primary:
    // F_phi = -(KEE * TFR_Q1 * R1 * R2 * a_phi1(0) / c^2) * sum_{k=0..N1-1, i=0..N-1} sin^2(phi_i) / r_{ki}^3
    // Induced EMF on secondary loop is E_2 = 2 * pi * R2 * (F_phi / E_CHARGE)
    var omega_tfr = (2 * Math.PI * TFR_CYCLES) / (stepsmax * timestep);
    var aphi1_0 = TFR_R1 * TFR_AMP_PHI1 * (omega_tfr ** 2);

    var sum_pri = 0;
    for (var k = 0; k < TFR_N1; k++) {
        var zk = (k - (TFR_N1 - 1) / 2) * TFR_DZ;
        for (var i = 0; i < TFR_N_PER_LOOP; i++) {
            var phi = (2 * Math.PI * i) / TFR_N_PER_LOOP;
            var dist = Math.sqrt(TFR_R1 ** 2 + TFR_R2 ** 2 - 2 * TFR_R1 * TFR_R2 * Math.cos(phi) + zk ** 2);
            sum_pri += (Math.sin(phi) ** 2) / (dist ** 3);
        }
    }
    var fphi_lab = -(KEE * TFR_Q1 * TFR_R1 * TFR_R2 * aphi1_0 / CSQUARED) * sum_pri;
    var fphi_sim = particles[0].prop_fy0;

    // Convert tangential force per electron into loop EMF (Volts): E_loop = 2 * pi * R2 * F_phi * JOULETOEV
    var v2_emf_lab = 2 * Math.PI * TFR_R2 * (-fphi_lab) * JOULETOEV;
    var v2_emf_sim = 2 * Math.PI * TFR_R2 * (-fphi_sim) * JOULETOEV;
    var v1_n1turn_equiv = TFR_N1 * (v2_emf_sim / TFR_N2);

    hud_log('F_Weber (' + TFR_N2 + '-loop Secondary):', trim(fphi_sim), 'N');
    hud_log('F_Faraday (field approx):', trim(fphi_lab), 'N');
    hud_log(TFR_N2 + '-loop Secondary EMF Weber:', trim(v2_emf_sim), 'V');
    hud_log(TFR_N2 + '-loop Secondary EMF Faraday approx:', trim(v2_emf_lab), 'V');
    hud_log('Diff:', pct_diff(v2_emf_lab, v2_emf_sim));
    hud_log('Voltage Ratio V2 (' + TFR_N2 + ' loop) / V1 (' + TFR_N1 + ' loops):', trim(v2_emf_sim / v1_n1turn_equiv), '(N2/N1 = ' + TFR_N2 + '/' + TFR_N1 + ')');
}

function post_parallel_wires() {
    // Analytical Ampère force on a neutral atom of Wire 2 (at y = +PAR_D0/2) and Wire 1 (at y = -PAR_D0/2):
    // Both wires carry equal current magnitude |I1| = |I2| = (PAR_Q / PAR_DX) * e * PAR_VD.
    // For finite wire length L: geom_amp = (L/2)/sqrt((L/2)^2 + d0^2) * [1 + d0^2 / ((L/2)^2 + d0^2)] -> 1 as L -> inf
    // F_y2 = -2 * KEE * (PAR_Q / PAR_DX) * PAR_VD * par_vd2_signed / (c^2 * PAR_D0) * geom_amp, F_y1 = -F_y2
    var halfL = PAR_L / 2;
    var r_end = Math.sqrt(halfL ** 2 + PAR_D0 ** 2);
    var geom_amp = (halfL / r_end) * (1 + (PAR_D0 ** 2) / (r_end ** 2));
    var fy2_lab_ampere = -2 * KEE * (PAR_Q / PAR_DX) * PAR_VD * par_vd2_signed / (CSQUARED * PAR_D0) * geom_amp;
    var fy1_lab_ampere = -fy2_lab_ampere;
    var fy1_sim_weber = particles[0].prop_fy_atom0;
    var fy2_sim_weber = particles[2].prop_fy_atom0;
    var dy1 = particles[0].prop_p.y - particles[0].prop_ystart;
    var dy2 = particles[2].prop_p.y - particles[2].prop_ystart;
    var current_amp = (PAR_Q / PAR_DX) * (1 / JOULETOEV) * PAR_VD;

    hud_log('Current magnitude |I1| = |I2|:', trim(current_amp), 'A');
    hud_log('Wire 1 F_Weber:', trim(fy1_sim_weber), 'N | Ampère approx:', trim(fy1_lab_ampere), 'N');
    hud_log('Wire 2 F_Weber:', trim(fy2_sim_weber), 'N | Ampère approx:', trim(fy2_lab_ampere), 'N');
    hud_log('Diff:', pct_diff(fy2_lab_ampere, fy2_sim_weber));
    hud_log('Distance moved |Δy1|:', trim(Math.abs(dy1)), 'm | |Δy2|:', trim(Math.abs(dy2)), 'm');
}

function calc_mag_biotsavart(d0, theta1, vd1, theta2, vd2) {
    // Analytical Biot-Savart / Ampère magnetic force F = sum (I2 dl2 x B1) and torque tau = sum r2 x (I2 dl2 x B1)
    var build_segs = function (x_cm, theta, vd) {
        var out = [];
        var cosT = Math.cos(theta), sinT = Math.sin(theta);
        for_each_mag_seg(function (k, i, s, phi) {
            var cosP = Math.cos(phi), sinP = Math.sin(phi);
            var dx = s * cosT - MAG_R * cosP * sinT;
            var dy = s * sinT + MAG_R * cosP * cosT;
            out.push({
                x: x_cm + dx, y: dy, z: MAG_R * sinP, dx: dx, dy: dy,
                qvx: -MAG_Q * vd * sinP * sinT,
                qvy: MAG_Q * vd * sinP * cosT,
                qvz: -MAG_Q * vd * cosP,
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
    var m1 = mag_states[0], m2 = mag_states[1];
    var bs = calc_mag_biotsavart(MAG_D0, 0, MAG_VD, 0, MAG_VD);
    var dx1 = m1.x - m1.x0;
    var dx2 = m2.x - m2.x0;

    hud_log('Magnet 1 F_x Weber:', trim(m1.fx0), 'N | Biot-Savart approx:', trim(-bs.fx2), 'N');
    hud_log('Magnet 2 F_x Weber:', trim(m2.fx0), 'N | Biot-Savart approx:', trim(bs.fx2), 'N');
    hud_log('Diff:', pct_diff(bs.fx2, m2.fx0));
    hud_log('Attraction |Δx1|:', trim(Math.abs(dx1)), 'm | |Δx2|:', trim(Math.abs(dx2)), 'm');
}

function post_magnets_repel() {
    var m1 = mag_states[0], m2 = mag_states[1];
    var bs = calc_mag_biotsavart(MAG_D0_REPEL, MAG_TILT, MAG_VD, MAG_TILT, -MAG_VD);
    var dtheta1_deg = ((m1.theta - m1.theta0) * 180) / Math.PI;
    var dtheta2_deg = ((m2.theta - m2.theta0) * 180) / Math.PI;
    var dr1 = Math.hypot(m1.x - m1.x0, m1.y);
    var dr2 = Math.hypot(m2.x - m2.x0, m2.y);

    hud_log('Magnet 2 F_x Weber:', trim(m2.fx0), 'N | Biot-Savart approx:', trim(bs.fx2), 'N');
    hud_log('Torque τ_z Weber:', trim(m2.tau0), 'N·m | Biot-Savart approx:', trim(bs.tau2), 'N·m');
    hud_log('Force Diff:', pct_diff(bs.fx2, m2.fx0), '| Torque Diff:', pct_diff(bs.tau2, m2.tau0));
    hud_log('Rotation Δθ1:', trim(dtheta1_deg), 'deg | Δθ2:', trim(dtheta2_deg), 'deg');
    hud_log('Repulsion |Δr1|:', trim(dr1), 'm | |Δr2|:', trim(dr2), 'm');
}

function post_magnet_coil() {
    // Analytical Faraday/Lenz law tangential force on stationary coil (radius TFR_R2) from bar magnet moving at MCOIL_VZ:
    // F_phi = (3 * KEE * MAG_Q * MAG_R * TFR_R2 * MCOIL_VZ * MAG_VD / c^2) * sum_{k, i} z_k * sin^2(phi_i) / r_{ki}^5
    var sum_mag = 0;
    for_each_mag_seg(function (k, i, s, phi) {
        var zk = MCOIL_Z0 + s;
        var dist = Math.sqrt(MAG_R ** 2 + TFR_R2 ** 2 - 2 * MAG_R * TFR_R2 * Math.cos(phi) + zk ** 2);
        sum_mag += zk * (Math.sin(phi) ** 2) / (dist ** 5);
    });
    var fphi_lab = (3 * KEE * MAG_Q * MAG_R * TFR_R2 * MCOIL_VZ * MAG_VD / CSQUARED) * sum_mag;
    var fphi_sim = particles[0].prop_fy0;
    var emf_sim = 2 * Math.PI * TFR_R2 * (-fphi_sim) * JOULETOEV;
    var emf_lab = 2 * Math.PI * TFR_R2 * (-fphi_lab) * JOULETOEV;
    var dphi_deg = (particles[0].prop_phi * 180) / Math.PI;

    log_force_cmp('F_Weber (Coil):', fphi_sim, 'F_Faraday (field approx):', fphi_lab);
    hud_log('Induced Coil EMF Weber:', trim(emf_sim), 'V | Faraday approx:', trim(emf_lab), 'V');
    hud_log('Coil Electron Rotation Δφ:', trim(dphi_deg), 'deg | Induced |v_φ|/c:', trim(Math.abs(particles[0].prop_vphi) / CSPEED));
}

function post_wire_magnet() {
    var m1 = mag_states[0];
    // Analytical Biot-Savart torque from straight DC wire on x-axis (z = 0) onto bar magnet at (0, 0, OERSTED_Z0):
    var b_fac = (2 * KEE * OERSTED_N * OERSTED_Q * OERSTED_VD) / (CSQUARED * OERSTED_L);
    var tau_bs = 0;
    for_each_mag_seg(function (k, i, s, phi) {
        var yp = MAG_R * Math.cos(phi), dz = OERSTED_Z0 + MAG_R * Math.sin(phi);
        var rho2 = yp * yp + dz * dz;
        var by = b_fac * dz / rho2;
        var bz = -b_fac * yp / rho2;
        var qvy = MAG_Q * MAG_VD * Math.sin(phi);
        var qvz = -MAG_Q * MAG_VD * Math.cos(phi);
        var dfx = qvy * bz - qvz * by;
        tau_bs += -yp * dfx;
    });
    var dtheta_deg = ((m1.theta - m1.theta0) * 180) / Math.PI;

    hud_log('Torque τ_z Weber:', trim(m1.tau0), 'N·m | Biot-Savart approx:', trim(tau_bs), 'N·m');
    hud_log('Diff:', pct_diff(tau_bs, m1.tau0));
    hud_log('Magnet Rotation Δθ:', trim(dtheta_deg), 'deg (twists toward 90° perp to wire)');
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
    wire: post_wire,
    coil: post_coil,
    ac_wire: post_ac_wire,
    coil_turnon: post_coil_turnon,
    transformer: post_transformer,
    parallel_wires: post_parallel_wires,
    antiparallel_wires: post_parallel_wires,
    magnets_attract: post_magnets_attract,
    magnets_repel: post_magnets_repel,
    magnet_coil: post_magnet_coil,
    wire_magnet: post_wire_magnet,
    // spherical_shell: post_spherical_shell
};

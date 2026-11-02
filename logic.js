var totaltime;
var steps;
var stepsmax;

const particles = [];

var timestep; // s. make smaller for larger atoms.
// timestep = 1e-20;
timestep = 2e-10;
// timestep = 5e-20;
// timestep = 2e-21;

var current_config = 'wire';

function reset_sim_state() {
    particles.length = 0;
    mag_states.length = 0;
    steps = 0;
    totaltime = 0;
    // stepsmax = 1000;
    // stepsmax = 1e6;
    // stepsmax = 5e4;
    stepsmax = 2000;
}

function logic_pre() {
    reset_sim_state();
    if (CONFIG_PRE[current_config]) {
        CONFIG_PRE[current_config]();
    }

    plotallparticles();

    // now that all positions are set and won't change before main loop.
    /*
    // For net-charged environments (Config 8 hollow sphere), iterate fixed-point relaxation at t = 0
    // so self-consistent acceleration converges before step 1 begins.
    var init_passes = (current_config === 'spherical_shell') ? 30 : 1;
    for (var pass = 0; pass < init_passes; pass++) {
        calc_forces();
        calc_each_accel();
    }
    */
    calc_forces();
    calc_each_accel();

    if (particles[0]) {
        particles[0].prop_fx0 = M_ELEC * particles[0].prop_a.x;
        particles[0].prop_fy0 = M_ELEC * particles[0].prop_a.y;
        particles[0].prop_fy_mid = particles[0].prop_fy0;
        if (particles[0].prop_is_par_lead) {
            particles[0].prop_fy_atom0 = M_ELEC * (particles[0].prop_a.y + particles[1].prop_a.y);
            particles[2].prop_fy_atom0 = M_ELEC * (particles[2].prop_a.y + particles[3].prop_a.y);
        }
        if (particles[0].prop_istension) {
            var f_mid0 = 0;
            for (var ts = TENSION_N / 2; ts < TENSION_N; ts++) {
                f_mid0 += particles[2 * ts].prop_f.x + particles[2 * ts + 1].prop_f.x;
            }
            particles[0].prop_f_mid0 = f_mid0;
            particles[0].prop_f_end0 = particles[2 * (TENSION_N - 1)].prop_f.x + particles[2 * (TENSION_N - 1) + 1].prop_f.x;
        }
    }

    for (var m = 0; m < mag_states.length; m++) {
        mag_states[m].fx0 = mag_states[m].fx;
        mag_states[m].fy0 = mag_states[m].fy;
        mag_states[m].tau0 = mag_states[m].tau;
    }

    /*
    var testelec_shell = particles.find(particle => particle.prop_id === 'elec_test_shell');
    if (testelec_shell) {
        testelec_shell.prop_ax0 = testelec_shell.prop_a.x;
    }
    */
}

// note: p, v, a are all set before the loop starts.
function logic_step() {
    // timestep = gettimestep();

    // POSITION CALCULATION

    calc_each_pos();

    // recalc forces
    calc_forces();

    totaltime += timestep;

    // VELOCITY CALCULATION

    calc_each_vel_and_accelnext();

    particles.forEach(particle => particle.prop_a.copy(particle.prop_anext)); // setup for next loop

    if (steps === Math.floor(stepsmax / 2) && particles[0] && particles[0].prop_id === 'elec_test_coil') {
        particles[0].prop_fy_mid = -M_ELEC * Math.hypot(particles[0].prop_a.x, particles[0].prop_a.y);
    }

    if (steps === COIL_TURNON_STEP && particles[0] && particles[0].prop_id === 'elec_test_coil_turnon') {
        hud_log('Voltage turned ON at step', steps);
        particles[0].prop_fy_turnon = M_ELEC * particles[0].prop_a.y;
        particles[0].prop_r0_turnon = particles[0].prop_p.x;
    }
}

function logic_post() {
    plotallparticles();

    if (CONFIG_POST[current_config]) {
        CONFIG_POST[current_config]();
    }
}

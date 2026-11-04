var totaltime;
var steps;
var stepsmax;

const particles = [];

var timestep; // s. make smaller for larger atoms.
// timestep = 1e-20;
timestep = 1e-9;
// timestep = 5e-20;
// timestep = 2e-21;

var current_config = 'line';

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

    // Record t = 0 force/torque snapshots used by CONFIG_POST analytical comparisons at the end of the run
    if (particles[0]) {
        particles[0].prop_fx0 = _fx[0];
        particles[0].prop_fy0 = _fy[0];
        particles[0].prop_fy_mid = _fy[0];
        if (particles[0].prop_kind === 'line_plll_lead') {
            particles[0].prop_fy_atom0 = _fy[0] + _fy[1];
            particles[2].prop_fy_atom0 = _fy[2] + _fy[3];
        }
        if (particles[0].prop_kind === 'line_tension') {
            var f_mid0 = 0;
            for (var ts = CFG.line_tension.N / 2; ts < CFG.line_tension.N; ts++) {
                f_mid0 += _fx[2 * ts] + _fx[2 * ts + 1];
            }
            particles[0].prop_f_mid0 = f_mid0;
            particles[0].prop_f_end0 = _fx[2 * (CFG.line_tension.N - 1)] + _fx[2 * (CFG.line_tension.N - 1) + 1];
        }
    }

    for (var mag of mag_states) {
        mag.fx0 = mag.fx;
        mag.tau0 = mag.tau;
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

    // VELOCITY & ACCELERATION CALCULATION

    calc_each_vel_and_accel();

    // Record mid-run force snapshots (coil center crossing in Config 2, turn-on step in Config 4) for CONFIG_POST
    if (steps === Math.floor(stepsmax / 2) && particles[0] && particles[0].prop_id === 'elec_test_coil') {
        particles[0].prop_fy_mid = -Math.hypot(_fx[0], _fy[0]);
    }

    if (steps === CFG.coil_turnon.STEP && particles[0] && particles[0].prop_id === 'elec_test_coil_turnon') {
        hud_log('Voltage turned ON:', trim(steps * timestep * 1e9), 'ns');
        particles[0].prop_fy_turnon = _fy[0];
        particles[0].prop_r0_turnon = particles[0].prop_p.x;
    }
}

function logic_post() {
    plotallparticles();

    if (CONFIG_POST[current_config]) {
        CONFIG_POST[current_config]();
    }
}

const stepsuntilrender = 10;
var is_running = false;
var config_list_el;
var desc_el;
var status_step_el;
var status_bar_el;
var log_el;
var btn_start;
var btn_start_label_el;
var btn_pause;
var config_btns = [];

function set_status(state_text) {
    var t_cur_ns = trim(steps * timestep * 1e9);
    var is_active = state_text === 'Running' || state_text === 'Paused';
    var is_done = state_text === 'Completed';
    btn_start.classList.toggle('is-active', is_active);
    btn_start.classList.toggle('is-done', is_done);
    status_step_el.textContent = (is_active || is_done) ? t_cur_ns + ' ns' : '';
    status_bar_el.style.width = (stepsmax > 0 ? (steps / stepsmax) * 100 : 0) + '%';
}

function update_config_scroll_fade() {
    var at_bottom = config_list_el.scrollTop + config_list_el.clientHeight >= config_list_el.scrollHeight - 2;
    config_list_el.classList.toggle('at-bottom', at_bottom);
}

function reset_and_prepare(config_name) {
    is_running = false;
    if (config_name) {
        current_config = config_name;
    }

    // Clear Three.js reference line segments from previous config (keep persistent point clouds)
    for (var i = scene.children.length - 1; i >= 0; i--) {
        var obj = scene.children[i];
        if (obj === points_nuke || obj === points_elec) continue;
        if (obj.geometry) obj.geometry.dispose();
        if (obj.material) obj.material.dispose();
        scene.remove(obj);
    }

    logic_pre();

    hud_clear();
    var lines = (config_descriptions[current_config] || '').split('\n');
    desc_el.innerHTML = lines.map(function (line) {
        return '<div class="desc-row"><span class="desc-bullet">•</span><span>' + line.replace(/^•\s*/, '') + '</span></div>';
    }).join('');
    log_el.innerHTML = '<span class="hud-placeholder">Click Start or press Space to run</span>';
    set_status('Ready');
    btn_start_label_el.textContent = 'Start';
    btn_pause.textContent = 'Pause';
    btn_pause.disabled = true;

    renderer.render(scene, camera);
}

function start_simulation() {
    // If a run already finished or progressed, reset cleanly before running
    if (steps > 0) {
        reset_and_prepare(current_config);
    }
    hud_clear();
    log_el.innerHTML = '<span class="hud-placeholder">Running simulation...</span>';
    is_running = true;
    btn_start_label_el.textContent = 'Restart';
    btn_pause.textContent = 'Pause';
    btn_pause.disabled = false;
    set_status('Running');
}

function select_config_btn(btn) {
    config_btns.forEach(function (b) { b.classList.remove('active'); });
    btn.classList.add('active');
    btn.scrollIntoView({ block: 'nearest' });
    update_config_scroll_fade();
    reset_and_prepare(btn.getAttribute('data-config'));
}

function animate_loop() {
    if (is_running) {
        for (var i = 0; i < stepsuntilrender && steps < stepsmax; i++) {
            logic_step();
            steps++;
        }
        plotallparticles();

        // >= not just ===, because stepsuntilrender*ninteger might be > stepsmax
        if (steps >= stepsmax) {
            is_running = false;
            btn_start_label_el.textContent = 'Start Again';
            set_status('Completed');
            btn_pause.disabled = true;
            logic_post();
        } else {
            set_status('Running');
        }
    }

    controls.update();
    renderer.render(scene, camera);
}

function init_viz() {
    config_list_el = document.getElementById('config-list');
    desc_el = document.getElementById('hud-desc');
    status_step_el = document.getElementById('status-step');
    status_bar_el = document.getElementById('status-bar-fill');
    log_el = document.getElementById('hud-log');
    btn_start = document.getElementById('btn-start');
    btn_start_label_el = document.getElementById('btn-start-label');
    btn_pause = document.getElementById('btn-pause');
    config_btns = Array.from(document.querySelectorAll('.config-btn'));

    config_list_el.addEventListener('scroll', update_config_scroll_fade);

    config_btns.forEach(function (btn) {
        btn.addEventListener('click', function () {
            btn.blur();
            select_config_btn(btn);
        });
    });

    btn_start.addEventListener('click', function () {
        btn_start.blur();
        start_simulation();
    });

    btn_pause.addEventListener('click', function () {
        btn_pause.blur();
        if (steps === 0 || steps >= stepsmax) return;
        is_running = !is_running;
        btn_pause.textContent = is_running ? 'Pause' : 'Resume';
        set_status(is_running ? 'Running' : 'Paused');
    });

    window.addEventListener('keydown', function (e) {
        if (e.metaKey || e.ctrlKey || e.altKey) return;
        if (e.code === 'Space') {
            e.preventDefault();
            if (document.activeElement && document.activeElement.blur) {
                document.activeElement.blur();
            }
            if (steps === 0 || steps >= stepsmax) {
                start_simulation();
            } else {
                btn_pause.click();
            }
        } else if (e.code === 'ArrowDown' || e.code === 'ArrowUp') {
            e.preventDefault();
            var idx = config_btns.findIndex(function (b) { return b.classList.contains('active'); });
            var next_idx = (idx + (e.code === 'ArrowDown' ? 1 : -1) + config_btns.length) % config_btns.length;
            select_config_btn(config_btns[next_idx]);
        }
    });

    var url_config = new URLSearchParams(window.location.search).get('config');
    if (url_config && CONFIG_PRE[url_config]) {
        current_config = url_config;
        config_btns.forEach(function (b) {
            var is_match = b.getAttribute('data-config') === current_config;
            b.classList.toggle('active', is_match);
            if (is_match) b.scrollIntoView({ block: 'center' });
        });
        update_config_scroll_fade();
    }

    reset_and_prepare(current_config);
    renderer.setAnimationLoop(animate_loop);
}

var scene;
var camera;
var renderer;
var controls;

async function setup() {
    scene = new THREE.Scene();
    // var scene_height = window.innerHeight * 0.8;
    var scene_height = window.innerHeight * 1;

    camera = new THREE.PerspectiveCamera(
        35,
        window.innerWidth / scene_height,
        // 0.1 * 10 ** -11,
        // 1000 * 10 ** -11
        0.001,
        10
    );

    // camera.position.set(0, 0, 6 * BOHR_X);
    camera.position.set(0, 0, 0.06);
    camera.lookAt(new THREE.Vector3(0, 0, 0));

    renderer = new THREE.WebGPURenderer({
        antialias: true,
        alpha: false,
        stencil: false,
        powerPreference: 'high-performance',
    });
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(window.innerWidth, scene_height);

    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;

    // call controls.update() after a manual change to the camera
    controls.update();

    document.querySelector('#thescene').appendChild(renderer.domElement);

    window.addEventListener('resize', () => {
        // scene_height = window.innerHeight * 0.8;
        scene_height = window.innerHeight * 1;

        camera.aspect = window.innerWidth / scene_height;
        camera.updateProjectionMatrix();
        renderer.setSize(window.innerWidth, scene_height);
    });

    await renderer.init();
    // Check this only AFTER renderer.init() resolves
    const isUsingWebGPU = renderer.backend.isWebGPUBackend;
    console.log(`Active Pipeline: ${isUsingWebGPU ? 'WebGPU' : 'WebGL 2 fallback'}`);

    init_viz();
}

if (typeof window !== 'undefined') {
    setup();
}

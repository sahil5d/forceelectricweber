# Force_electric_Weber

Simulate electromagnetic and magnetic phenomena using **only** Wilhelm Weber's direct charge-to-charge force law (1846):

$$F = \frac{k q_1 q_2}{r^2}\left(1 - \frac{\dot{r}^2}{2c^2} + \frac{r\ddot{r}}{c^2}\right) = \frac{k q_1 q_2}{r^2}\left(1 - \frac{v_r^2}{2c^2} + \frac{v_\perp^2}{c^2} + \frac{r a_r}{c^2}\right)$$

Every wire (line), coil, and bar magnet in the scene is modeled explicitly as neutral pairs of positive ions (protons) and mobile/circulating electrons. (Note: wires (lines), coils, and magnets are modeled as macroscopic charge segments rather than at the individual atomic level; electrons are constrained to drift along the protons of their conductor, and internal forces within the same wire (line), coil, or magnet are not computed.) At each timestep, the simulation sums only the pairwise central Weber force along the line connecting each pair of charges on interacting bodies—reproducing magnetic deflection, Faraday induction, transformer action, Ampère wire forces, permanent-magnet forces/torques, and Ørsted compass alignment **without** intermediate $\vec{B}$ fields, $\vec{A}$ potentials, or separate Lorentz/Faraday rules.

---

### Pairwise Weber Kinematics (`calc_forces.js`)

For each pair of charges $q_1, q_2$ with relative position $\vec{r} = \vec{r}_2 - \vec{r}_1$, relative velocity $\vec{v} = \vec{v}_2 - \vec{v}_1$, and relative acceleration $\vec{a} = \vec{a}_2 - \vec{a}_1$:

- **Distance & direction:** $r = |\vec{r}|$, $\quad \hat{r} = \vec{r} / r$
- **Range rate ($\dot{r}$):** $\dot{r} = \hat{r} \cdot \vec{v}$ (scalar radial velocity)
- **Second derivative of distance ($\ddot{r}$):**
  $$\ddot{r} = \frac{v^2 - \dot{r}^2}{r} + \hat{r} \cdot \vec{a}$$
  Note that $\ddot{r} = \frac{d^2 r}{dt^2}$ is a scalar and differs from the radial projection of acceleration $\hat{r}\cdot\vec{a}$ by the transverse centripetal/flyby term $v_\perp^2 / r = (v^2 - \dot{r}^2)/r$.

---

### Phenomena Simulated

| # | Configuration | Scene Setup | Emergent Weber Behavior | Textbook Field Law Replaced |
| :--- | :--- | :--- | :--- | :--- |
| **1** | **Line** | Lone electron ($v = 0.1c$) moving parallel to a $600\text{ m}$ neutral DC wire | Deflects perpendicular to the wire via relative-velocity terms ($\dot{r}^2/2c^2$, $v_\perp^2/c^2$) | Lorentz force $\vec{F} = q\vec{v} \times \vec{B}$ around a straight wire |
| **2** | **Coil** | Lone electron ($v = 0.05c$) passing through the center of a coil ($R = 8\text{ m}$) | Deflects sideways in the coil plane from pairwise forces with circulating coil electrons | Biot-Savart coil field $\vec{B}$ + $\vec{F} = q\vec{v} \times \vec{B}$ |
| **3** | **Line AC** | Stationary electron near a $120\text{ m}$ AC wire with oscillating electrons | Accelerates parallel/antiparallel to the wire via the $r(\hat{r}\cdot\vec{a})/c^2$ term | Faraday inductive electric field $\vec{E} = -\partial\vec{A}/\partial t$ |
| **4** | **Coil Startup** | Stationary electron inside a circular coil ($R = 8\text{ m}$) whose current turns on at step 600 | Receives a tangential inductive kick opposite the ramping electron acceleration | Faraday/Lenz changing-flux induction $\nabla \times \vec{E} = -\partial\vec{B}/\partial t$ |
| **5** | **Coils AC (Transformer)** | 3-loop AC primary coil ($R_1 = 8\text{ m}$) coaxial with a 1-loop secondary coil ($R_2 = 6.5\text{ m}$) | Primary electron acceleration directly drives secondary coil current with EMF ratio $V_2/V_1 \approx 1/3$ | Mutual inductance & magnetic flux linkage ($V_2/V_1 = N_2/N_1$) |
| **6** | **Lines Parallel** | Two $120\text{ m}$ neutral DC wires at $y = \pm 2\text{ m}$ carrying currents in the same direction | Both wires attract and move symmetrically toward $y = 0$ | Ampère's parallel-wire attraction $\vec{F} = I\vec{L} \times \vec{B}$ |
| **7** | **Lines Parallel (Opp Current)** | Two $120\text{ m}$ neutral DC wires at $y = \pm 2\text{ m}$ carrying opposite currents | Both wires repel and move symmetrically apart | Ampère's antiparallel-wire repulsion |
| **8** | **Bar Magnets Attract** | Two 10-loop Ampèrian bar magnets ($R = 1.4\text{ m}$) with aligned loop currents (N–S facing) | Attract each other along their shared axis | Magnetic pole / dipole-dipole attraction |
| **9** | **Bar Magnets Repel** | Two 10-loop bar magnets with opposing currents (N–N facing), tilted $6^\circ$ | Simultaneously repel apart and rotate via pairwise torque | Magnetic dipole repulsion & torque $\vec{\tau} = \vec{\mu} \times \vec{B}$ |
| **10** | **Bar Magnet and Coil (Induction)** | 10-loop bar magnet ($v_z = 0.06c$) moving along $+z$ through a conducting coil ($R = 6.5\text{ m}$) | Relative motion ($\dot{r}^2/2c^2$, $r\ddot{r}/c^2$) induces loop current and EMF in the stationary coil | Faraday moving-magnet induction $\mathcal{E} = -d\Phi_B/dt$ |
| **11** | **Line and Bar Magnet (Compass)** | 10-loop bar magnet at $(0,0,-2.8\text{ m})$ initially parallel to a $45\text{ m}$ DC wire on the $x$-axis | Pure pairwise Weber torque twists the bar magnet toward $90^\circ$ (perpendicular to the wire) | Ørsted's experiment / compass alignment in a wire's $\vec{B}$ field |

---

### Running & Code Structure

Open `index.html` directly in any browser, select a configuration in the top-right panel, and click **Start**. At completion (2000 steps), the bottom-right HUD compares the simulated Weber force, torque, or induced EMF against the standard analytical Lorentz/Biot-Savart/Faraday/Ampère formulas.

- **`calc_forces.js`** — Core $O(N^2)$ pairwise Weber force calculation (`calc_forces()`).
- **`calc.js`** — Velocity-Verlet integration (`calc_each_pos()`, `calc_each_accel()`, `calc_each_vel_and_accelnext()`) for free particles, constrained wire line or coil electrons, and rigid-body bar magnets.
- **`configs_pre.js`** — Particle geometry, initial conditions (`CONFIG_PRE`), and HUD descriptions (`config_descriptions`) for all 11 configurations.
- **`configs_post.js`** — Post-simulation analytical field-theory comparisons (`CONFIG_POST`) logged to the HUD.
- **`logic.js`** — Simulation lifecycle (`logic_pre()`, `logic_step()`, `logic_post()`).
- **`lib.js`** — Physical constants, particle factory (`add_particle()`), and Three.js rendering helpers.

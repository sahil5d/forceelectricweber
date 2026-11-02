// Accumulate pairwise Weber forces directly onto each particle's prop_f
var dir = new vec();
var pdiff = new vec();
var vdiff = new vec();
var adiff = new vec();
var force = new vec(); // the initial -1 is so -1*oppchargeproduct will be 1, -1*samechargeproduct will be -1.
var minr; // save smallest pairwise distance

// call this before loop, and every loop
function calc_forces() {
    minr = 1000; // 1000 meters

    for (var i = 0; i < particles.length; i++) {
        if (!particles[i].prop_iswire) {
            particles[i].prop_f.set(0, 0, 0);
        }
    }

    for (var i = 0; i < particles.length - 1; i++) {
        var pi = particles[i];
        if (pi.prop_iswire) {
            continue;
        }

        // note: could change iloop max from '< .length-1' to '< .length' if need to do a check on the last element in array (like max r from origin). causes no issue on pairdata (inner jloop won't run when i=.length-1)

        for (var j = i + 1; j < particles.length; j++) {
            var pj = particles[j];

            // Ignore internal forces within the same wire, coil, or bar magnet
            if ((pi.prop_group && pi.prop_group === pj.prop_group) ||
                (pi.prop_is_par_lead && pj.prop_is_par_lead)) {
                continue;
            }

            pdiff.subVectors(pj.prop_p, pi.prop_p);
            var distsq = pdiff.lengthSq(); // better than calling .distanceToSquared() which repeats .sub operation
            var dist = Math.sqrt(distsq);
            dir.copy(pdiff).multiplyScalar(1 / dist);

            // kqq/r^2 Coulomb
            force.copy(dir).multiplyScalar(-1 * KEE * pi.prop_charge * pj.prop_charge / distsq);

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

            vdiff.subVectors(pj.prop_v, pi.prop_v);
            adiff.subVectors(pj.prop_a, pi.prop_a);
            var rdot = dir.dot(vdiff);
            var rddot = (vdiff.lengthSq() - rdot ** 2) / dist + dir.dot(adiff);
            var weberfactor = 1 - (rdot ** 2) / (2 * CSQUARED) + (dist * rddot) / CSQUARED;
            force.multiplyScalar(weberfactor);

            pi.prop_f.add(force);
            if (!pj.prop_iswire) {
                pj.prop_f.sub(force);
            }

            // optional
            if (dist < minr) {
                minr = dist;
            }
        }
    }
}

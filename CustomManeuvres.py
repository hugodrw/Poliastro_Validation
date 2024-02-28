"""Custom Orbital maneuvers.

"""
from astropy import units as u

from poliastro.core.maneuver import (
    hohmann as hohmann_fast,
)

from CustomLowLevel import hohmann_any_angle


from poliastro.maneuver import Maneuver
from poliastro.twobody.orbit import Orbit
from InPlanePhysics import delta_u
from poliastro.util import norm, wrap_angle
from poliastro.core.elements import coe_rotation_matrix, rv2coe, rv_pqw


def hohmann_with_phasing(orbit_i: Orbit, orbit_f: Orbit):
    r"""Compute a Hohmann transfer with correct phasing to a target debris.
    For circular orbits only.

    Parameters
    ----------

    """
    # Calculate transfer time for delta_u
    r_f = orbit_f.a
    r_f = r_f.to_value(u.m)
    rv = orbit_i.rv()
    rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))
    k = orbit_i.attractor.k
    _, _, t_trans = hohmann_any_angle(k, rv, r_f)

    # Calculate delta at which the burn should be applied
    target_delta = delta_u(t_trans, orbit_f)
    print('Target Delta: ' , target_delta)

    # Calulate the current delta
    current_delta =  orbit_f.nu -  orbit_i.nu << u.deg
    if current_delta < 0:
        current_delta = 360 * u.deg + current_delta # wrap to 360
    print('current_delta : ' , current_delta)

    # Calculate the angular velocities
    w_i = orbit_i.n.to(u.deg / u.s)
    w_f = orbit_f.n.to(u.deg / u.s)

    # Calculate the time to the first burn
    dist = current_delta - target_delta
    print('dist before: ' , dist)
    if dist < 0:
        dist = 360 * u.deg + dist # wrap to 360
    print('dist after: ' , dist)
    t_1 = dist / (w_i - w_f)
    print('t_1: ' , t_1)

    # Propagate to the first burn
    orbit_i = orbit_i.propagate(t_1)
    orbit_f = orbit_f.propagate(t_1)
    print('orbit_i.nu: ' , orbit_i.nu << u.deg)
    print('new delta: ' , orbit_f.nu - orbit_i.nu << u.deg)


    # Compute delta_v vectors from first burn location
    r_f = orbit_f.a
    r_f = r_f.to_value(u.m)
    rv = orbit_i.rv()
    rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))

    # Calculate hohmann DV and Transfer Time from the first burn location
    k = orbit_i.attractor.k
    dv_a, dv_b, t_trans = hohmann_any_angle(k, rv, r_f)
    dv_a, dv_b, t_trans = dv_a * u.m / u.s, dv_b * u.m / u.s, t_trans * u.s
    t_2 = t_trans

    return Maneuver(
        (t_1.decompose(), dv_a.decompose()),
        (t_2.decompose(), dv_b.decompose()),
    )
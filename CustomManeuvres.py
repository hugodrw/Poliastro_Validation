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

    # Calculate delta at which the burn should be applied
    target_delta = delta_u(orbit_i.a, orbit_f.a)
    print('Target Delta: ' , target_delta)

    # Calulate the current delta
    current_delta =  orbit_f.nu -  orbit_i.nu << u.deg
    if current_delta < 0:
        current_delta = 360 * u.deg + current_delta # wrap to 
    print('current_delta : ' , current_delta)

    # Calculate the angular velocities
    w_i = orbit_i.n.to(u.deg / u.s)
    w_f = orbit_f.n.to(u.deg / u.s)

    # Calculate the time to the first burn
    t_1 = (current_delta - target_delta) / (w_i - w_f)
    print('t_1: ' , t_1)

    # Propagate to the first burn
    orbit_i = orbit_i.propagate(t_1)
    orbit_f = orbit_f.propagate(t_1)
    print('orbit_i.nu: ' , orbit_i.nu << u.deg)
    print('new delta: ' , orbit_f.nu - orbit_i.nu << u.deg)

    # Format for hohmann_fast
    r_f = orbit_f.a
    r_f = r_f.to_value(u.m)
    rv = orbit_i.rv()
    rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))

    # Calculate hohmann DV and Transfer Time from the first burn location
    k = orbit_i.attractor.k
    
    # Debug
    print('rv: ' , rv)

    _, ecc, inc, raan, argp, nu = rv2coe(k, *rv)
    rot_matrix = coe_rotation_matrix(inc, raan, argp)
    print('rot_matrix: ' , rot_matrix)


    dv_a, dv_b, t_trans = hohmann_any_angle(k, rv, r_f)
    dv_a, dv_b, t_trans = dv_a * u.m / u.s, dv_b * u.m / u.s, t_trans * u.s

    t_2 = t_trans


    return Maneuver(
        (t_1.decompose(), dv_a.decompose()),
        (t_2.decompose(), dv_b.decompose()),
    )

    

    # # Calculate Delta-U
    # d_u = delta_u(orbit_i.a, orbit_f.a)

    # print('Delta_U: ' , d_u)
    # print('orbit_i.nu: ' , orbit_i.nu)

    # # Find location of first burn
    # first_burn = orbit_f.nu - d_u

    # # 


    # # t_1 = orbit_i.time_to_anomaly(orbit_f.nu + d_u << u.)
    # # orbit_i = orbit_i.propagate_to_anomaly(((orbit_f.a + d_u) ) * u.deg)


    # k = orbit_i.attractor.k

    # rv = orbit_i.rv()
    # rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))

    # k = k.to_value(u.m**3 / u.s**2)
    # r_f = r_f.to_value(u.m)

    # dv_a, dv_b, t_trans = hohmann_fast(k, rv, r_f)
    # dv_a, dv_b, t_trans = dv_a * u.m / u.s, dv_b * u.m / u.s, t_trans * u.s

    # # t_1 is the time until the first burn

    # return Maneuver(
    #     (t_1.decompose(), dv_a.decompose()),
    #     (t_trans.decompose(), dv_b.decompose()),
    # )



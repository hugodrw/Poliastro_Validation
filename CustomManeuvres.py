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
import numpy as np


def hohmann_with_phasing(orbit_i: Orbit, orbit_f: Orbit, debug=True):
    r"""Compute a Hohmann transfer with correct phasing to a target debris.
    For circular orbits only.

    Parameters
    ----------

    """
    # Downwards Hohmann
    down = False

    # Calculate transfer time for delta_u
    r_f = orbit_f.a
    r_f = r_f.to_value(u.m)
    rv = orbit_i.rv()
    rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))
    k = orbit_i.attractor.k
    _, _, t_trans = hohmann_any_angle(k, rv, r_f)

    # Calculate delta at which the burn should be applied
    target_delta = delta_u(t_trans, orbit_f)
    if target_delta < 0:
        down = True
        target_delta = 360 * u.deg + target_delta # wrap to 360
    print('Target Delta: ' , target_delta)

    # Calulate the current delta
    mean_anomaly_i = (orbit_i.nu) << u.deg
    mean_anomaly_f = (orbit_f.nu) << u.deg
    current_delta =  mean_anomaly_f - mean_anomaly_i << u.deg
    if current_delta < 0:
        current_delta = 360 * u.deg + current_delta # wrap to 360
    print('current_delta : ' , current_delta) if debug else None

    # Calculate the angular velocities
    w_i = orbit_i.n.to(u.deg / u.s)
    w_f = orbit_f.n.to(u.deg / u.s)

    # Calculate the time to the first burn
    dist = current_delta - target_delta if not down else target_delta - current_delta
    print('dist: ', dist)
    if dist < 0:
        dist = 360 * u.deg + dist # wrap to 360
    t_1 = dist / np.abs((w_i - w_f))
    print('t_1: ' , t_1) if debug else None

    # Propagate to the first burn
    orbit_i = orbit_i.propagate(t_1)
    orbit_f = orbit_f.propagate(t_1)
    
    if debug:
        mean_anomaly_i = (orbit_i.nu) << u.deg
        mean_anomaly_f = (orbit_f.nu) << u.deg

        print('mean_anomaly_i: ' , mean_anomaly_i)
        print('new delta: ' , mean_anomaly_f - mean_anomaly_i << u.deg)

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

def simple_inc_change(orbit_i: Orbit, orbit_f: Orbit, debug=True):
    r"""Compute thrust vectors and phase time needed for an inclination change.

    Parameters
    ----------

    """
    # 
    mean_anomaly_i = orbit_i.nu << u.deg

    # Propagate to thrust location (0, 180)
    thrust_location = 0 * u.deg if mean_anomaly_i <= 0 else 179.999 * u.deg
    time_to_thrust = orbit_i.time_to_anomaly(thrust_location)
    if debug:
        print('currrent_anomaly', mean_anomaly_i)
        print('thrust_location', thrust_location)
        print('time_to_thrust', time_to_thrust)
    orbit_i = orbit_i.propagate_to_anomaly(thrust_location)

    # Calculate the thrust value
    v = norm(orbit_i.v << u.m / u.s)
    inc_i = orbit_i.inc << u.rad
    inc_f = orbit_f.inc << u.rad
    inc_delta = inc_f - inc_i
    thrust_norm = 2*v*np.sin((inc_delta << u.rad)/2)

    # Calculate the thrust vector
    y_thrust = np.sin(inc_delta/2)*thrust_norm
    z_thrust = -np.cos(inc_delta/2)*thrust_norm
    # Inverse for 0 degrees
    if thrust_location == 0 * u.deg:
        print('zero deg transformation')
        y_thrust = -y_thrust
        z_thrust = -z_thrust
    
    if debug:
        print('thrust_norm', thrust_norm)
        print('y_thrust', y_thrust)
        print('z_thrust', z_thrust)

    thrust_vector = np.array([0 ,y_thrust.value,z_thrust.value]) * u.m / u.s

    # Use rotation matrix to go from orbital to general referential
    k = orbit_i.attractor.k
    rv = orbit_i.rv()
    rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))
    _, ecc, inc, raan, argp, nu = rv2coe(k, *rv)
    rot_matrix = coe_rotation_matrix(inc, raan, argp)
    thrust_vector = rot_matrix @ thrust_vector

    return Maneuver(
        (time_to_thrust.decompose(), thrust_vector.decompose())
    )
from astropy import units as u

def delta_u(r1, r2):
    # Assuming r1 and r2 are in units of length (e.g., meters)
    r1 = r1.to(u.m)
    r2 = r2.to(u.m)
    
    # Calculate the angle
    angle = u.rad * (1 - (1 - (r1 + r2) / (2 * r2))**0.5)
    
    return angle << u.deg
import math
from datetime import datetime

from scipy.integrate import ode

from igrf12py import igrf12syn as igrf12

from ..gpstk import PyPosition


def differential(t, s):
    r_km, theta, phi = s
    Bt, Bp, Br, B = igrf12(0,
                           date_dt.year,
                           2,  # geocentric
                           r_km,
                           90 - math.degrees(theta),
                           math.degrees(phi) % 360)
    Br *= -1
    dr = Br / B
    dp = Bp / B
    dt = Bt / B
    return [dr,
            dt / r_km,
            dp / (r_km * math.cos(theta))]


def line(date_dt, pos0, dt=1e2):
    """
    ???

    https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field--igrf--model/content/igrfline.m
    """
    r_km0 = pos0.radius / 1e3
    theta0 = math.radians(pos0.geocentricLatitude)
    phi0 = math.radians(pos0.longitude)
    y0, t0 = [r_km0, theta0, phi0], 0
    tracer = ode(differential).set_integrator('dopri5')
    tracer.set_initial_value(y0, t0)
    pos_path = [PyPosition(math.degrees(theta0),
                           math.degrees(phi0),
                           r_km0 * 1e3,
                           PyPosition.CoordinateSystem['geocentric'])]
    while True:
        assert tracer.successful()
        r_km_i, theta_i, phi_i = tracer.integrate(tracer.t + dt)
        pos = PyPosition(math.degrees(theta_i),
                         math.degrees(phi_i),
                         r_km_i * 1e3,
                         PyPosition.CoordinateSystem['geocentric'])
        #print(pos.height)
        if pos.height < 0:
            break
        pos_path.append(pos)
    return pos_path


if __name__ == '__main__':
    date_dt = datetime(2007, 7, 17, 6, 30)
    lat = -60
    lon = 180
    alt = 0
    pos0 = PyPosition(lat, lon, alt,
                      PyPosition.CoordinateSystem['geodetic'])
    pos_path = line(date_dt, pos0)
    print(len(pos_path))

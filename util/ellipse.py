import math

import numpy as NP
from numpy.linalg import eig, inv


"""
From http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html

Reference:

Fitzgibbon, Pilu and Fischer in Fitzgibbon, A.W., Pilu, M., and
Fischer R.B., Direct least squares fitting of ellipsees, Proc. of the
13th Internation Conference on Pattern Recognition, pp 253-257,
Vienna, 1996.
"""


def fit_ellipse(x,y):
    x = x[:,NP.newaxis]
    y = y[:,NP.newaxis]
    D =  NP.hstack((x*x, x*y, y*y, x, y, NP.ones_like(x)))
    S = NP.dot(D.T,D)
    C = NP.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(NP.dot(inv(S), C))
    n = NP.argmax(NP.abs(E))
    a = V[:,n]
    return a


def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return NP.array([x0,y0])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*NP.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*NP.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*NP.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=NP.sqrt(up/down1)
    res2=NP.sqrt(up/down2)
    return NP.array([res1, res2])


def ellipse_angle_of_rotation2( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return NP.pi/2
    else:
        if a > c:
            return NP.arctan(2*b/(a-c))/2
        else:
            return NP.pi/2 + NP.arctan(2*b/(a-c))/2


if __name__ == '__main__':
    arc = 0.8
    R = NP.arange(0,arc*NP.pi, 0.01)
    x = 1.5*NP.cos(R) + 2 + 0.1*NP.random.rand(len(R))
    y = NP.sin(R) + 1. + 0.1*NP.random.rand(len(R))

    a = fit_ellipse(x,y)
    center = ellipse_center(a)
    #phi = ellipse_angle_of_rotation(a)
    phi = ellipse_angle_of_rotation2(a)
    axes = ellipse_axis_length(a)

    print('center = {}'.format(center))
    print('angle of rotation = {}'.format(phi))
    print('axes = {}'.format(axes))

    a, b = axes
    xx = center[0] + a*NP.cos(R)*NP.cos(phi) - b*NP.sin(R)*NP.sin(phi)
    yy = center[1] + a*NP.cos(R)*NP.sin(phi) + b*NP.sin(R)*NP.cos(phi)

    import pylab as PL
    from matplotlib.patches import Ellipse

    ellipse = Ellipse(center,
                      2 * axes[0],
                      2 * axes[1],
                      math.degrees(phi),
                      color='g',
                      alpha=0.2)

    fig = PL.figure()
    ax = PL.subplot(111)
    PL.scatter(x, y, marker='x')
    PL.plot(xx,yy, color = 'red')
    ax.add_artist(ellipse)
    PL.show()

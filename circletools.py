#!/usr/bin/env python3

import numpy as np
from numpy import exp, pi, sqrt, sin, cos, log, arccosh
from collections import defaultdict

def circle(T) :
    z1,z2,z3 = T
    ax, ay = z1.real, z1.imag
    bx, by = z2.real, z2.imag
    cx, cy = z3.real, z3.imag
    a = bx - ax
    b = by - ay
    c = cx - ax
    d = cy - ay
    e = a * (ax + bx) + b*(ay + by)
    f = c * (ax + cx) + d*(ay + cy)
    g = 2 * (a * (cy - by) - b * (cx - bx))
    if g == 0. :
        return (z1, z2)
    else :
        center = (d*e - b*f)/g + 1j *(a*f - c*e)/g
        radius = sqrt((ax - center.real)**2 + (ay - center.imag)**2)
        return (center, radius)

def shear_graft(s,x,y,z) :
    # Multiplicative shear s and graft pi/2 accross (x,y) edge.
    # Returns 4th point, opposite z.
    return ((1j + s)*x*y - 1j*x*z - s*y*z)/(s*x + 1j*y - (1j + s)*z)

def shear(s,x,y,z) :
    # Multiplicative shear s accross (x,y) edge.
    # Returns 4th point, opposite z.
    return ((1 + s)*x*y - (x + s*y)*z)/(s*x + y - (1 + s)*z)

def L(T,s) :
    return (shear(s,T[0],T[1],T[2]),T[1],T[0])

def R(T,s) :
    return (shear(s,T[2],T[0],T[1]),T[0],T[2])

def LGr(T,s) :
    return (shear_graft(s,T[0],T[1],T[2]),T[1],T[0])

def RGr(T,s) :
    return (shear_graft(s,T[2],T[0],T[1]),T[0],T[2])

def rot_L(T) :
    return (T[2],T[0],T[1])

def rot_R(T) :
    return (T[1], T[2], T[0])

def shear_up(s,t) :
    return [t, (-1+2*(s-t)*t)/(s-2*t+2*s*t**2), (s-2*t+2*s*t**2)/(-1+2*(s-t)*t), t, 1]

def shear_right(s,t) :
    return [1/t, s, 1/s, 1/t, 1]

def trace_one(s,t) :
    return (-1j - s + 2*t - 2*(s - 1j)*t**2)/(s*t)

def trace_two(s,t) :
    return (-1j - s + 2*1j*s*t - 4*1j*s*t**3 - 4*(s - 1j)*t**4)/(t*(s - 2*t + 2*s*t**2))

def conj_param_squared(s,t) :
    return (s - 2*1j*s*t**2 + 1j*(1 + (1 + 1j)*t)**2)/(1 - 1j*s - (2 - 2*1j)*t + 2*(s - 1j)*t**2)

# We can assume from the packing combinatorics that |Im(c)| <= Pi and
# |Im(c tau) | <= Pi because the big circle is self tangent so its translates
# can't go too far away. The equations below are adapted to that the imaginary
# satify the restrictions.

def c_param(s,t) :
    return 2*arccosh(trace_one(s,t)/2)

# unused
def c_tau_param(s,t) :
    return 2*arccosh(trace_two(s,t)/2)

def tau_param(s,t) :
    # slightly different from Mathematica. Possibly an difference in arccosh implementation
    if t > 1/sqrt(2) :
        return -arccosh(trace_two(s,t)/2)/arccosh(trace_one(s,t)/2)
    else :
        return arccosh(trace_two(s,t)/2)/arccosh(trace_one(s,t)/2)

def fern(W,H, sR, sU) :
    # 2W and 2H will be the total grid of the 4 pt circle
    # s_ are shear coords seqeunces for right and up generators.
    # we start with the first key triangles, the first one will always be the top one.
    # half coordinates mean intermediate triangles, integer are pt circles
    # always "top, bottom" triangle
    all_T = defaultdict(list)
    all_T[(0,0)] = [(exp(3*pi*1j/4), exp(-3*pi*1j/4), exp(pi*1j/4)),
                    (exp(-pi*1j/4), exp(pi*1j/4),exp(-3*pi*1j/4))]
    paths = {('vertical',1) : [RGr,RGr,LGr,RGr,L], ('vertical',-1) : [RGr,LGr,RGr,LGr,L],
            ('horizontal',1) : [LGr,LGr,RGr,LGr,R], ('horizontal',-1) : [LGr,RGr,LGr,RGr,R]}
    s_vertical = {1 : sU, -1 : sU[-2::-1] + [sU[-1]], 1/2 : sU[0], -1/2 : sU[-2]}
    s_horizontal = {1 : sR, -1 : sR[-2::-1] + [sU[-1]]}
    
    for sgn in [-1,1] :
        path = paths[('vertical',sgn)]
        shears = s_vertical[sgn]
        for m in range(H) :
            T = all_T[(0, m*sgn)][(1-sgn)//2]
            new = []
            for k in range(len(path)) :
                T = path[k](T, shears[k])
                if k != 1 :
                    new.append(T)
            all_T[(0,(m + 1/2)*sgn)] = new[:3]
            all_T[(0,(m + 1)*sgn)] = [rot_L(new[-2]), new[-1]][::-sgn]

    for m in range(-H,H+1) :
        for sgn in [-1,1] :
            path = paths[('horizontal',sgn)]
            shears = s_horizontal[sgn]
            for n in range(W) :
                T = all_T[(n*sgn,m)][(1+sgn)//2]
                new = []
                for f,s in zip(path,shears) :
                    T = f(T,s)
                    new.append(T)
                all_T[((n + 1/2)*sgn, m)] = new[:3]
                all_T[((n + 1)*sgn, m)] = [rot_R(new[-2]), new[-1]][::sgn]
                for ud in [-1,1] :
                    T = all_T[((n + 1)*sgn, m)][(1-ud)//2]
                    all_T[((n + 1)*sgn, m + ud/2)].append(RGr(T, s_vertical[ud/2]))

    return all_T

def all_circles(all_T, include_dual, x = None) :
    all_circles = {'small' : [], 'large' : [], 'dual' : []}
    for p, v in all_T.items() :
        x_whole = (p[0] - int(p[0]) == 0)
        y_whole = (p[1] - int(p[1]) == 0)
        if x is not None :
            v = tuple(map(lambda z : (x+z)/(x-z), v))
        if x_whole and y_whole :
            all_circles['small'].append(circle(v[0]))
        else :
            if y_whole and not x_whole :
                all_circles['large'].append(circle(v[1]))
                if include_dual :
                    all_circles['dual'].extend((circle(v[0]),circle(v[2])))
            elif include_dual :
                all_circles['dual'].extend(map(circle,v))

    return all_circles

from scipy.linalg import solve_banded

def bezier_control_points(data) :
    # Assume data is a numpy arry of shape (k,2)
    count = len(data) - 1
    if count < 1 :
        raise Exception("Requires at least two points")
    
    # Our equations are  :
    # 2 P1_0 + P1_1 = D_0 + 2 D_1
    # P1_{i-1} + 4 P1_i + 1 P1_{i+1} = 4 D_i + 2 D_{i+1}
    # 2 P1_{n-2} + 7 P1_{n-1} = 8 D_{n-1} + D_n
    
    diag = np.tile([1,4,1],(count,1))
    diag[0] = [0,2,1]
    diag[-1] = [1,7,0]
    diag[-2,2] = 2
    
    rhs = 4 * data[:count] + 2 * data[1:]
    rhs[0] -= 3 * data[0]
    rhs[-1] += 4 * data[count-1] - data[count]
    
    P1 = solve_banded((1,1), diag.transpose(), rhs,
                      overwrite_ab = True, overwrite_b = True,
                      check_finite = False)
    # Equations for P2
    # P2_i = 2 D_{i+1} - P1_{i+1}
    # P_{n-1} = (D_n + P1_{n-1})/2
    P2 = np.zeros((count,2))
    P2[:-1] = 2 * data[1:-1] - P1[1:]
    P2[-1] = (data[-1] + P1[-1])/2
                      
    return (P1, P2)

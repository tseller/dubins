# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 16:35:33 2013

@author: tim
"""

import numpy

def perp(v):
    #if not 0, return the unit perpendicular vector in the counterclockwise direction
    return numpy.array([-v[1], v[0]])/numpy.linalg.norm(v) if numpy.linalg.norm(v) else v

def normalize(v):
    return v / numpy.linalg.norm(v) if numpy.linalg.norm(v) else v
    

class Circle:
    def __init__(self, _c, _r, _o):
        self.center = _c
        self.radius = _r
        self.orientation = _o
        
    def print_data(self):
        print 'center: %s, radius: %s, orientation: %s' % (self.center, self.radius, self.orientation)
        

def calcArcLength(p,q,radius, orientation):
    gamma = numpy.dot(p,q) / (numpy.linalg.norm(p) * (numpy.linalg.norm(q)))

    theta = numpy.arccos(gamma)     

    if gamma > 1:
        theta = 0
    if gamma < -1:
        theta = numpy.pi

    return radius * theta if numpy.cross(p,q) * orientation >= 0 else radius*(2 * numpy.pi - theta)

def calcDubinsPaths(T0, T1, r):
    dubins_length = float('inf')
    
    c = []
    d = []
    # T0 and T1 are points in the tangent bundle.
    # thus, they have the form (p,x) where p is a point in the plane
    # and x is a tangent vector

    #c[0] and c[1] are the centers of the circles, radius r,
    # to which T0['x'] is tangent. Similarly for d[0] and d[1].
    c.append(Circle(T0['p'] + r * perp(T0['x']), r, 1))
    c.append(Circle(T0['p'] - r * perp(T0['x']), r, -1))
    d.append(Circle(T1['p'] + r * perp(T1['x']), r, 1))
    d.append(Circle(T1['p'] - r * perp(T1['x']), r, -1))

    for i in c:
        for j in d:
            if i.orientation == j.orientation:
                # if the centers are coincident, there will be just one arc.
                # so, let the SR (or SL) components be degenerate.
                if numpy.linalg.norm(i.center - j.center) == 0:
                    # L, R
                    int1 = T1['p'] - i.center
                    int2 = T1['p'] - j.center
                    
                    s = calcArcLength(T0['p'] - i.center, T1['p'] - i.center, r, i.orientation)
                    path = '%s for %s' % ('L' if i.orientation == 1 else 'R', s)
 
                    if s < dubins_length:
                        dubins_length = s
                        dubins_path = path

                    #print path, s
                else:
                    # RSR, LSL, S
                    int1 = i.center - i.orientation * r * perp(j.center - i.center)
                    int2 = j.center - i.orientation * r * perp(j.center - i.center)

                    s1 = calcArcLength(T0['p'] - i.center, int1 - i.center, r, i.orientation)
                    path = '%s for %s' %('L' if i.orientation == 1 else 'R', s1)
                    s2 = numpy.linalg.norm(i.center - j.center)
                    path += ', S for %s' %(s2)
                    s3 = calcArcLength(int2 - j.center, T1['p'] - j.center, r, j.orientation)
                    path += ', %s for %s' %('L' if j.orientation == 1 else 'R', s3)
                    
                    if s1 + s2 + s3 < dubins_length:
                        dubins_length = s1 + s2 + s3
                        dubins_path = path

                    #print path, s1 + s2 + s3
                                   
                    # RLR, LRL
                    if numpy.linalg.norm(i.center - j.center) < 4 * r:
                        e = Circle((i.center+j.center)/2 +  i.orientation * numpy.sqrt(4 * numpy.square(r) - numpy.square(numpy.linalg.norm((j.center-i.center)/2))) * perp(j.center-i.center), r, -1 * i.orientation)
                
                        int1 = (i.center + e.center)/2
                        int2 = (e.center + j.center)/2

                        s1 = calcArcLength(T0['p'] - i.center, int1 - i.center, r, i.orientation)
                        path = '%s for %s' %('L' if i.orientation == 1 else 'R', s1)
                        s2 = calcArcLength(int1 - e.center, int2 - e.center, r, -i.orientation)
                        path += ', %s for %s' %('L' if e.orientation == 1 else 'R', s2)
                        s3 = calcArcLength(int2 - j.center, T1['p'] - j.center, r, j.orientation)
                        path += ', %s for %s' %('L' if j.orientation == 1 else 'R', s3)

                        if s1 + s2 + s3 < dubins_length:
                            dubins_length = s1 + s2 + s3
                            dubins_path = path

                        #print path, s1 + s2 + s3
            # RSL, LSR
            else:
                if numpy.linalg.norm(i.center - j.center) >= 2 * r:
                    c_perp = r * numpy.sqrt(numpy.square(numpy.linalg.norm(j.center - i.center)/2) - numpy.square(r)) / numpy.linalg.norm((j.center - i.center)/2)
                    c_parallel = numpy.linalg.norm((j.center - i.center)/2) - numpy.square(r) / numpy.linalg.norm((j.center - i.center)/2)

                    int1 = (i.center + j.center)/2 - i.orientation * c_perp * perp(j.center - i.center) - c_parallel * normalize(j.center - i.center)
                    int2 = (i.center + j.center)/2 + i.orientation * c_perp * perp(j.center - i.center) + c_parallel * normalize(j.center - i.center)
                    
                    s1 = calcArcLength(T0['p'] - i.center, int1 - i.center, r, i.orientation)
                    path = '%s for %s' %('L' if i.orientation == 1 else 'R', s1)
                    s2 = numpy.linalg.norm(int1 - int2)
                    path += ', S for %s' %(s2)
                    s3 = calcArcLength(int2 - j.center, T1['p'] - j.center, r, j.orientation)
                    path += ', %s for %s' %('L' if j.orientation == 1 else 'R', s3)

                    if s1 + s2 + s3 < dubins_length:
                        dubins_length = s1 + s2 + s3
                        dubins_path = path

                    #print path, s1 + s2 + s3

    return dubins_length
    
#sample input
'''
T0 = {'p':numpy.array([0,0]), 'x':numpy.array([1,0])}
T1 = {'p':numpy.array([0,0]), 'x':numpy.array([0,1])}
print calcDubinsPaths(T0, T1, 1)
'''    
#!/usr/bin/env python
from iges.entity import Entity
import os

class Line(Entity):
    """Straight line segment"""

    def add_parameters(self, parameters):
        self.x1 = float(parameters[1])
        self.y1 = float(parameters[2])
        self.z1 = float(parameters[3])
        self.x2 = float(parameters[4])
        self.y2 = float(parameters[5])
        self.z2 = float(parameters[6])

    def __str__(self):
        s = '--- Line ---' + os.linesep
        s += Entity.__str__(self) + os.linesep
        s += "From point {0}, {1}, {2} {3}".format(self.x1, self.y1, self.z1, os.linesep)
        s += "To point {0}, {1}, {2}".format(self.x2, self.y2, self.z2)
        return s

class CompositeCurveEntity(Entity):
    """ Composite Curve
    iges spec v5.3 p. 69 section 4.4
    """

    def add_parameters(self, parameters):
        self.N = int(parameters[1])
        
        self.DE = []

        for i in range(self.N):
            self.DE.append(int(parameters[i + 2]))
            
    def __str__(self):
        s = '--- Composite Curve ---' + os.linesep
        # s += Entity.__str__(self) + os.linesep
        s += str(self.N) + os.linesep
        s += str(self.DE)

        return s
        
    
        
class RationalBSplineCurve(Entity):
    """rational b-spline curve
    iges spec v5.3 p. 123 section 4.23
    see also appendix b, p. 545
    """

    def add_parameters(self, parameters):
        self.K = int(parameters[1])
        self.M = int(parameters[2])
        self.prop1 = int(parameters[3]) # planar
        self.prop2 = int(parameters[4]) # closed
        self.prop3 = int(parameters[5]) # polynonial
        self.prop4 = int(parameters[6]) # periodic
        
        self.N = 1 + self.K - self.M
        self.A = self.N + 2 * self.M

        # Knot sequence
        self.T = []
        for i in range(7, 7 + self.A + 1):
            self.T.append(float(parameters[i]))

        # Weights
        self.W = []
        for i in range(self.A + 8, self.A + self.K + 8 + 1):
            self.W.append(float(parameters[i]))

        # Control points
        self.control_points = []
        for i in range(9 + self.A + self.K, 9 + self.A + 4*self.K + 1, 3):
            point = (float(parameters[i]), float(parameters[i+1]), float(parameters[i+2]))
            self.control_points.append(point)

        # Parameter values
        self.V0 = float(parameters[12 + self.A + 4 * self.K])
        self.V1 = float(parameters[13 + self.A + 4 * self.K])

        # Unit normal (only for planar curves)
        if len(parameters) > 14 + self.A + 4 * self.K + 1:
            self.planar_curve = True
            self.XNORM = float(parameters[14 + self.A + 4 * self.K])
            self.YNORM = float(parameters[15 + self.A + 4 * self.K])
            self.ZNORM = float(parameters[16 + self.A + 4 * self.K])
        else:
            self.planar_curve = False

    def __str__(self):
        s = '--- Rational B-Spline Curve ---' + os.linesep
        s += Entity.__str__(self) + os.linesep
        s += "K: {}".format(self.K) + os.linesep
        s += "degree: {}".format(self.M) + os.linesep
        s += "pieces: {}".format(self.N) + os.linesep
        s += "ncpt: {}".format(len(self.control_points)) + os.linesep
        s += "knots: {}".format(str(self.T)) + os.linesep
        s += str(self.W) + os.linesep
        s += str(self.control_points) + os.linesep
        s += "Parameter: v(0) = {0}    v(1) = {1}".format(self.V0, self.V1) + os.linesep
        if self.planar_curve:
            s += "Unit normal: {0} {1} {2}".format(self.XNORM, self.YNORM, self.ZNORM)
        return s


class RationalBSplineSurface(Entity):
    """Rational B-Spline Surface
    IGES Spec v5.3 p. 126 Section 4.24
    """

    def add_parameters(self, parameters):
        self.K1 = int(parameters[1])
        self.K2 = int(parameters[2])
        self.M1 = int(parameters[3])
        self.M2 = int(parameters[4])
        
        self.prop1 = int(parameters[5])
        self.prop2 = int(parameters[6])
        self.prop3 = int(parameters[7])
        self.prop4 = int(parameters[8])
        self.prop5 = int(parameters[9])
        
        self.N1 = 1 + self.K1 - self.M1
        self.N2 = 1 + self.K2 - self.M2
        
        self.A = self.N1 + 2 * self.M1
        self.B = self.N2 + 2 * self.M2
        self.C = (1 + self.K1) * (1 + self.K2)

        # Knot sequence
        self.T1 = []
        for i in range(10, 10 + self.A + 1):
            self.T1.append(float(parameters[i]))
            
        self.T2 = []
        for i in range(11 + self.A, 11 + self.A + self.B + 1):
            self.T2.append(float(parameters[i]))

        # Weights
        # 1+K1 * 1+K2
        self.W = []
        for i in range(self.A + self.B + 12, 11 + self.A + self.B + self.C + 1):
            self.W.append(float(parameters[i]))

        # Control points
        self.control_points = []
        for i in range(12 + self.A + self.B + self.C, 9 + self.A + self.B + 4*self.C + 1, 3):
            point = (float(parameters[i]), float(parameters[i+1]), float(parameters[i+2]))
            self.control_points.append(point)

        # Parameter values
        self.U0 = float(parameters[12 + self.A + self.B + 4 * self.C])
        self.U1 = float(parameters[13 + self.A + self.B + 4 * self.C])
                        
        self.V0 = float(parameters[14 + self.A + self.B + 4 * self.C])
        self.V1 = float(parameters[15 + self.A + self.B + 4 * self.C])
                        
        # # Unit normal (only for planar curves)
        # if len(parameters) > 14 + self.A + 4 * self.K + 1:
        #     sf.planar_curve = True
        #     self.XNORM = float(parameters[14 + self.A + 4 * self.K])
        #     self.YNORM = float(parameters[15 + self.A + 4 * self.K])
        #     self.ZNORM = float(parameters[16 + self.A + 4 * self.K])
        # else:
        #     self.planar_curve = False

    def __str__(self):
        s = '--- Rational B-Spline Surface ---' + os.linesep
        s += Entity.__str__(self) + os.linesep
        s += "knots U:" + os.linesep + str(self.T1) + os.linesep
        s += "knots V:" + os.linesep + str(self.T2) + os.linesep
        s += str(self.W) + os.linesep
        # s += str(self.control_points) + os.linesep
        s += "Degree: u - {}, v - {}".format(self.M1, self.M2) + os.linesep
        s += "Rational/Polynomial: {}".format(self.prop3) + os.linesep
        s += "Periodic: u - {}, v - {}".format(self.prop4, self.prop5) + os.linesep
        s += "Parameter: u(0) = {0}    u(1) = {1}".format(self.U0, self.U1) + os.linesep
        s += "Parameter: v(0) = {0}    v(1) = {1}".format(self.V0, self.V1) + os.linesep
        return s
        
class BoundaryEntity(Entity):
    """ Boundary Entity
    IGES Spec v5.3 p. 157 Section 4.31
    """

    def add_parameters(self, parameters):
        self.TYPE = int(parameters[1])
        self.PREF = int(parameters[2])
        self.SPTR = int(parameters[3])
        self.N = int(parameters[4])

        self.PSCPT = []

        curr_idx = 5
        
        for i in range(self.N):
            CRVPT = int(parameters[curr_idx])
            curr_idx += 1

            SENSE = int(parameters[curr_idx])
            curr_idx += 1
            
            num_curves = int(parameters[curr_idx])
            curr_idx += 1

            curves = []
            for j in range(num_curves):
                curves.append(int(parameters[curr_idx]))
                curr_idx += 1
                
            self.PSCPT.append((CRVPT, SENSE, curves))

        
    def __str__(self):
        s = '--- Boundary ---' + os.linesep
        # s += Entity.__str__(self) + os.linesep
        s += str(self.PSCPT)
        # s += str(self.control_points) + os.linesep
        return s
                
            

class BoundedSurfaceEntity(Entity):
    """ Bounded Surface
    IGES Spec v5.3 p. 165 Section 4.33
    """

    def add_parameters(self, parameters):
        self.TYPE = int(parameters[1])
        self.SPTR = int(parameters[2])
        self.N = int(parameters[3])

        self.BDPT = []

        for i in  range(4, 4 + self.N):
            self.BDPT.append(int(parameters[i]))
            

    def __str__(self):
        s = '--- Bounded Surface ---' + os.linesep
        # s += Entity.__str__(self) + os.linesep
        s += str(self.SPTR) + os.linesep
        s += str(self.BDPT)
        # s += str(self.control_points) + os.linesep
        return s
        
class ParametericCurveEntity(Entity):
    """ Curve on a parametrix surface
    IGES Spec v5.3 p. 163 Section 4.32
    """

    def add_parameters(self, parameters):
        self.CRTN = int(parameters[1])
        self.SPTR = int(parameters[2])
        self.BPTR = int(parameters[3])
        self.CPTR = int(parameters[4])
        self.PREF = int(parameters[5])


    def __str__(self):
        s = '--- Parameteric Curve ---' + os.linesep
        # s += Entity.__str__(self) + os.linesep
        s += "CRTN " + str(self.CRTN) + os.linesep
        s += "SPTR " + str(self.SPTR) + os.linesep
        s += "BPTR " + str(self.BPTR) + os.linesep
        s += "CPTR " + str(self.CPTR) + os.linesep
        s += "PREF " + str(self.PREF)
        # s += str(self.control_points) + os.linesep
        return s
    
class TrimmedSurfaceEntity(Entity):
    """ Trimmed Surface
    IGES Spec v5.3 p. 167 Section 4.34
    """

    def add_parameters(self, parameters):
        self.PTS = int(parameters[1])
        self.N1 = int(parameters[2])
        self.N2 = int(parameters[3])
        self.PTO = int(parameters[4])

        self.PTI = []
        for i in range(self.N2):
            self.PTI.append(int(parameters[i+5]))


    def __str__(self):
        s = '--- Trimmed Surface ---' + os.linesep
        # s += Entity.__str__(self) + os.linesep
        s += str(self.PTS) + os.linesep
        s += "N1: " + str(self.N1) + os.linesep
        s += "N2: " + str(self.N2) + os.linesep
        s += "PTO: " + str(self.PTO) + os.linesep
        s += "PTI: " + str(self.PTI)
        # s += str(self.control_points) + os.linesep
        return s

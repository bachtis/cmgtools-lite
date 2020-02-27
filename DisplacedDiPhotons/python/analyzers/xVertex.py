import math
import ROOT
from sympy import Symbol, nsolve

class xVertex(object):
    def __init__(self, v1, v2, e1, e2, mass):
        self.v1 = v1
        self.v2 = v2
        self.e1 = e1
        self.e2 = e2
        self.mass = mass

    # Given two ecal photon positions, return angle between photon plane and xy plane
    def getRotAngle(self):
        v1 = self.v1
        v3 = v1.Cross(self.v2)
        v3.SetMag(1.0)
        return math.acos(v3[2])

    # Get axis of rotation for photon - xy plane rotation
    def getRotAxis(self):
        v1 = self.v1
        v3 = v1.Cross(self.v2)
        v3.SetMag(1.0)
        ez = ROOT.TVector3(0,0,1.0)
        return v3.Cross(ez)

    # Given mass and two ecal energies, find inscribed angle phi
    def getPhi(self):
        m = self.mass
        e1 = self.e1
        e2 = self.e2
        return math.acos(1.0-m**2/(2*e1*e2))

    ########FOLLOWING METHODS ARE FOR 2D WHEN GETTING VERTEX, AFTER WE ROTATE THE ECAL VECTORS###################

    # Get radius of circles given inscribed angle and two points
    def getRadius(self, x1,x2,y1,y2,phi):
        return math.sqrt(((x1-x2)/2)**2+((y1-y2)/2)**2)/math.sin(phi)

    # Get center of circles given inscribed angle and two points
    def getCenters(self, x1,x2,y1,y2,phi):
        r =self.getRadius(x1,x2,y1,y2,phi)
        q = math.sqrt((x2-x1)**2+(y2-y1)**2)
        x3 = (x1+x2)/2
        y3 = (y1+y2)/2
        xp = x3+math.sqrt(r**2-(q/2)**2)*(y1-y2)/q
        yp = y3+math.sqrt(r**2-(q/2)**2)*(x2-x1)/q
        xm = x3-math.sqrt(r**2-(q/2)**2)*(y1-y2)/q
        ym = y3-math.sqrt(r**2-(q/2)**2)*(x2-x1)/q
        return [[xp,yp],[xm,ym]]

    # Return pt of two points
    def getPt(self, x,y,x1,x2,y1,y2):
        e1 = self.e1
        e2 = self.e2
        cosTheta1 = (x*(x1-x)+y*(y1-y))/(math.sqrt(x**2+y**2)*math.sqrt((x1-x)**2+(y1-y)**2))
        cosTheta2 = (x*(x2-x)+y*(y2-y))/(math.sqrt(x**2+y**2)*math.sqrt((x2-x)**2+(y2-y)**2))
        theta1 = math.acos(cosTheta1)
        theta2 = math.acos(cosTheta2)
        pt1 = e1*math.sin(theta1)
        pt2 = e2*math.sin(theta2)
        return pt1+pt2

    # Get phi from vertex of (x,y) going to (x1,y1) (x2,y2) - to remove points directly between ecal hits
    def getPhiFromPoints(self, x,y,x1,x2,y1,y2):
        cosTheta = ((x1-x)*(x2-x)+(y1-y)*(y2-y))/(math.sqrt((x1-x)**2+(y1-y)**2)*math.sqrt((x2-x)**2+(y2-y)**2))
        theta = math.acos(cosTheta)
        return theta


    # Put everything together, return [best point, pt at that point]
    def getVertex(self):
        v1 = self.v1
        v2 = self.v2
        axis = self.getRotAxis()
        theta = self.getRotAngle()
        v1.Rotate(theta, axis)
        v2.Rotate(theta, axis)
        phi = self.getPhi()
        radius = self.getRadius(v1[0], v2[0], v1[1], v2[1], phi)
        centers = self.getCenters(v1[0], v2[0], v1[1], v2[1], phi)
        steps = 1000
        delta = 2*math.pi/steps
        points = []
        for c in centers:
            for i in range(steps):
                angle = i*delta
                x = c[0]+radius*math.cos(angle)
                y = c[1]+radius*math.sin(angle)
                if abs(self.getPhiFromPoints(x,y,v1[0],v2[0],v1[1],v2[1]) - phi) > 0:
                    continue
                points.append([[x,y],self.getPt(x,y,v1[0],v2[0],v1[1],v2[1])])
        if len(points) == 0:
            return [ROOT.TVector3(-999,-999,-999),999]
        best = min(points, key = lambda x: abs(x[1]))
        coord = ROOT.TVector3(best[0][0], best[0][1], 0)
        coord.Rotate(-theta, axis)
        return [coord, best[1]]

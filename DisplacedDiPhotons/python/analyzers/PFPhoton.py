import math
import ROOT
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi

class PFPhoton(object):
    def __init__(self, cand):
        self.packedCand = cand
        self.LV = cand.p4()
        self.pdgId = cand.pdgId()
        self.caloPosition = self.propagateToCalo()

    def p4(self):
        return self.LV

    def pdgId(self):
        return self.pdgId

    def m(self):
        return self.LV.mass()

    def propagateToCalo(self):
        vertex = self.packedCand.vertex()
        x,y,z = 0,0,0
        vx = vertex.x()
        vy = vertex.y()
        vz = vertex.z()
        px = self.LV.px()
        py = self.LV.py()
        pz = self.LV.pz()
        pt = self.LV.pt()
        phi = self.LV.phi()
        # If endcap, propagate to disk at +-z=324 cm
        if abs(self.LV.eta()) > 1.442:
            z = math.copysign(324., pz)
            x = vx + (z-vz)*(px/pz)
            y = vy + (z-vy)*(py/pz)
        # If barrel, propagate to cylinder with r=137 cm
        # To first approx assume only dz significant
        else:
            x = 137.*math.cos(phi)
            y = 137.*math.sin(phi)
            z = vz + 137.*(pz/pt)
        return ROOT.TVector3(x,y,z)
        
    def __getattr__(self, name):
        return getattr(self.LV, name)

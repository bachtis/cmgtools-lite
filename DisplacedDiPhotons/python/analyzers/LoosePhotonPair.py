import math
import ROOT
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from CMGTools.DisplacedDiPhotons.analyzers.xVertex import *

class LoosePhotonPair(object):
    def __init__(self, leg1, leg2, pdg = 0):
        self.leg1 = leg1
        self.leg2 = leg2
        self.pdg = pdg
        self.LV = leg1.p4() + leg2.p4()
        self.vertex10 = XVertex(leg1.caloPosition, leg2.caloPosition, leg1.energy, leg2.energy, 10)
        self.vertex20 = XVertex(leg1.caloPosition, leg2.caloPosition, leg1.energy, leg2.energy, 20)
        self.vertex30 = XVertex(leg1.caloPosition, leg2.caloPosition, leg1.energy, leg2.energy, 30)
        self.vertex40 = XVertex(leg1.caloPosition, leg2.caloPosition, leg1.energy, leg2.energy, 40)
        self.vertex50 = XVertex(leg1.caloPosition, leg2.caloPosition, leg1.energy, leg2.energy, 50)
        self.vertex60 = XVertex(leg1.caloPosition, leg2.caloPosition, leg1.energy, leg2.energy, 60)

    def p4(self):
        return self.LV

    def m(self):
        return self.LV.mass()

    def pdgId(self):
        return self.pdg

    def deltaPhi(self):
        return abs(deltaPhi(self.leg1.phi(), self.leg2.phi()))

    def deltaR(self):
        return abs(deltaR(self.leg1.eta(), self.leg1.phi(), self.leg2.eta(), self.leg2.phi()))

    def __getattr__(self, name):
        return getattr(self.LV, name)

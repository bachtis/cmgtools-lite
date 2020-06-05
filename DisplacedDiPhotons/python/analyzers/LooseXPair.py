import math
import ROOT
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi

class XPair(object):
    def __init__(self, x1, x2):
        self.x1 = x1
        self.x2 = x2
        self.LV = x1.p4()+x2.p4()

    def rawP4(self):
        return self.x1.p4()+self.x2.p4()

    def p4(self):
        return self.LV

    def m(self):
        return self.LV.mass()

    def m12(self):
        return self.x1.m()

    def m34(self):
        return self.x2.m()

    def m13(self):
        return (self.x1.leg1.p4()+self.x2.leg1.p4()).mass()

    def m14(self):
        return (self.x1.leg1.p4()+self.x2.leg2.p4()).mass()

    def m23(self):
        return (self.x1.leg2.p4()+self.x2.leg1.p4()).mass()

    def m24(self):
        return (self.x1.leg2.p4()+self.x2.leg2.p4()).mass()

    def deltaM(self):
        return abs(self.x1.m()-self.x2.m())

    def __getattr__(self, name):
        return getattr(self.LV, name)

class ZXX(object):
    def __init__(self, Z, XPair):
        self.Z = Z
        self.XX = XPair
        self.LV = Z.p4()+XPair.p4()

    def rawP4(self):
        return self.Z.p4() + self.XX.p4()

    def p4(self):
        return self.LV

    def  __getattr__(self, name):
        return getattr(self.LV, name)


class WXX(object):
    def __init__(self, W, XPair):
        self.W = W
        self.XX = XPair
        self.LV = W.p4() + XPair.p4()

    def rawP4(self):
        return self.W.p4() + self.XX.p4()

    def p4(self):
        return self.LV

    def  __getattr__(self, name):
        return getattr(self.LV, name)

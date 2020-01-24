from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import *


import ROOT
import os
import math

import  itertools
from CMGTools.VVResonances.tools.Pair import *
from CMGTools.DisplacedDiPhotons.analyzers.PhotonPair import *
from CMGTools.DisplacedDiPhotons.analyzers.XPair import *

class VHGGBuilder(Analyzer):

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(VHGGBuilder, self).__init__(cfg_ana, cfg_comp, looperName)


    def declareHandles(self):
        super(VHGGBuilder, self).declareHandles()

    def makeW(self,leptons,met):
        Ws=[]
        # Pair leptons to met
        for l in leptons:
            W=Pair(l,met,l.charge()*24)
            Ws.append(W)
        return Ws

    def makeZ(self,leptons):
        Zs=[]
        # Pair leptons if near Z mass
        for l1,l2 in itertools.combinations(leptons,2):
            if (l1.charge()+l2.charge())!=0:
                continue
            Z = Pair(l1,l2,23)
            mass = Z.p4().mass()
            if mass<50 or mass>140:
                continue
            Zs.append(Z)
        return Zs

    def makeX(self,photons):
        Xs=[]
        # Pair photons (TODO: Conversions)
        for g1,g2 in itertools.combinations(photons,2):

            X = PhotonPair(g1,g2,1995)
            Xs.append(X)
            
        return Xs

    def makeXPair(self, photons):
        XXs = []
        if len(photons) < 4:
            return []
        for g1,g2,g3,g4 in itertools.combinations(photons, 4):
            X12 = PhotonPair(g1, g2)
            X34 = PhotonPair(g3, g4)
            X13 = PhotonPair(g1, g3)
            X24 = PhotonPair(g2, g4)
            X23 = PhotonPair(g2, g3)
            X14 = PhotonPair(g1, g4)
            
            X1234 = XPair(X12, X34)
            X1324 = XPair(X13, X24)
            X1423 = XPair(X14, X23)
            XXs.append(X1234)
            XXs.append(X1324)
            XXs.append(X1423)
        return XXs

    def process(self, event):
        self.readCollections(event.input)
        event.Z = []
        event.ZX=[]
        event.WX=[]
        event.ZXX = []
        event.WXX = []
        goodLeptons = filter(lambda x: x.pt()>0,event.selectedLeptons)
        goodPhotons = filter(lambda x: x.pt()>0,event.selectedPhotons)

        Zs = self.makeZ(goodLeptons)
        Ws = self.makeW(goodLeptons,event.met)
        Xs = self.makeX(goodPhotons)
        XXs = self.makeXPair(goodPhotons)

        if len(Zs)>0:
            bestZ = max(Zs,key=lambda x: x.leg1.pt()+x.leg2.pt())
            event.Z.append(bestZ)
            goodXs=[]
            for X in Xs:
                overlap=False
                for l in [bestZ.leg1,bestZ.leg2]:
                    for g in [X.leg1,X.leg2]:
                        if deltaR(l.eta(),l.phi(),g.eta(),g.phi())<0.5:
                            overlap=True
                            break
                if not overlap:
                    goodXs.append(X)
            
            if len(goodXs)>0:
                bestX = max(goodXs,key=lambda x: x.leg1.pt()+x.leg2.pt())
                bestZX = Pair(bestZ,bestX)
                bestZX.otherLeptons = len(goodLeptons)-2
                event.ZX.append(bestZX)

            goodXXs = []
            for XX in XXs:
                overlap = False
                for l in [bestZ.leg1, bestZ.leg2]:
                    for X in XX:
                        for g in [X.leg1, X.leg2]:
                            if deltaR(l.eta(), l.phi(), g.eta(), g.phi()) < 0.5:
                                overlap = True
                                break
                if not overlap:
                    goodXXs.append(XX)
            if len(goodXXs) > 0:
                bestXX = min(goodXXs, key = lambda x: x.deltaM())
                bestZXX = ZXX(bestZ, bestXX)
                event.ZXX.append(bestZXX)

        elif len(Ws) > 0:
            bestW = max(Ws, key = lambda x: x.leg1.pt())
            goodXs = []
            for X in Xs:
                for W in Ws:
                    overlap = False
                    for g in [X.leg1, X.leg2]:
                        if deltaR(W.leg1.eta(), W.leg1.phi(), g.eta(), g.phi()) < 0.5:
                            overlap = True
                            break
                if not overlap:
                    goodXs.append(X)
                        
            if len(goodXs) > 0:
                bestX = max(goodXs, key=lambda x: x.leg1.pt() + x.leg2.pt())
                bestWX = Pair(bestW, bestX)
                event.WX.append(bestWX)
            goodXXs = []
            for XX in XXs:
                overlap = False
                for W in Ws:
                    for X in XX:
                        for g in [X.leg1, X.leg2]:
                            if deltaR(W.leg1.eta(), W.leg1.phi(), g.eta(), g.phi()) < 0.5:
                                overlap = True
                                break
                goodXXs.append(XX)

            if len(goodXXs) > 0:
                bestXX = min(goodXXs, key = lambda x: x.deltaM())
                bestWXX = WXX(bestW, bestXX)
                event.WXX.append(bestWxx)

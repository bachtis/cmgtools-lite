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
        # Loop over 4 photons to make XXs
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

    def checkFSR_ZX(self, Z,X):
        flag = 0
        mz = Z.p4().mass()
        for g in [X.leg1, X.leg2]:
            p4New = Z.p4() + g.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91):
                flag = 1

        p4Total = Z.p4() + X.p4()
        if abs(p4Total.mass()-91) < abs(mz-91):
            flag = 1

        return flag

    def checkFSR_ZXX(self, Z,XX):
        flag = 0
        mz = Z.p4().mass()
        for g in [XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2]:
            p4New = Z.p4() + g.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91):
                flag =  1
        for g1, g2 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],2):
            p4New = Z.p4() + g1.p4(2) + g2.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91):
                flag = 1
        for g1, g2, g3 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],3):
            p4New = Z4.p4() + g1.p4(2) + g2.p4(2) + g3.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91):
                flag = 1
        for g1, g2, g3, g4 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],4):
            p4New = Z4.p4() + g1.p4(2) + g2.p4() + g3.p4(2) + g4.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91):
                flag = 1

        return flag

    def process(self, event):
        self.readCollections(event.input)

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


        # Make ZX/ZXX Pairs first
        nZPairs = 0
        if len(Zs)>0:
            bestZ = max(Zs,key=lambda x: x.leg1.pt()+x.leg2.pt())

            goodXs=[]
            # Select non overlapping Xs
            for X in Xs:
                overlap=False
                for l in [bestZ.leg1,bestZ.leg2]:
                    for g in [X.leg1,X.leg2]:
                        if deltaR(l.eta(),l.phi(),g.eta(),g.phi())<0.5:
                            overlap=True
                            break
                if not overlap:
                    goodXs.append(X)

            # Pair Z to best X
            if len(goodXs)>0:
                bestX = max(goodXs,key=lambda x: x.leg1.pt()+x.leg2.pt())
                bestZX = Pair(bestZ,bestX)
                bestZX.otherLeptons = len(goodLeptons)-2
                bestZX.hasFSR = self.checkFSR_ZX(bestZ, bestX)
                bestZX.deltaPhi_g1 = min(abs(deltaPhi(bestZ.leg1.phi(),bestX.leg1.phi())),abs(deltaPhi(bestZ.leg2.phi(),bestX.leg1.phi())))
                bestZX.deltaPhi_g2 = min(abs(deltaPhi(bestZ.leg1.phi(),bestX.leg2.phi())),abs(deltaPhi(bestZ.leg2.phi(),bestX.leg2.phi())))
                event.ZX.append(bestZX)
                nZPairs+=1

            # Repeat for XX pairs
            goodXXs = []
            for XX in XXs:
                overlap = False
                for l in [bestZ.leg1, bestZ.leg2]:
                    for X in [XX.x1, XX.x2]:
                        for g in [X.leg1, X.leg2]:
                            if deltaR(l.eta(), l.phi(), g.eta(), g.phi()) < 0.5:
                                overlap = True
                                break
                if not overlap:
                    goodXXs.append(XX)
            #Pair best XX to best Z
            if len(goodXXs) > 0:
                bestXX = max(goodXXs, key = lambda x: x.x1.pt()+x.x2.pt())
                bestZXX = ZXX(bestZ, bestXX)
                bestZXX.otherLeptons = len(goodLeptons) - 2
                bestZXX.hasFSR = self.checkFSR_ZXX(bestZ, bestXX)
                event.ZXX.append(bestZXX)
                nZPairs+=1


        # If no Zs, search for Ws
        if len(Ws) > 0 and nZPairs == 0:
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
                bestWX.otherLeptons = len(goodLeptons) -1
                event.WX.append(bestWX)

            goodXXs = []
            for XX in XXs:
                overlap = False
                for W in Ws:
                    for X in [XX.x1, XX.x2]:
                        for g in [X.leg1, X.leg2]:
                            if deltaR(W.leg1.eta(), W.leg1.phi(), g.eta(), g.phi()) < 0.5:
                                overlap = True
                                break
                goodXXs.append(XX)

            if len(goodXXs) > 0:
                bestXX = max(goodXXs, key = lambda x: x.x1.pt()+x.x2.pt())
                bestWXX = WXX(bestW, bestXX)
                bestWXX.otherLeptons = len(goodLeptons) - 1
                event.WXX.append(bestWXX)

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
from CMGTools.DisplacedDiPhotons.analyzers.xVertex import *

debug = False

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
            # Did this need to be added? Wasn't in and couldn't find any instances where it wasn't true
            if l1.pdgId()+l2.pdgId()!=0:
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
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1

        p4Total = Z.p4() + X.p4()
        if abs(p4Total.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
            flag = 1

        return flag

    def checkFSR_ZXX(self, Z,XX):
        flag = 0
        mz = Z.p4().mass()
        for g in [XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2]:
            p4New = Z.p4() + g.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag =  1
        for g1, g2 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],2):
            p4New = Z.p4() + g1.p4(2) + g2.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1
        for g1, g2, g3 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],3):
            p4New = Z.p4() + g1.p4(2) + g2.p4(2) + g3.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1
        for g1, g2, g3, g4 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],4):
            p4New = Z.p4() + g1.p4(2) + g2.p4(2) + g3.p4(2) + g4.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1

        return flag

        
    def log(self, signal, leptons):
        for S in signal:
            print "S vertex=({:.2f}, {:.2f}, {:.2f}), pt={:.2f}, eta={:.2f}, phi={:.2f}, ID={}".format(S.vx(),S.vy(),S.vz(),S.pt(),S.eta(),S.phi(), S.pdgId())
            for d in range(S.numberOfDaughters()):
                print "    Daughter vertex=({:.2f}, {:.2f}, {:.2f}), pt={:.2f}, eta={:.2f}, phi={:.2f}, ID={}".format(S.daughter(d).vx(),S.daughter(d).vy(),S.daughter(d).vz(),S.daughter(d).pt(),S.daughter(d).eta(),S.daughter(d).phi(), S.daughter(d).pdgId())
        for l in leptons:
            print "Lepton pt={}, eta={}, phi={}, ID={}".format(l.pt(), l.eta(), l.phi(), l.pdgId())

    def process(self, event):
        self.readCollections(event.input)

        event.ZX=[]
        event.WX=[]
        event.ZXX = []
        event.WXX = []
        leptons = filter(lambda x: x.pt() > 0, event.selectedLeptons)
        goodLeptons = filter(lambda x: x.relIso03 < .1, leptons)
        photons = filter(lambda x: x.pt()>0,event.selectedPhotons)
        gen = event.genParticles
        genPhotons = filter(lambda x: x.pdgId() == 22, gen)
        genWs = filter(lambda x: x.pdgId() == 24 and x.numberOfDaughters()==2, gen)
#        import pdb
#        pdb.set_trace()
        goodPhotons = []
        for x in photons:
            overlap = False
            for l in goodLeptons:
                if deltaR(l.eta(), l.phi(), x.eta(), x.phi()) < 0.3:
                    overlap = True
            if not overlap:
                goodPhotons.append(x)

        Zs = self.makeZ(goodLeptons)
        Ws = self.makeW(goodLeptons,event.met)
        Xs = self.makeX(goodPhotons)
        XXs = self.makeXPair(goodPhotons)
        genLeptons = filter(lambda x: abs(x.pdgId())==11 or abs(x.pdgId())==13, gen)
        signal = filter(lambda x: abs(x.pdgId())==9000006, gen)
        signalPhotons = filter(lambda x: abs(x.mother().pdgId())==9000006, genPhotons)
        event.GenPhoton = signalPhotons
        # Make ZX/ZXX Pairs first
        nZPairs = 0
        if len(Zs)>0:
            bestZ = max(Zs,key=lambda x: x.leg1.pt()+x.leg2.pt())
            goodXs = Xs
#            goodXs=[]
#             Select non overlapping Xs
#            for X in Xs:
#                overlap=False
#                for l in [bestZ.leg1,bestZ.leg2]:
#                    for g in [X.leg1,X.leg2]:
#                        if deltaR(l.eta(),l.phi(),g.eta(),g.phi())<0.5:
#                            overlap=True
#                            break
#                if not overlap:
#                    goodXs.append(X)

            # Pair Z to best X
            if len(goodXs)>0:
               # bestX = max(goodXs,key=lambda x: x.leg1.pt()+x.leg2.pt())
                bestX = min(goodXs, key = lambda x: (x.p4()+bestZ.p4()).pt())
                bestZX = Pair(bestZ,bestX)
                bestZX.otherLeptons = len(goodLeptons)-2
                bestZX.hasFSR = self.checkFSR_ZX(bestZ, bestX)
                bestZX.deltaPhi_g1 = min(abs(deltaPhi(bestZ.leg1.phi(),bestX.leg1.phi())),abs(deltaPhi(bestZ.leg2.phi(),bestX.leg1.phi())))
                bestZX.deltaPhi_g2 = min(abs(deltaPhi(bestZ.leg1.phi(),bestX.leg2.phi())),abs(deltaPhi(bestZ.leg2.phi(),bestX.leg2.phi())))
                event.ZX.append(bestZX)
                nZPairs+=1

                if debug:
                    self.log(signal, genLeptons)
                    print "ZX Vertex: ({:.2f}, {:.2f}, {:.2f}) valid={}".format(bestZX.leg2.vertex15.vertex[0],bestZX.leg2.vertex15.vertex[1],bestZX.leg2.vertex15.vertex[2], bestZX.leg2.vertex15.valid)
                    print "X Photons:\n   ",
                    print bestZX.leg2.leg1
                    print "   ",
                    print bestZX.leg2.leg2
                    print "Gen Photons"
                    for g in genPhotons:
                        print "    vertex=({:.2f}, {:.2f}, {:.2f}), pt={:.2f}, eta={:.2f}, phi={:.2f}, energy={:.2f}, mother={}".format(g.vx(),g.vy(),g.vz(),g.pt(),g.eta(),g.phi(),g.energy(),g.mother().pdgId())
                    print "Reco photons"
                    for p in photons:
                        print "    ",
                        print p
                    
                    import pdb
                    pdb.set_trace()

                
            goodXXs = XXs
            # Repeat for XX pairs
#            goodXXs = []
#            for XX in XXs:
#                overlap = False
#                for l in [bestZ.leg1, bestZ.leg2]:
#                    for X in [XX.x1, XX.x2]:
#                        for g in [X.leg1, X.leg2]:
#                            if deltaR(l.eta(), l.phi(), g.eta(), g.phi()) < 0.5:
#                                overlap = True
#                                break
#                if not overlap:
#                    goodXXs.append(XX)
            #Pair best XX to best Z
            if len(goodXXs) > 0:
                sortedXXs = sorted(goodXXs, key = lambda x: (x.p4()+bestZ.p4()).pt())[0:3]
                #                bestXX = max(goodXXs, key = lambda x: x.x1.pt()+x.x2.pt())
                bestXX = min(sortedXXs,  key = lambda x: deltaR(x.x1.leg1.eta(),x.x1.leg1.phi(),x.x1.leg2.eta(),x.x1.leg2.phi()) + deltaR(x.x2.leg1.eta(), x.x2.leg1.phi(), x.x2.leg2.eta(), x.x2.leg2.phi()))
                bestZXX = ZXX(bestZ, bestXX)
                bestZXX.otherLeptons = len(goodLeptons) - 2
                bestZXX.hasFSR = self.checkFSR_ZXX(bestZ, bestXX)
                event.ZXX.append(bestZXX)
                nZPairs+=1


        # If no Zs, search for Ws
        if len(Ws) > 0 and nZPairs == 0:
            bestW = max(Ws, key = lambda x: x.leg1.pt())
            goodXs = Xs
#            goodXs = []
#            for X in Xs:
#                for W in Ws:
#                    overlap = False
#                    for g in [X.leg1, X.leg2]:
#                        if deltaR(W.leg1.eta(), W.leg1.phi(), g.eta(), g.phi()) < 0.5:
#                            overlap = True
#                            break
#                if not overlap:
#                    goodXs.append(X)
                        
            if len(goodXs) > 0:
                #                bestX = max(goodXs, key=lambda x: x.leg1.pt() + x.leg2.pt())
                bestX = min(goodXs, key = lambda x: (x.p4()+bestW.p4()).pt())
                bestWX = Pair(bestW, bestX)
                bestWX.otherLeptons = len(goodLeptons) -1
                bestWX.deltaPhi_g1 = deltaPhi(bestW.leg1.phi(),bestX.leg1.phi())
                bestWX.deltaPhi_g2 = deltaPhi(bestW.leg1.phi(),bestX.leg2.phi())


                # Check for electrons from z's misID'd as photons
                misID = 0
                if abs(bestW.leg1.pdgId()) == 11:
                    for l in leptons:
                        if abs(l.pdgId()) != 11 or bestW.leg1.charge()+l.charge()!=0:
                            continue
                        for x in [bestX.leg1, bestX.leg2]:
                            newP4 = bestW.leg1.p4()+l.p4() + x.p4(2)
                            mass = newP4.mass()
                            if mass < 100 and mass > 80:
                                misID = 1
                        newP4 = bestW.leg1.p4() + l1.p4() + bestX.leg1.p4(2) + bestX.leg2.p4(2)
                        mass = newP4.mass()
                        if mass < 100 and mass > 80:
                            misID = 1
                    for x in [bestX.leg1, bestX.leg2]:
                        newP4 = bestW.leg1.p4() + x.p4(2)
                        mass = newP4.mass()
                        if mass < 100 and mass > 80:
                            misID = 1
                    newP4 = bestW.leg1.p4() + bestX.leg1.p4(2) + bestX.leg2.p4(2)
                    mass = newP4.mass()
                    if mass < 100 and mass > 80:
                        misID = 1
                bestWX.misID = misID
                event.WX.append(bestWX)
                if debug:
                    self.log(signal, genLeptons)
                    print "WX Vertex: ({:.2f}, {:.2f}, {:.2f}) valid={}".format(bestWX.leg2.vertex15.vertex[0],bestWX.leg2.vertex15.vertex[1],bestWX.leg2.vertex15.vertex[2], bestWX.leg2.vertex15.valid)
                    print "X Photons:\n   ",
                    print bestWX.leg2.leg1
                    print "   ",
                    print bestWX.leg2.leg2
                    print "Gen Photons"
                    for g in genPhotons:
                        print "    vertex=({:.2f}, {:.2f}, {:.2f}), pt={:.2f}, eta={:.2f}, phi={:.2f}, energy={:.2f}, mother={}".format(g.vx(),g.vy(),g.vz(),g.pt(),g.eta(),g.phi(),g.energy(),g.mother().pdgId())
                    print "Reco photons"
                    for p in photons:
                        print "    ",
                        print p
                    
                    import pdb
                    pdb.set_trace()

            goodXXs = XXs

#            goodXXs = []
#            for XX in XXs:
#                overlap = False
#                for W in Ws:
#                    for X in [XX.x1, XX.x2]:
#                        for g in [X.leg1, X.leg2]:
#                            if deltaR(W.leg1.eta(), W.leg1.phi(), g.eta(), g.phi()) < 0.5:
#                                overlap = True
#                                break
#                if not overlap:
#                    goodXXs.append(XX)

            if len(goodXXs) > 0:
                #bestXX = max(goodXXs, key = lambda x: x.x1.pt()+x.x2.pt())
                #bestWXX = WXX(bestW, bestXX)
                sortedXXs = sorted(goodXXs, key = lambda x: (x.p4()+bestW.p4()).pt())[0:3]
                bestXX = min(sortedXXs,  key = lambda x: deltaR(x.x1.leg1.eta(),x.x1.leg1.phi(),x.x1.leg2.eta(),x.x1.leg2.phi()) + deltaR(x.x2.leg1.eta(), x.x2.leg1.phi(), x.x2.leg2.eta(), x.x2.leg2.phi()))
                bestWXX = WXX(bestW, bestXX)
                bestWXX.otherLeptons = len(goodLeptons) - 1
                bestWXX.deltaPhi_X1_g1 = deltaPhi(bestW.leg1.phi(), bestXX.x1.leg1.phi())
                bestWXX.deltaPhi_X1_g2 = deltaPhi(bestW.leg1.phi(), bestXX.x1.leg2.phi())
                bestWXX.deltaPhi_X2_g1 = deltaPhi(bestW.leg1.phi(), bestXX.x2.leg1.phi())
                bestWXX.deltaPhi_X2_g2 = deltaPhi(bestW.leg1.phi(), bestXX.x2.leg2.phi())
                event.WXX.append(bestWXX)
                if debug:
                    self.log(signal, genWs)
                    print "X1 Vertex: ({:.2f}, {:.2f}, {:.2f}) valid={}".format(bestWXX.XX.x1.vertex15.vertex[0],bestWXX.XX.x1.vertex15.vertex[1],bestWXX.XX.x1.vertex15.vertex[2], bestWXX.XX.x1.vertex15.valid)
                    print "X2 Vertex: ({:.2f}, {:.2f}, {:.2f}) valid={}".format(bestWXX.XX.x2.vertex15.vertex[0],bestWXX.XX.x2.vertex15.vertex[1],bestWXX.XX.x2.vertex15.vertex[2], bestWXX.XX.x2.vertex15.valid)
                    print "X1 Photons:\n   ",
                    print bestWXX.XX.x1.leg1
                    print "   ",
                    print bestWXX.XX.x2.leg2
                    print "X2 Photons:\n   ",
                    print bestWXX.XX.x2.leg1
                    print "   ",
                    print bestWXX.XX.x1.leg2
                    print "Gen Photons"
                    for g in genPhotons:
                        print "    vertex=({:.2f}, {:.2f}, {:.2f}), pt={:.2f}, eta={:.2f}, phi={:.2f}, energy={:.2f}, mother={}".format(g.vx(),g.vy(),g.vz(),g.pt(),g.eta(),g.phi(),g.energy(),g.mother().pdgId())
                    print "Reco photons"
                    for p in photons:
                        print "    ",
                        print p
                    
                    import pdb
                    pdb.set_trace()

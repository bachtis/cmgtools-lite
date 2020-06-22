from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import *
from CMGTools.DisplacedDiPhotons.analyzers.vhGGTypes  import *
import PhysicsTools.HeppyCore.framework.config as cfg


vhGGTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='vhGGTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
#     PDFWeights = PDFWeights,
     globalVariables = [
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"),
         NTupleVariable("Flag_badChargedHadronFilter", lambda ev: ev.badChargedHadron, help="bad charged hadron filter decision"),
         NTupleVariable("Flag_badMuonFilter", lambda ev: ev.badMuon, help="bad muon filter decision"),
     ],
     globalObjects =  {
            "met" : NTupleObject("met", metType, help="PF E_{T}^{miss}, after type 1 corrections"),
     },

     collections = {
        "ZX" : NTupleCollection("ZX",ZXType ,5, help="ZX candidates"),
        "WX" : NTupleCollection("WX",WXType ,5, help="WX candidates"),
        "ZXX": NTupleCollection("ZXX", ZXXType, 5, help = "ZXX candidates"),
        "WXX": NTupleCollection("WXX", WXXType, 5, help = "WXX candidates"),
        "GenPhoton": NTupleCollection("GenPhoton", genType, 6, help = "MC photons"),
        "looseZX": NTupleCollection("looseZX", looseZXType, 5, help = "Loose ZX Candidates"),
        "looseWX": NTupleCollection("looseWX", looseWXType, 5, help = "Loose WX Candidates"),
        "looseZXX": NTupleCollection("looseZXX", looseZXXType, 5, help = "Loose ZXX Candidates"),
        "looseWXX": NTupleCollection("looseWXX", looseWXXType, 5, help = "Loose WXX Candidates"),
     }
)

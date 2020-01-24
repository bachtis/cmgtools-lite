from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import * 

import ROOT
dummyLV=ROOT.math.XYZTLorentzVector(0.0,0.0,0.0001,0.0001)

conversionType = NTupleObjectType("conversion",baseObjectTypes = [], variables = [
    NTupleVariable("convX", lambda x: x.conversionVertex().position().x(), float , help = "Reco Conversion Vertex X"),
    NTupleVariable("convY", lambda x: x.conversionVertex().position().y(), float , help = "Reco Conversion Vertex Y"),
    NTupleVariable("convZ", lambda x: x.conversionVertex().position().z(), float , help = "Reco Conversion Vertex Z"),
])


ZType = NTupleObjectType("PairType", baseObjectTypes=[fourVectorType], variables = [
    NTupleVariable("deltaPhi",   lambda x : x.deltaPhi(), float),       
    NTupleVariable("deltaR",   lambda x : x.deltaR(), float),       
    NTupleSubObject("l1",  lambda x : x.leg1,leptonType),
    NTupleSubObject("l2",  lambda x : x.leg2,leptonType),    
])

XType = NTupleObjectType("PhotonPair", baseObjectTypes=[fourVectorType], variables = [
    NTupleVariable("deltaPhi",   lambda x : x.deltaPhi(), float),       
    NTupleVariable("deltaR",   lambda x : x.deltaR(), float),       
    NTupleSubObject("g1",  lambda x : x.leg1,photonType),
    NTupleSubObject("g2",  lambda x : x.leg2,photonType),    
#    NTupleSubObject("conv1", lambda x: x.conv1, conversionType),
#    NTupleSubObject("conv2", lambda x: x.conv2, conversionType),
])

WType = NTupleObjectType("PairType", baseObjectTypes=[fourVectorType], variables = [
    NTupleVariable("deltaPhi",   lambda x : x.deltaPhi(), float),       
    NTupleVariable("deltaR",   lambda x : x.deltaR(), float),       
    NTupleSubObject("l1",  lambda x : x.leg1,leptonType),
    NTupleSubObject("l2",  lambda x : x.leg2,metType),    
])


ZXType = NTupleObjectType("ZXType", baseObjectTypes=[], variables = [
    NTupleSubObject("Z",  lambda x : x.leg1,ZType),
    NTupleSubObject("X",  lambda x : x.leg2,XType),
    NTupleVariable("otherLeptons",   lambda x : x.otherLeptons, int),          
])

XXType = NTupleObjectType("XXType", baseObjectTypes=[], variables = [
    NTupleSubObject("X1", lambda x : x.x1, XType),
    NTupleSubObject("X2", lambda x : x.x2, XType),
    NTupleVariable("m13", lambda x: x.m13, float),
    NTupleVariable("m24", lambda x: x.m24, float),
    NTupleVariable("m14", lambda x: x.m14, float),
    NTupleVariable("m23", lambda x: x.m23, float),
])

WXType = NTupleObjectType("WXType", baseObjectTypes=[], variables = [
    NTupleSubObject("W",  lambda x : x.leg1,WType),
    NTupleSubObject("X",  lambda x : x.leg2,XType),
])

ZXXType = NTupleObjectType("ZXXType", baseObjectTypes=[], variables = [
    NTupleSubObject("Z", lambda x: x.Z, ZType),
    NTupleSubObject("XX", lambda x: x.XX, XXType)
])

WXXType = NTupleObjectType("WXXType", baseObjectTypes=[], variables = [
    NTupleSubObject("W", lambda x: x.W, WType),
    NTupleSubObject("XX", lambda x: x.XX, XXType)
])

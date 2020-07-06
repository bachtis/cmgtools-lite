#!/usr/bin/env python

import ROOT
from array import array
import os, sys, re, optparse,pickle,shutil,json

def returnString(func):
    if func.GetName().find("pol")!=-1:
        st='0'
        for i in range(0,func.GetNpar()):
            st=st+"+("+str(func.GetParameter(i))+")"+("*MH"*i)
        return st    
    elif func.GetName().find("llog")!=-1:
        return str(func.GetParameter(0))+"+"+str(func.GetParameter(1))+"*log(MH)"
    if func.GetName().find("laur")!=-1:
        st='0'
        for i in range(0,func.GetNpar()):
            st=st+"+("+str(func.GetParameter(i))+")"+"/MH^"+str(i)
        return st    

    else:
        return ""

parser = optparse.OptionParser()
parser.add_option("-g","--graphs",dest="graphs",default='',help="Comma   separated graphs and functions to fit  like MEAN:pol3,SIGMA:pol2")
parser.add_option("-o","--output",dest="output",help="Output JSON",default='')
parser.add_option("-d","--debugFile",dest="debugFile",help="Output debug file",default='')
parser.add_option("-m","--min",dest="min",type=float, help="minimum x",default=0)
parser.add_option("-M","--max",dest="max",type=float, help="maximum x",default=0)


(options,args) = parser.parse_args()
#define output dictionary


rootFile=ROOT.TFile(args[0])


graphStr= options.graphs.split(',')
parameterization={}



ff=ROOT.TFile(options.debugFile,"RECREATE")
ff.cd()
for string in graphStr:
    comps =string.split(':')      
    graph=rootFile.Get(comps[0])
    if comps[1].find("pol")!=-1:
        func=ROOT.TF1(comps[1],comps[1],0,13000)
    elif  comps[1]=="llog":
        func=ROOT.TF1("llog","[0]+[1]*log(x)",1,13000)
        func.SetParameters(1,1)
    elif  comps[1].find("laur")!=-1:
        order=int(comps[1].split("laur")[1])
        st='0'
        for i in range(0,order):
            st=st+"+["+str(i)+"]"+"/x^"+str(i)
        print 'Laurent String',st    
        func=ROOT.TF1(comps[1],st,1,13000)
        for i in range(0,order):
            func.SetParameter(i,0)

    ## hardcoded fix for pathological slopes:
    if comps[0]=="slope":
        N=graph.GetN()
        Y=graph.GetY()
        for i in range(N):
            if Y[i]>0 or Y[i]<-0.08:
                graph.RemovePoint(i)
    
    graph.Fit(func,"","",options.min,options.max)
    parameterization[comps[0]]=returnString(func)
    graph.Write(comps[0])
    func.Write(comps[0]+"_func")

ff.Close()
f=open(options.output,"w")
json.dump(parameterization,f)
f.close()


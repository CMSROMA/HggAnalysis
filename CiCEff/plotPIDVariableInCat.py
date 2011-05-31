#! /usr/bin/env python

from ROOT import *

import sys
import os
import math


                               
gROOT.SetStyle("Plain")
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
#gROOT.SetBatch()

variables=["isosumoet","isosumoetbad","trkisooet","sieie","hovere","r9","drtotk_25_99","pixel"]
categories=["0","1","2","3"]
#f=TFile("WJetsMADGRAPH_Wenu_preselected.root");
f=TFile("../redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6_00.root")

histogram = {} 
for variable in variables:
    for cat in categories:
####EoPClass_hadrons_subdet2_ptbin1_class1
        print variable+"_cat"+cat+"_"
        histo=variable+"_cat"+cat+"_"
        histogram[histo]=f.Get(histo)
        print histogram[histo].Integral()

canvas = TCanvas("c","c",1)


for variable in variables:
    print variable
    maximum=-9999.;
    maximumX=-9999.;
    minimumX=9999.;
    leg=TLegend(0.8,0.8,0.99,0.99)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    for cat in categories:
        histoName=variable+"_cat"+cat+"_"
        #put underflow and overflow in histoRange and renormalize
        histogram[histoName].SetBinContent(histogram[histoName].GetNbinsX(),histogram[histoName].GetBinContent(histogram[histoName].GetNbinsX())+histogram[histoName].GetBinContent(histogram[histoName].GetNbinsX()+1))
        histogram[histoName].SetBinContent(1,histogram[histoName].GetBinContent(0)+histogram[histoName].GetBinContent(1))
        histogram[histoName].Scale(1./histogram[histoName].Integral())
        if (histogram[histoName].GetMaximum()>=maximum):
            maximum = histogram[histoName].GetMaximum()
        if (histogram[histoName].GetXaxis().GetXmax()>=maximumX):
            maximumX = histogram[histoName].GetXaxis().GetXmax()
        if (histogram[histoName].GetXaxis().GetXmin()<=minimumX):
            minimumX = histogram[histoName].GetXaxis().GetXmin()
    print maximum,maximumX,minimumX
    maximum=maximum*1.8
    if (variable=="isosumoet" or
        variable == "isosumoetbad" or
        variable == "trkisooet" or
        variable ==  "hovere" or
        variable == "drtotk_25_99"):
        maximum=maximum*30.
    if (variable=="OneOverEMinusOneOverP"):
        minimumX=-0.04
        maximumX=0.04
    a=TH2F("a","a",10,minimumX,maximumX,10,0.0000001,maximum)
    a.GetXaxis().SetTitle(variable)
    a.Draw()
    i=0
    for cat in categories:
        if (i<4):
            icolor = i+1
        elif (icolor <9):
            icolor = i+2
        else:
            icolor = i+3
        histoName=variable+"_cat"+cat+"_"
        histogram[histoName].SetMarkerColor(icolor)
        histogram[histoName].SetMarkerSize(1.2)
        histogram[histoName].SetLineColor(icolor)
        #                    histogram[histoName].SetLineStyle(icolor)
        histogram[histoName].SetLineWidth(3)
        leg.AddEntry(histogram[histoName],"cat"+cat,"l")
        #                    histogram[histoName].SetMaximum(maximum)
        #                    histogram[histoName].SetMinimum(0.0000001)
        #                    histogram[histoName].GetXaxis().SetRangeUser(minimumX,maximumX)
        histogram[histoName].Draw("HSAME")
        i=i+1
        if (variable=="isosumoet" or
            variable == "isosumoetbad" or
            variable == "trkisooet" or
            variable ==  "hovere" or
            variable == "drtotk_25_99"):
            c.SetLogy(1)
        else:
            c.SetLogy(0)

    leg.Draw()
    c.SaveAs(variable+".png")
                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                        


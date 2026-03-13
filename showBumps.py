# Show most significant bumps 
from globalAux import *
import numpy as np
import pyBumpHunter as BH
from math import log
import math,os,json
import ROOT
from ROOT import TCanvas, TPostScript, TFile, TLegend, gPad, TF1, TRandom3, TH1D, TMath
from array import array

myinput="interactive"
# trigger type
if (len(sys.argv) ==2):
   myinput =sys.argv[1]

# -----------------------------
# Configuration
# -----------------------------
CHANNELS = ["jj", "jb", "bb", "je", "jm", "jg", "be", "bm", "bg"]

def trigger_settings(trig_type: int) -> tuple[int, str]:
    """Return (xmin, label) based on trigger type."""
    settings = {
        1: (600, r"T1:\; MET"),
        2: (300, r"T2:\; 1 \ell"),
        3: (300, r"T3:\; 2 \ell"),
        4: (500, r"T4:\; 1 γ"),
        5: (300, r"T5:\; 2 γ"),
        6: (1500, r"T6:\; 1 jet"),
        7: (600, r"T7:\; 4 jets"),
    }
    return settings.get(trig_type, (300, "1 lep"))


# some default  min and max values X, Y ranges
XMAX = 9000
XMIN = 300
YMIN = 0.81
YMAX = 100000

# write into file
rootfile="figs/bumps.root"
hfile=TFile(rootfile)
hfile.ls()



hname="histo_4_jg_bump0"
bnam="back_4_jg_bump0"

chanks=hname.split("_")
TRIG_TYPE=int(chanks[1])

XMIN, TLABEL = trigger_settings(TRIG_TYPE)

hhD=hfile.Get(hname)
hhBak=hfile.Get(bnam)


figdir="figs/"
name=os.path.basename(__file__)
name=name.replace(".py","")
epsfig=figdir+name+".eps"
ps1 = TPostScript( epsfig,113)

# print("Integral=",hhD.Integrate())
c1=TCanvas("c_massjj","BPRE",10,10,500,500);
c1.cd(1);
gPad.SetLogy(1)
gPad.SetLogx(1)
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.12)
gPad.SetLeftMargin(0.14)
gPad.SetRightMargin(0.04)
hhD.SetAxisRange(YMIN, YMAX,"y");
hhD.SetAxisRange(XMIN, XMAX,"x");
hhD.SetTitle("")
hhD.SetMarkerColor(1)
hhD.SetMarkerSize(0.8)
hhD.SetMarkerStyle(20)
hhD.SetStats(0)
hhD.Draw("pe")
hhBak.SetLineColor(2)
hhBak.Draw("l same")
hhD.GetXaxis().SetTitle( "Mass [GeV]" );
hhD.GetYaxis().SetTitle( "Pseudo Events" );

leg2=TLegend(0.5, 0.55, 0.89, 0.84);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.04);
leg2.SetHeader("BUMP")
leg2.AddEntry(hhBak,"Fit template","lp")
leg2.AddEntry(hhD,"Obtained sample","lp")
leg2.Draw("same");
hhD.Draw("pe same")

print (epsfig)
ps1.Close()

c1.Update()
if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);


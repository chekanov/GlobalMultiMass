from math import log
import math,os,sys,json
import ROOT
from ROOT import TCanvas, TPostScript, TLegend, gPad, TF1, TRandom3, TH1D, TMath
from array import array

myinput="interactive"
# trigger type
if (len(sys.argv) ==2):
   myinput =sys.argv[1]

# -----------------------------
# Configuration
# -----------------------------
CHANNELS = ["jj", "jb", "bb", "je", "jm", "jg", "be", "bm", "bg"]

# we assume that some fraction of jj events ends up in as jb, bb, etc (in decreasing oder)
# For example, "jb" = 0.4 means that 40% of events from "jj" are in "jb".  
# Any other possiblity would require samplings from SM Monte Carlo, which is impossible to do
# since we need billions of events from Monte Carlo. This also assumes that overlaping shape repeats
# the shape of jj (which is a good assamption since all such masses are basically p5 with very similar shapes
# For example, jb:0.4 means that 40% of events originate from jj.
#              je:0.2 means that 20% of events originte from jj etc 
# Trigger-dependent overlap values 
# Should be modified using Wasikul's plots
DEFALT_OVERLAP1={"jj":0.0, "bb":0.41, "jb":0.63, "jm":0.34, "je":0.34, "jg":0.06, "be":0.20, "bm":0.2, "bg":0.01}  
DEFALT_OVERLAP2={"jj":0.0, "bb":0.27, "jb":0.37, "jm":0.53, "je":0.53, "jg":0.01, "be":0.24, "bm":0.25, "bg":0.01}
DEFALT_OVERLAP3={"jj":0.0, "bb":0.01, "jb":0.28, "jm":0.53, "je":0.52, "jg":0.01, "be":0.16, "bm":0.19, "bg":0.01}
DEFALT_OVERLAP4={"jj":0.0, "bb":0.13, "jb":0.14, "jm":0.02, "je":0.02, "jg":0.99, "be":0.01, "bm":0.01, "bg":0.24}
DEFALT_OVERLAP5={"jj":0.0, "bb":0.03, "jb":0.16, "jm":0.02, "je":0.01, "jg":0.99, "be":0.01, "bm":0.01, "bg":0.16}
DEFALT_OVERLAP6={"jj":0.0, "bb":0.31, "jb":0.47, "jm":0.04, "je":0.05, "jg":0.12, "be":0.02, "bm":0.02, "bg":0.02}
DEFALT_OVERLAP7={"jj":0.0, "bb":0.43, "jb":0.57, "jm":0.04, "je":0.05, "jg":0.06, "be":0.03, "bm":0.03, "bg":0.02}
# This is overlap for different triggers
DEFALT_OVERLAP_TRIGGER={1:DEFALT_OVERLAP1,2:DEFALT_OVERLAP2, 3:DEFALT_OVERLAP3, 
                        4:DEFALT_OVERLAP4, 5:DEFALT_OVERLAP5, 6:DEFALT_OVERLAP6,7:DEFALT_OVERLAP7} 


# Maximum pseudo-experiments for global p-value
# for 5 sigmal local, set to maximum value, like 10^6!
MaxEvents=10


## Expected local significance as in BumpHunter
ExpectedLocalZvalue=4  



def p_to_z_value (p, excess) :
  """the function normal_quantile converts a p-value into a significance,
  i.e. the number of standard deviations corresponding to the right-tail of 
  a Gaussian"""
  if excess :
    zval = ROOT.Math.normal_quantile(1-p,1);
  else :
    zval = ROOT.Math.normal_quantile(p,1);
  return zval


def z_to_p_value(z_value):
    """
    Converts a Z-score to a one-sided p-value using ROOT::Math.
    """
    # normal_cdf_c(z, sigma=1, x0=0) computes 1 - Phi(z) directly
    p_value = ROOT.Math.normal_cdf_c(z_value)
    return p_value


#“Asimov” (profile likelihood) significance (recommended)
# when S<<D
# Often used in HEP; behaves better, especially when SS is not ≪D≪D:
def asimov_significance(S, D):
    """
    Asimov (profile likelihood) significance for counting experiment:
        Z = sqrt( 2 * [ (S + D) * ln(1 + S/D) - S ] )

    S: expected signal (>= 0)
    D: expected background (> 0)
    returns: Z in "sigma" units (float)
    """
    S = float(S)
    D = float(D)

    if S < 0:
        raise ValueError("S must be >= 0")
    if D <= 0:
        raise ValueError("D must be > 0")

    z2 = 2.0 * ((S + D) * math.log1p(S / D) - S)  # log1p(x) = ln(1+x)
    if z2 < 0:  # guard against tiny negative due to rounding
        z2 = 0.0
    return math.sqrt(z2)



# get maximum X of a histogram 
def getMaxHistoX(h1, xmin):
    xaxis = h1.GetXaxis()
    Ntot = xaxis.GetNbins()+1
    xmax=0
    for i in range(Ntot):
        y1 = h1.GetBinContent(i+1)
        x1 = h1.GetBinCenter(i+1)
        if (y1 < xmin):
             xmax= x1
             break
    return xmax 




print("Searching for bumps with Z=",ExpectedLocalZvalue," which is ",z_to_p_value(ExpectedLocalZvalue)," p-value")



# Make Bin=30 to fluctuate by 500 events for debugging!
# Comment this out to remove this fluctuation
FluctuateBin={}
# FluctuateBin={30: 500}

# default trigger. In future we will add other independent triggers
TRIG_TYPE = 2

# default histogram Bins
mjjBinsL = [99,112,125,138,151,164,177,190, 203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156];
mjjBins = array("d", mjjBinsL)


# Get random numbers
r=TRandom3()


# some X , and Y ranges 
XMAX = 9000
XMIN = 300
YMIN = 0.81
YMAX = 100000 

# counter for significant events
NrFound=0


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


XMIN, TLABEL = trigger_settings(TRIG_TYPE)


# -----------------------------
# ROOT fit background templates extraced from 1% of data 
# -----------------------------
class FiveParam2015:
    """Standard P5 background function."""

    def __call__(self, x, par):
        ecm = 13000.0
        xx = x[0] / ecm

        ff1 = par[0] * TMath.Power((1.0 - xx), par[1])
        ff2 = TMath.Power(xx, (par[2] + par[3] * log(xx) + par[4] * log(xx) * log(xx)))
        return ff1 * ff2


class FiveParam2015Gauss:
    """P5 background + Gaussian peak."""

    def __call__(self, x, par):
        ecm = 13000.0
        xx = x[0] / ecm

        ff1 = par[0] * TMath.Power((1.0 - xx), par[1])
        ff2 = TMath.Power(xx, (par[2] + par[3] * log(xx) + par[4] * log(xx) * log(xx)))
        background = ff1 * ff2

        sigma = par[7]
        gauss = par[5] * TMath.Gaus(xx, par[6], sigma) if sigma > 0 else 0.0
        return background + gauss

# -----------------------------
# Main loop to create backgrond function and make random histograms 
# -----------------------------
# map with collection of histograms
Histograms={}
Histograms_OriginalEvents={}
Histograms_FromOverlap={}
BackgroudFunction={}
Histograms_fromJJ={}

# use it for debug
# CHANNELS=["jj"]


# Just count histograms..
Ntot=0

## fill from the files parameters
print("Read all paramters from JSON") 
mypar={}
for TRIG_TYPE in range(1, 8):
     for channel in CHANNELS:
        fitfile = f"fits/fitme_p5_t{TRIG_TYPE}_{channel}.json"
        if os.path.isfile(fitfile) is False:
                continue
        with open(fitfile, "r") as jfile:
               data = json.load(jfile)
               mypar[fitfile] = data

print("Loop over events",MaxEvents)
for event in range(MaxEvents):

    BumpFound=False # so far no bump found for this run


    if (event %1000 == 0 ): print("Event=",event)

    # loop over 7 triggers
    for TRIG_TYPE in range(1, 8):
        # print("TRIGGER=", TRIG_TYPE)
        DEFALT_OVERLAP = DEFALT_OVERLAP_TRIGGER[TRIG_TYPE] # get overlap for this trigger

        # Loop over each mass for a channel
        for channel in CHANNELS:

            # get saved paramters
            fitfile = f"fits/fitme_p5_t{TRIG_TYPE}_{channel}.json"
            if fitfile not in mypar: continue
            data=mypar[fitfile]
            fit_min, fit_max = XMIN, XMAX
            # Background-only TF1 (5 params)
            name = f"{TRIG_TYPE}_{channel}"
            mback=FiveParam2015()
            back = TF1("back_" + name, mback, fit_min, fit_max, 5)

            # Background + signal TF1 (8 params)
            mbacksig=FiveParam2015Gauss()
            backsig = TF1(f"sig_{name}", mbacksig, fit_min, fit_max, 8)

            parameters = data["parameters"]
            nom_func = data["name"]
            errors = data["errors"]
            ndf = int(data["ndf"])
            chi2 = float(data["chi2"])
            fit_min = float(data["fmin"])
            fit_max = float(data["fmax"])

            """
            chi2_ndf = chi2 / ndf if ndf else float("inf")
              print(
                "Nominal Fit chi2/ndf=",
                chi2_ndf,
                " Fit parameters=",
                parameters,
                " errors=",
                errors,
            )
            """

            for i, value in enumerate(parameters):
                value = float(value)
                #print(f"p{i}={value}")
                back.SetParameter(i, value)
                if value == 0.0:
                    back.FixParameter(i, value)

            # we only do this once for JJ, since overlap for other channels
            # will be obtained from this one. We fluctuate bins according to Poisson
            # the fit function
            if channel == "jj":
                hbackJJ_name = f"histoJJ_{TRIG_TYPE}_{channel}"
                hbackJJ = TH1D(hbackJJ_name, hbackJJ_name, len(mjjBins) - 1, mjjBins)
                hbackJJ.SetTitle(hbackJJ_name)
                hbackJJ.SetName(hbackJJ_name)
                hbackJJ.SetDirectory(0)
                #print("Created JJ template")

                residuals=[]
                for i in range(hbackJJ.GetNbinsX() - 1):
                    center = hbackJJ.GetBinCenter(i + 1)
                    B = back.Eval(center)
                    pseudo = r.PoissonD(B)

                    # only for debugging.. Make it fluctuate 
                    if (FluctuateBin != None): 
                                    if (i in FluctuateBin): pseudo = pseudo+FluctuateBin[i] # just for debuggin.. Make this bin to fluctuate! 


                    residuals.append( (pseudo - B) / B) # we keep relative deviation due to Poisson statistics 
                    if pseudo > 0:
                        hbackJJ.SetBinContent(i + 1, pseudo)
                        hbackJJ.SetBinError(i + 1, TMath.Sqrt(pseudo))
                    else:
                        hbackJJ.SetBinContent(i + 1, 0)
                        hbackJJ.SetBinError(i + 1, 0)


            # create overlap histogram as in JJ data
            # we use JJ as a template, but scale according to overlaps..
            # this way we keep info of exact Poisson fluctuations
            hback1_name = f"histo1_{TRIG_TYPE}_{channel}"
            hback1 = hbackJJ.Clone()
            hback1.SetTitle(hback1_name)
            hback1.SetName(hback1_name)
            hback1.SetDirectory(0)

            # Now we create 2nd histogram using the background function for this channel 
            # do not fluctauate it as we do this later
            hback2_name = f"histo2_{TRIG_TYPE}_{channel}"
            hback2 = hbackJJ.Clone() # we clone it from JJ 
            hback2.SetTitle(hback2_name)
            hback2.SetName(hback2_name)
            hback2.SetDirectory(0)
            for i in range(hback2.GetNbinsX() - 1):
                center = hback2.GetBinCenter(i + 1)
                B = back.Eval(center)  # get this value from function
                hback2.SetBinContent(i + 1, B)

            # we will keep track of the previous histogram
            hback3_name = f"histo3_{TRIG_TYPE}_{channel}"
            hback3 = hback2.Clone() # we clone it from JJ 
            hback3.SetTitle(hback1_name)
            hback3.SetName(hback1_name)
            hback3.SetDirectory(0)
 
            # now make sure fraction comes from JJ
            #F1=hback1.Integral()
            #F2=hback2.Integral()
            # No our JJ histogram has same normalization as hback2 * DEFALT_OVERLAP[channel]
            #hback1.Scale(F2/F1)
            #hback1.Scale(DEFALT_OVERLAP[channel])
            # now we get fluctuations, and scale them..
            # they represent deviations from backround after scaling them down to the expected
            # fraction..
            #newresiduals=[]
            #for res in residuals:
            #     newresiduals.append( res* DEFALT_OVERLAP[channel] * (F2/F1))


            # now we refill overlap histogram using DEFALT_OVERLAP[channel] asamption
            # but we will keep same shape. Use residuals from JJ to modify fluctuations
            for i in range(hback2.GetNbinsX() - 1):
                center=hback2.GetBinCenter(i + 1)
                events2=hback2.GetBinContent(i + 1)
                expected = (events2 * DEFALT_OVERLAP[channel])

                ## keep track expected without fluctuations 
                hback3.SetBinContent(i + 1, expected)
                hback3.SetBinError(i + 1, TMath.Sqrt(expected))


                events=expected + expected * residuals[i] # fluctuate according to relative fluctuations of JJ 
                hback1.SetBinContent(i + 1, events)
                hback1.SetBinError(i + 1, TMath.Sqrt(events))

                # keep the main histogram random
                B = back.Eval(center)
                pseudo = r.PoissonD(B)
                # also keep random 
                if pseudo > 0:
                        hback2.SetBinContent(i + 1, pseudo)
                        hback2.SetBinError(i + 1, TMath.Sqrt(pseudo))
                else:
                        hback2.SetBinContent(i + 1, 0)
                        hback2.SetBinError(i + 1, 0)
  
                


            # and build remaining part of the histogram subtracting hback1 from hback2
            # This subtracted entries are truly independent JJ, so we fluctuate them..
            for i in range(hback2.GetNbinsX() - 1):
                center = hback2.GetBinCenter(i + 1)
                RemainingContent = hback2.GetBinContent(i + 1) - hback3.GetBinContent(i + 1)
                if RemainingContent < 0:
                    RemainingContent = 0

                # fluctuate according to Poisson statistics
                pseudo = r.PoissonD(RemainingContent)
                if pseudo > 0:
                    hback2.SetBinContent(i + 1, pseudo)
                    hback2.SetBinError(i + 1, TMath.Sqrt(pseudo))
                else:
                    hback2.SetBinContent(i + 1, 0)
                    hback2.SetBinError(i + 1, 0)

            # Finally, we can combine 2 histograms
            # This histogram would have 2 parts, with overlap and correlated fluctuations
            hback_name = f"histo_{TRIG_TYPE}_{channel}"
            hback = hback2.Clone()
            hback.SetTitle(hback_name)
            hback.SetName(hback_name)
            hback.SetDirectory(0)
            hback.Add(hback1)

            # finally, clean-up all these things..
            # adding 2 histogram where 1 is a template may have content <1. Fix
            for i in range(hback.GetNbinsX() - 1):
                center = hback.GetBinCenter(i + 1)
                y=hback.GetBinContent(i + 1)
                if (y<1.0):
                    hback.SetBinContent(i + 1, 0)
                    hback.SetBinError(i + 1, 0)
                y=hback1.GetBinContent(i + 1)
                if (y<1.0):
                    hback1.SetBinContent(i + 1, 0)
                    hback1.SetBinError(i + 1, 0)



            # collect histograms
            BackgroudFunction[hback_name] = back  
            Histograms[hback_name] = hback
            Histograms_fromJJ[hback_name] = hbackJJ 
            Histograms_FromOverlap[hback1_name] = hback1     # events coming from overlap
            Histograms_OriginalEvents[hback2_name] = hback2  # original events

            # Now, search for 3-sigma local bump in the histogram "hback", and count the event.
            # How many such bumps are found for the maximum number of events ?
            # The deviations are searched agains the backgrond function "back"
            ########################## Use Bump Hunter algorithm #########
            ## Run over difference data-fit, find maximum, and them loook at left and right devitations
            ## untill you see signal above > ExpectedLocalSignificance  


            # This is for debug only! 
            # if (channel != "jb" or TRIG_TYPE != 2): continue

            Ntot=Ntot+1

            ## Let's mimic BumpHunt for now. Find residuals from the fit and just look at significance of a single bin
            ## Note: Real BumpHunt also adds adjusted bins
            sign=0;
            XmaxVal=getMaxHistoX(hback,1)
            #print("XmaxVal=",XmaxVal)
            for i in range(hback.GetNbinsX() - 1):
                center = hback.GetBinCenter(i + 1)
                if (center<XMIN or center>XmaxVal): continue
                ydata = hback.GetBinContent(i + 1)
                yfit= back.Eval(center) 
                deviation=ydata-yfit
                # also need to find left and right bins! Do it for 1 bin now! 
                significance = 0
                if (deviation>0): significance=asimov_significance(deviation, yfit) 
                # print(center, significance)
                if (significance>sign): 
                                sign=significance;
                                sign_center=center

            if sign>ExpectedLocalZvalue:
                    BumpFound=True 
                    print("Found bump with significance=",sign," and postion=",sign_center)
    if (BumpFound): NrFound= NrFound+1

# the probability that background fluctuations alone (the null hypothesis) could produce a 
# result as extreme as, or more extreme than, the observed experimental data
print()
print("Total events=", MaxEvents)
print("Found events=", NrFound)
pvalue=float(NrFound)/MaxEvents
print("Global p-value=",pvalue)


figdir="figs/"
name=os.path.basename(__file__)
name=name.replace(".py","")
epsfig=figdir+name+".eps"
ps1 = TPostScript( epsfig,113)


# Plot one histogram for debugging
TRIG_TYPE=2
channel="jb"
Xmin=XMIN
Xmax=XMAX
Ymin=YMIN
Ymax=YMAX
hhD=Histograms[f"histo_{TRIG_TYPE}_{channel}"]
hhBak=BackgroudFunction[f"histo_{TRIG_TYPE}_{channel}"]
# parts of this histograms
hhD1=Histograms_OriginalEvents[ f"histo2_{TRIG_TYPE}_{channel}" ]
hhD2=Histograms_FromOverlap[ f"histo1_{TRIG_TYPE}_{channel}" ]


# plot original JJ
channel2="jj"
hhJJ=Histograms_fromJJ[f"histo_{TRIG_TYPE}_{channel2}"]
hhJJ.SetAxisRange(Ymin, Ymax,"y");
hhJJ.SetAxisRange(Xmin, Xmax,"x");
hhJJ.SetMarkerColor(4)
hhJJ.SetMarkerSize(0.5)
hhJJ.SetMarkerStyle(21)


# print("Integral=",hhD.Integrate())
c1=TCanvas("c_massjj","BPRE",10,10,500,500);
c1.cd(1);
gPad.SetLogy(1)
gPad.SetLogx(1)
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.12)
gPad.SetLeftMargin(0.14)
gPad.SetRightMargin(0.04)
hhD.SetAxisRange(Ymin, Ymax,"y");
hhD.SetAxisRange(Xmin, Xmax,"x");
hhD.SetTitle("")
hhD.SetMarkerColor(1)
hhD.SetMarkerSize(0.8)
hhD.SetMarkerStyle(20)
hhD.SetStats(0)
hhD.Draw("pe")
hhBak.SetLineColor(1)
hhBak.Draw("l same")

hhD.GetXaxis().SetTitle( "Mass [GeV]" );
hhD.GetYaxis().SetTitle( "Pseudo Events" );



hhD1.SetAxisRange(Xmin, Xmax,"x");
hhD1.SetMarkerColor(2)
hhD1.SetMarkerSize(0.8)
hhD1.SetMarkerStyle(20)
hhD1.Draw("pe same")

hhD2.SetAxisRange(Xmin, Xmax,"x");
hhD2.SetMarkerColor(3)
hhD2.SetMarkerSize(0.8)
hhD2.SetMarkerStyle(20)
hhD2.Draw("pe same")

# original JJ
hhJJ.Draw("pe same")


leg2=TLegend(0.5, 0.55, 0.89, 0.84);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.04);
leg2.SetHeader("T="+str(TRIG_TYPE)+" "+str(channel))
leg2.AddEntry(hhBak,"Fit template","lp")
leg2.AddEntry(hhD,"Obtained sample","lp")
leg2.AddEntry(hhD2,"From JJ overlap","lp")
leg2.AddEntry(hhD1,"Expected from "+str(channel),"lp")
leg2.AddEntry(hhJJ,"JJ distribution","lp")
leg2.Draw("same");
hhD.Draw("pe same")


# only for debugging.. Make it fluctuate 
if (len(FluctuateBin)>0):
            print("We fluctuated events in overlap! Check FluctuateBin")

print (epsfig)
ps1.Close()

c1.Update()
if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);

# This code runs multiple pseudo-experiments for 63 jet+X masses.
# You may adjust pyBumpHunter since it uses this package..

############### USER SETTING ######################


## Expected local significance as in BumpHunter
ExpectedLocalZvalue=6 


# Maximum pseudo-experiments for global p-value
# for 6 sigma Z (local), set to maximum value to  a large number
MaxEvents=20


# Do not process histograms with less than 50 entries 
MinEntries = 50

# if you do not want overlap, set to True:
noOverlap=False 

# Make Bin=30 to fluctuate by 500 events for debugging!
# Comment this out to remove this fluctuation
FluctuateBin={}
# FluctuateBin={30: 500}

# CM energy for fit functions
cms=13000.


# THIS IS FOR PLOTTING and debugging.. 
# NOT REFLECTED ON P-Value 
PLOT_TRIGGER_TYPE=2 
PLOT_CHANNEL="jb"


############# END USER SETTING ######################


########### DO NOT CHANGE ANYTHING BELOW ############

import sys
sys.path.append("./pyBumpHunter/")
# Old
# sys.path.append("../pyBumpHunter/")

# also import some functions...
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


if (noOverlap==True):
  print("No overlap requested")
  for d in DEFALT_OVERLAP_TRIGGER.values():
    for k in d:
        d[k] = 0
  print(DEFALT_OVERLAP_TRIGGER)

if (noOverlap==False):
          print("Running with overlaps")
else:
          print("Running without overlaps")  

print("The number of pseudo-experiments=",MaxEvents)
print("Searching for bumps with Z=",ExpectedLocalZvalue," which is ",z_to_p_value(ExpectedLocalZvalue)," p-value")
print("Min number of entries=",MinEntries)

# some default  min and max values X, Y ranges
XMAX = 9000
XMIN = 300
YMIN = 0.81
YMAX = 100000


# default histogram Bins
mjjBinsL = [99,112,125,138,151,164,177,190, 203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156];
mjjBins = array("d", mjjBinsL)


# Get random numbers
r=TRandom3()


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

# collect bumps for further analysois
BumpCollector=[] 

# Just count histograms..
Ntot=0

## fill from the files parameters
print("Read all paramters from JSON") 
mypar={}
mypar_alt={}
for TRIG_TYPE in range(1, 8):
     for channel in CHANNELS:
        fitfile = f"fits/fitme_p5_t{TRIG_TYPE}_{channel}.json"
        if os.path.isfile(fitfile) is False:
                continue
        with open(fitfile, "r") as jfile:
               data = json.load(jfile)
               mypar[fitfile] = data
        fitfile_alt = f"fits/fitme_p5alt_t{TRIG_TYPE}_{channel}.json"
        if os.path.isfile(fitfile_alt) is False:
                continue
        with open(fitfile_alt, "r") as jfile:
               data_alt = json.load(jfile)
               mypar_alt[fitfile] = data_alt 

print("Loop over events",MaxEvents)
for event in range(MaxEvents):

    BumpFound=False # so far no bump found for this run


    if (event %1000 == 0 ): print("Event=",event)

    # loop over 7 triggers
    for TRIG_TYPE in range(1, 8):

        XMIN, TLABEL = trigger_settings(TRIG_TYPE)

        # print("TRIGGER=", TRIG_TYPE)
        DEFALT_OVERLAP = DEFALT_OVERLAP_TRIGGER[TRIG_TYPE] # get overlap for this trigger

        # Loop over each mass for a channel
        for channel in CHANNELS:

            # get saved paramters
            fitfile = f"fits/fitme_p5_t{TRIG_TYPE}_{channel}.json"
            if fitfile not in mypar: continue
            data_fit=mypar[fitfile]
            fit_min, fit_max = XMIN, XMAX
            # Background-only TF1 (5 params)
            name = f"{TRIG_TYPE}_{channel}"
            mback=FiveParam2015()
            back = TF1("back_" + name, mback, fit_min, fit_max, 5)

            # Background + signal TF1 (8 params)
            mbacksig=FiveParam2015Gauss()
            backsig = TF1(f"sig_{name}", mbacksig, fit_min, fit_max, 8)

            parameters = data_fit["parameters"]
            nom_func = data_fit["name"]
            errors = data_fit["errors"]
            ndf = int(data_fit["ndf"])
            chi2 = float(data_fit["chi2"])
            fit_min = float(data_fit["fmin"])
            fit_max = float(data_fit["fmax"])

            # alternative function for systematics  
            data_fit_alt=mypar_alt[fitfile]
            parameters_alt = data_fit_alt["parameters"]

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

            TotalEvents=hback.Integral()
            if (TotalEvents<MinEntries):
                               continue # ignore low multiplicities 

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

            sign=0;
            sign_center=0
            #XmaxVal=getMaxNonzero(hback,XMIN,2.0)
            #print("XMIN=",XMIN, " XMAX=",XmaxVal) 
            # get NuPy arrays

            ## make smaller by 2*sigma
            fit_min_x=fit_min+fit_min*0.15
            fit_max_x=fit_max-fit_max*0.15
            data_x_center, data_bin_width, data_y, hist_x = get_input(hback,fit_min=fit_min,fit_max=fit_max)
            bkg_x_center = data_x_center
            p1=parameters[0]
            p2=parameters[1]
            p3=parameters[2]
            p4=parameters[3]
            p5=parameters[4]
            integral_data = None
            #bkg_function_nom = FiveParam(cms, bkg_x_center, p1, p2, p3, p4, p5)
            #bkg_sample_nom, weight_nom = construct_bkg_sample(bkg_function_nom, bkg_x_center, integral_data)
            #bkg_sample_nom=bkg_function_nom

            bkg_sample_nom=tf1_to_numpy(bkg_x_center, back)
           
            # Do not use it for now..
            #p1_alt=parameters_alt[0]
            #p2_alt=parameters_alt[1]
            #p3_alt=parameters_alt[2]
            #p4_alt=parameters_alt[3]
            #p5_alt=parameters_alt[4]
            #bkg_function_alt =FiveParam_alt(cms, bkg_x_center, p1_alt, p2_alt, p3_alt, p4_alt, p5_alt)
            #bkg_sample_alt, weight_alt = construct_bkg_sample(bkg_function_alt, bkg_x_center, integral_data)

            
            ## Construct BumpHunter and weights and alternatiev fuctions (local for Jet+X) 
            # hunter = BH.BumpHunter1D( rang=None, width_min=2, width_max=None, width_step=1, scan_step=1, npe=100, seed=666, bins=hist_x, weights=weight_nom, weights_alt=weight_alt)
            # Default from GIT. Just se 1 psedo-experiment since intrested in local significance 
            hunter = BH.BumpHunter1D( rang=None, width_min=2, width_max=30, width_step=1, scan_step=1, npe=0, seed=666, bins=hist_x)
            ## Bump Scan with systematics from Jet+X analysis 
            ##hunter.bump_scan( data = data_y, bkg = bkg_sample_nom, bkg_alt = bkg_sample_alt, do_pseudo = True, stat_only = True)
            hunter.bump_scan( data = data_y, bkg = bkg_sample_nom, is_hist=True, do_pseudo = False)

            #print("Range=",fit_min,fit_max)
            #hback.Print("All")
            #print(data_y)
            #print(bkg_function_nom)

            # local_pvalue = hunter.min_Pval_ar  # Array of minimum p-values from pseudo-experiments
            # The observed local p-value for the most significant bump
            #observed_local_pvalue = hunter.min_Pval_data
            #print(observed_local_pvalue)
            # print([attr for attr in dir(hunter) if not attr.startswith('_')])
            # Local p-value (per scan window)
            # print(hunter.signal_eval)  # Array of local p-values for each scan window

            # Print all results
            # hunter.print_bump_info()
            sign_center=0;
            # local significance ? 
            local_p = hunter.min_Pval_ar[0]
            Zval=0
            if (local_p>0 and local_p<1):
                         Zval=ROOT.RooStats.PValueToSignificance(local_p)
                         #print(hunter.bump_info(data_y))

            # also use p_back = ROOT.RooStats.SignificanceToPValue(Z)
            #print("Local significance=",local_sign)
            #print("bumploc =", getattr(hunter, "bumploc", None))
            #print("bumpwidth =", getattr(hunter, "bumpwidth", None))
            #print("Local p-value:", hunter.min_Pval)
            #print(hunter.__dict__.keys())
            # Inspect all non-private attributes to find p-value fields
            #results_attrs = {attr: getattr(hunter, attr) for attr in dir(hunter)
            #       if not attr.startswith('_') and not callable(getattr(hunter, attr))}
            #for k, v in results_attrs.items():
            #       print(f"{k}: {v}")
            # Print bump
            #state=hunter.save_state()
            #print(state["min_Pval_ar"])
            #local_p_value = hunter.p_val
            ## where local significance??
            ## save 
            #state=hunter.save_state()
            #print(state)
            #result=[]
            #result.append(state["global_Pval"])
            #result.append(state["significance"])
            #print("Result=",result) 

            if Zval>ExpectedLocalZvalue:
                    BumpFound=True 
                    print("Bump with local Z=","{:.1f}".format(Zval)," and pos=",sign_center," chan=",channel," T=",TRIG_TYPE)
                    print(hunter.bump_info(data_y))
                    # collect bumps for visual inspection 
                    Bump=[hback,back]
                    BumpCollector.append(Bump)
    if (BumpFound): NrFound= NrFound+1

# the probability that background fluctuations alone (the null hypothesis) could produce a 
# result as extreme as, or more extreme than, the observed experimental data
print()
print("###### RESULT ###### ")
print(" Total events requested =", MaxEvents)
print(" Found events with bumps=", NrFound)
pvalue=float(NrFound)/MaxEvents
print(" Expected Z=",ExpectedLocalZvalue)
Zval=ROOT.RooStats.PValueToSignificance(pvalue)
print(" Found global p-value=",pvalue, " or Z =",Zval)
print("###### END RESULT ###### ")


figdir="figs/"
name=os.path.basename(__file__)
name=name.replace(".py","")
epsfig=figdir+name+".eps"
ps1 = TPostScript( epsfig,113)


# Plot one histogram for debugging
TRIG_TYPE=PLOT_TRIGGER_TYPE 
channel=PLOT_CHANNEL 
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


# write into file
rootfile="figs/bumps.root"
hfile=TFile(rootfile,"RECREATE","signatures")
for bum  in range(len( BumpCollector ) ):
    peak=BumpCollector[bum]

    hhh=peak[0].GetTitle()
    fun=peak[1].GetTitle()
   
    peak[0].SetTitle(peak[0].GetTitle()+"_bump"+str(bum))
    peak[0].SetName(peak[0].GetName()+"_bump"+str(bum))
   
    peak[1].SetTitle(peak[1].GetTitle()+"_bump"+str(bum))
    peak[1].SetName(peak[1].GetName()+"_bump"+str(bum))

    peak[0].Write()
    peak[1].Write()
hfile.Close()
print("Write ",rootfile, " with bumps above Z=",ExpectedLocalZvalue)


# only for debugging.. Make it fluctuate 
if (len(FluctuateBin)>0):
            print("We fluctuated events in overlap! Check FluctuateBin")

print (epsfig)
ps1.Close()

c1.Update()
if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);


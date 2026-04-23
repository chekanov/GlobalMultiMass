# Estimation of global statistical significance for multiple histograms using a toy pseudoexperiment
# You may adjust pyBumpHunter since it uses this package..
#  Authors: 
#  Sergei V.Chekanov (ANL)
#  Edison J. Weik 

############### USER SETTING ######################

## Expected local significance as in BumpHunter
ExpectedLocalZvalue=5 


# Maximum pseudo-experiments for global p-value
# for 6 sigma Z (local), set to maximum value to  a large number
MaxEvents=100

# Do multiplicative factor for function 
MinEntries = 0.0011 
# Do multiplicative factor for function 
# MinEntries = 11


# Do you want also re-fit template with the same function?
# This will be much slower, use it for the production mode.
fitAgain=False 


# CM energy for fit functions
cms=13000.


# THIS IS FOR PLOTTING and debugging.. 
# NOT REFLECTED ON P-Value 
PLOT_TRIGGER_TYPE=1 
PLOT_CHANNEL="jj"


############# END USER SETTING ######################


import argparse
# ----------------- Command line argument parser -----------------
parser = argparse.ArgumentParser(description="Run multiple pseudo-experiments for jet+X masses.")

parser.add_argument("--ExpectedLocalZvalue", type=float, default=ExpectedLocalZvalue,
                    help="Expected local significance as in BumpHunter (default: 7)")
parser.add_argument("--MaxEvents", type=int, default=MaxEvents,
                    help="Maximum pseudo-experiments for global p-value (default: 1000)")
parser.add_argument("--MinEntries", type=float, default=MinEntries,
                    help="Float value to multiply p0 in function function")
parser.add_argument("--doFit", type=lambda x: (str(x).lower() == 'true'), default=fitAgain,
                    help="Set to True to refit random template data. Will suppress fluctuations for low statistics histograms! (default: False)")
parser.add_argument("--interactive", type=lambda x: (str(x).lower() == 'true'), default=True,
                    help="Set to True to show interactive canvas (default: True)")

args = parser.parse_args()

# ----------------- Use arguments in your code -----------------
ExpectedLocalZvalue = float(args.ExpectedLocalZvalue) 
MaxEvents = int(args.MaxEvents) 
MinEntries = float(args.MinEntries) 
isInteractive=bool(args.interactive) 
fitAgain=bool(args.doFit)

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
from ROOT import TCanvas, TVectorD, TPostScript, TFile, TLegend, gPad, TF1, TRandom3, TH1D, TMath
from array import array
import psutil

process = psutil.Process(os.getpid())


# Current 28 configuration
CHANNELS = [
    "jj",
    "jb",
    "je+",
    "je-",
    "jm+",
    "jm-",
    "jg",
    "bb",
    "be+",
    "be-",
    "bm+",
    "bm-",
    "bg",
    "e+e+",
    "e+e-",
    "e+m+",
    "e+m-",
    "e+g",
    "e-e-",
    "e-m+",
    "e-m-",
    "e-g",
    "m+m+",
    "m+m-",
    "m+g",
    "m-m-",
    "m-g",
    "gg",
]

# 280 numbers
# CHANNELS = ["jj"]+[str(i) for i in range(1, 280)]


print("############## START ################")
print("The number of pseudo-experiments=",MaxEvents)
print("Searching for bumps with Z=",ExpectedLocalZvalue," which is ",z_to_p_value(ExpectedLocalZvalue)," p-value")
print("Float value to multiply p0=",MinEntries)
print("Is interactive mode ? =",isInteractive)
print("Do you refit pseudo-data? =",fitAgain)


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


# -----------------------------
# Main loop to create backgrond function and make random histograms 
# -----------------------------
# map with collection of histograms
Histograms={}
BackgroudFunction={}

# use it for debug
# CHANNELS=["jj"]

# collect bumps for further analysois
BumpCollector=[] 

## fill from the files parameters
mypar={}
myfunc={}
TRIG_TYPE=1
fit_min = 300
fit_max = 9000
print("Creating all functions ...")
for channel in CHANNELS:
        name = f"{TRIG_TYPE}_{channel}"
        parameters = [0.0016723174676486156, 0.5818046917532329, -6.57507710139123, -0.8003038477144173, 0.01]
        mback=FiveParam2015()
        back = TF1("back_" + name, mback, fit_min, fit_max, 5)
        mypar[name]=parameters
        myfunc[name]=back


channel="jj"
name = f"{TRIG_TYPE}_{channel}"
data_pars=mypar[name]
for i in range(len(parameters)):
                value = parameters[i]
                if (i == 0): value=MinEntries*value
                back.SetParameter(i, value)
                if value == 0.0:
                    back.FixParameter(i, value)
Sum=back.Integral(fit_min, fit_max);
print("Sum of integral=",Sum);

# sys.exit()


# make empty hemplate histogram
hbackJJ_name = f"histoJJ_{TRIG_TYPE}_{channel}"
hbackJJ = TH1D(hbackJJ_name, hbackJJ_name, len(mjjBins) - 1, mjjBins)
hbackJJ.SetTitle(hbackJJ_name)
hbackJJ.SetName(hbackJJ_name)
hbackJJ.SetDirectory(0)
                


NTOT=0
print("\n#######Loop over events",MaxEvents)
for event in range(MaxEvents):

    BumpFound=False # so far no bump found for this run

    if (event %50 == 0 and event>1 ): 
                    print("##  Pseudo Event=",event)
                    pvalue=float(NrFound)/NTOT
                    Zval=ROOT.RooStats.PValueToSignificance(pvalue)
                    mem = process.memory_info()
                    print("    So far global p-value=",pvalue, " or Z =",Zval, " evt found=",NrFound, " f(x) Integral=",Sum)
                    print(f"   RSS={mem.rss / 1024**2:.1f} MB", f"VMS={mem.vms / 1024**2:.1f} MB")


    # loop over 1 triggers
    if TRIG_TYPE==1:

        XMIN=300
        TLABEL = "jj"


        # Loop over each mass for a channel
        for channel in CHANNELS:

            # get the function from the map of functions
            name = f"{TRIG_TYPE}_{channel}"
            if (name in myfunc):  back = myfunc[name]
            else:
                 print("Cannot find the function",name) 
                 continue
         
            for i in range(len(parameters)):
                value = parameters[i]
                back.SetParameter(i, value)
                if value == 0.0:
                    back.FixParameter(i, value)


            # -----------------------------------------------------------------
            # CLEAN REWRITE: Single-pass loop for all bin fluctuations
            # -----------------------------------------------------------------
            hback_name = f"histo_{TRIG_TYPE}_{channel}"
            hback = hbackJJ.Clone()
            hback.SetTitle(hback_name)
            hback.SetName(hback_name)
            hback.SetDirectory(0)

         
            for i in range(hback.GetNbinsX() - 1):
                center = hback.GetBinCenter(i + 1)
                B_total = back.Eval(center)

               
                # Purely independent fluctuations
                pseudo_total = r.PoissonD(B_total) if B_total > 0 else 0

                hback.SetBinContent(i + 1, pseudo_total)
                hback.SetBinError(i + 1, TMath.Sqrt(pseudo_total) if pseudo_total > 0 else 0)

                

            ### finish. Fix some ranges and run BumpHunter
            TotalEvents=hback.Integral(hback.FindBin(fit_min), hback.FindBin(fit_max))

            

            # this will be slower, but could be more precise...
            if (fitAgain):
                      backFIT=back.Clone()
                      ChiMax=2.0
                      print("Fit event=",event, " T=",TRIG_TYPE, " ch=",channel)
                      fitr =  hback.Fit(backFIT,"ISMRQ0")
                      chi2ndf=100
                      if (backFIT.GetNDF()>0): chi2ndf=backFIT.GetChisquare()/backFIT.GetNDF()
                      if (chi2ndf > ChiMax): 
                         backFIT.Delete()
                         continue 
                      else:
                         back=backFIT.Clone()
                         backFIT.Delete()
                         print(" -> fit with chi2/ndf={:.2f} ".format(chi2ndf))

            # collect histograms
            BackgroudFunction[hback_name] = back  
            Histograms[hback_name] = hback
           

            # This is for debug only! 
            # if (channel != "jb" or TRIG_TYPE != 2): continue

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
            #p1=parameters[0]
            #p2=parameters[1]
            #p3=parameters[2]
            #p4=parameters[3]
            #p5=parameters[4]
            #integral_data = None
            #bkg_function_nom = FiveParam(cms, bkg_x_center, p1, p2, p3, p4, p5)
            #bkg_sample_nom, weight_nom = construct_bkg_sample(bkg_function_nom, bkg_x_center, integral_data)
            #bkg_sample_nom=bkg_function_nom

            # convert  ROOT background to NumPy array
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
            hunter = BH.BumpHunter1D( rang=None, width_min=2, width_max=30, width_step=1, scan_step=1, npe=10, seed=666, bins=hist_x)
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



            if Zval>ExpectedLocalZvalue:
                    # find bump postion 
                    x1,x2=0,0
                    if hasattr(hunter, "min_loc_ar") and hasattr(hunter, "min_width_ar"):
                       i = hunter.min_loc_ar[0] if np.ndim(hunter.min_loc_ar) else hunter.min_loc_ar
                       w = hunter.min_width_ar[0] if np.ndim(hunter.min_width_ar) else hunter.min_width_ar
                       # do not consider bumps at the tail 
                       if ((i+w)>len(bkg_x_center)-1):  continue
                       x1 = bkg_x_center[i]
                       x2 = bkg_x_center[i + w]
                    BumpFound=True
                    print("Bump with local Z=","{:.1f}".format(Zval)," and pos=",sign_center," chan=",channel," T=",TRIG_TYPE)
                    print(hunter.bump_info(data_y))
                    print("Best interval:", x1, x2 )
                    print("Best width:", x2-x1)
                    print("Record bump=",hback.GetTitle()) 
                    # collect bumps for visual inspection 
                    Bump=[hback,back, [x1, x2, Zval] ]
                    BumpCollector.append(Bump)

        # clean up
        try:
          
          del hback 
          del hunter
          del data_y
          del bkg_sample_nom
          del bkg_x_center
          del back
        except NameError:
             pass

    NTOT += 1
    if (BumpFound): NrFound += 1 
                 

# the probability that background fluctuations alone (the null hypothesis) could produce a 
# result as extreme as, or more extreme than, the observed experimental data
print()
print("###### RESULT ###### ")
print(" Total events requested =", NTOT)
print(" Found events with bumps=", NrFound)
pvalue=float(NrFound)/NTOT
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


leg2=TLegend(0.5, 0.55, 0.89, 0.84);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.04);
leg2.SetHeader("T="+str(TRIG_TYPE)+" "+str(channel))
leg2.AddEntry(hhBak,"Fit template","lp")
leg2.AddEntry(hhD,"Obtained sample","lp")
leg2.Draw("same");
hhD.Draw("pe same")


print (epsfig)
ps1.Close()

c1.Update()
if (isInteractive):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);


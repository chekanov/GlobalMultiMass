# Some useful functions..
# S.Chekanov (ANL)

import numpy as np

from math import log
import math,os,sys,json
import ROOT
from ROOT import TCanvas, TPostScript, TLegend, gPad, TF1, TRandom3, TH1D, TMath
from array import array
import uproot as upr  ## Used to read data from a root file


def getMaxNonzero(h1, xmin, ycut=0.5):
    """ Find last X value .."""
    xaxis = h1.GetXaxis()
    Ntot = xaxis.GetNbins()
    xmax = xmin
    for i in range(Ntot, 0, -1):
        y1 = h1.GetBinContent(i)
        x1 = h1.GetBinCenter(i)
        err = h1.GetBinWidth(i)
        if (x1 + err) < xmin:
            continue
        if y1 >= ycut:
            xmax = xaxis.GetBinUpEdge(i)
            break
    if (xmax - xmin < 100):
         print("XMIN and XMAX are too close!")
         h1.Print("All")

    return xmax



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


# ## Bkg sampling from function
#
# - Five parameter function form for dijet mass spectrum:
#
# $f(x) = p_1(1-x)^{p_2} x^{p_3+p_4\ln x + p_5\ln^2 x}$
#
# - alternative:
#
# $f(x) = p_1(1-x)^{p_2} x^{p_3+p_4\ln x + p_5/\sqrt{x}}$
#
# [Fit results](https://gitlab.cern.ch/dijetpluslepton/anomalydetection/ana/-/blob/master/figs_bh/ARpval_t2.txt?ref_type=heads#L8)
def FiveParam(Ecm, x_center, p1, p2, p3, p4, p5, bumphunter_implementation=False):
    #print("BumpHunter implementation=",bumphunter_implementation)
    #print('Use p5 with:',Ecm, p1, p2, p3, p4, p5)
    x = x_center / Ecm
    nlog=np.log(x)
    if bumphunter_implementation:
        fun = p1 * np.power((1.0 - x), p2) * 1.0 / np.power(x, (p3 + p4 * nlog +  p5 * nlog * nlog ))
    else:
        fun = p1 * np.power((1.0 - x), p2) * np.power(x, (p3 + p4 * nlog + p5 * nlog * nlog ))
    return fun

def FiveParam_alt(Ecm, x_center, p1, p2, p3, p4, p5):
    #print('Use p5_alt with:',Ecm, p1, p2, p3, p4, p5)
    x = x_center / Ecm
    nlog=np.log(x)
    fun = p1 * np.power((1.0 - x), p2) * np.power(x, (p3 + p4 * nlog + p5 / np.sqrt(x)))
    #fun = p1 * np.power((1.0 - x), p2) * 1 / np.power(x, (p3 + p4 * np.log(x)+  p5 * np.log(x) * np.log(x) ))
    return fun


# get input parameters: hist with data and ranges. Note that fit ranges taken from json files
def get_input(hist,fit_min,fit_max):
    # Extract bin contents and bin edges from a ROOT TH1 histogram
    nbins = hist.GetNbinsX()
    hist_y = np.array([hist.GetBinContent(i) for i in range(1, nbins + 1)])
    hist_x = np.array([hist.GetBinLowEdge(i) for i in range(1, nbins + 2)])
    #print('Initial bins:', hist_y.size)

    fit_range = (0,0)
    fit_range = (fit_min, fit_max)
    #print("fmin=",fit_min," fmax=",fit_max);

    hist_y = hist_y[(hist_x[:-1] >= fit_range[0]) & (hist_x[1:] <= fit_range[1])]
    hist_x = hist_x[(hist_x >= fit_range[0]) & (hist_x <= fit_range[1])]
    data_x_center = (hist_x[1:] + hist_x[:-1])/2
    data_bin_width = hist_x[1:] - hist_x[:-1]
    #print('Fit bins:', hist_y.size)
    data_y = np.repeat(data_x_center, hist_y.astype(int))
    #print('Obs data yields:', data_y.shape[0])
    return data_x_center, data_bin_width, data_y, hist_x


# Generate a sample of background points given a histogram and bins
def construct_bkg_sample(bkg_y_value, bkg_x_center, integral=None):
    '''
    bkg_y_value:
        background histogram (bin contents)
    bkg_x_center:
        bin center, bkg_x_center.size = bkg_x_center.size
    integral_data:
        yields to scale to (can scale to data yields?)
    '''

    # First generate points for each bkg bin to be bin content
    bkg_y_discrete = np.round(bkg_y_value * 0.1).astype(int)
    bkg_y_discrete[bkg_y_discrete == 0] = 1
    bkg_sample = np.repeat(bkg_x_center, bkg_y_discrete)

    # Scale to functional yields from discreted yields
    scale_bkg_discrete = bkg_y_value / bkg_y_discrete
    weights = np.repeat(scale_bkg_discrete, bkg_y_discrete)

    # Scale to integral, for properly computed p1 parameter, scaling is not expected
    if integral is not None:
        weights *= integral / np.sum(weights)
    else:
        print('Skip global scaling')
    return bkg_sample, weights


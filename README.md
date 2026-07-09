# Calculation of Global Significance for Multiple Invariant Masses Using Event Overlaps

This repository contains code to compute the global p-value (or global significance) in situations with multiple invariant mass distributions, including the possibility of event overlap between masses.
 
The code reads functional templates obtained from fits and performs multiple pseudo-experiments. In each pseudo-experiment, 63 histograms are randomly generated from the analytic functions. These histograms are organized into 9 different invariant masses across 7 independent trigger channels.

To better reproduce the behavior observed in real data, overlaps between invariant mass distributions within a given trigger stream can be introduced using the correlations observed in data. An alternative background hypothesis can also be enabled.

Once the random templates are generated for a pseudo-experiment, the code runs pyBumpHunter to search for statistically significant bumps above a specified local significance threshold Z. It then counts how many pseudo-experiments contain such significant local excesses.

Finally, using the fraction of pseudo-experiments that produce at least one significant bump, the code estimates the probability of observing such an excess under the background-only hypothesis, thereby determining the corresponding global p-value (or global Z-value).


## Running the Code

First, copy https://github.com/scikit-hep/pyBumpHunter to the external directory, and install it.

```bash
git clone https://github.com/scikit-hep/pyBumpHunter.git
cd pyBumpHunter
pip install .
```

Then you may adjust the statement "sys.path.append("pyBumpHunter/")" in the 1st line of "globalPvalue.py"

To run the analysis:

```bash
source setup.sh # Only if needed! 
python globalPvalue.py
```
If you use the command line arguments, use this: 

```bash
source setup.sh # Only if needed!
nohup python globalPvalue.py --ExpectedLocalZvalue 5 --MaxEvents 25000  --noOverlap false --interactive false --doFit false > globalPvalue5.log 2>&1 &
nohup python globalPvalue.py --ExpectedLocalZvalue 5 --MaxEvents 25000  --noOverlap true --interactive false --doFit false > globalPvalue5no.log 2>&1 &
```

which runs 25,000 pseudo-experiments with 9 invariant masses in 7 independent triggers, to estimate the 
global significance for the expected Z=5 sigma (local) from BumpHunter in any distribution. 
You may monitor the program by doing `grep 'So far' globalPvalue5.log`.


The file with most significant bumps is stored in "figs/bumps.root". Use the  showBumps.py to plot them for debugging.


## Benchmark results 

This table summarizes benchmark results. The first column gives the required local significance in any of the 63 histograms, for any width. The second column shows the corresponding global significance obtained. The third column shows the case with no mass overlaps, i.e. when all 63 histograms are treated as independent.


|Required local Z| Found global Z (overlap) | Found global Z (no overlap) |
|----------------|--------------------------|-----------------------------|
| 4              | 1.50  (p-value=0.066)    | 0.80 (p-value=0.213)        |
| 5              | 3.21  (p-value=0.00064)  | 2.67 (p-value=0.0037)       | 
| 6              | Need to run longer       | Need to run longer          |
| 7              | Need to run longer       | Need to run longer          | 

This table was created using using `"--doFit false"` option. 

The “overlap” case yields smaller global p-values (or larger global Z values) than the “no overlap” case. 
Thus the "no overlap" case represents the most conservative assumption, compared to the overlap case.
Generally, this is consistent with the observation that the toy experiment should yield a value of 
Z between two extremes: one in which all 63 histograms are treated as independent (leading to the smallest Z), 
and the other in which the 63 histograms are grouped into 7 independent groups, with 9 histograms in each group 
having 100% overlap (i.e., being identical repeated histograms).
The introduced correlation between overlapping part and non-overlapping part of the distribution 
is expected to be negative: when Poisson fluctuates  overlap fraction of the  histogram, the remaining part should go down to satisfy the expected template. This reduces fluctuations in the "overlap" case, leading to larger Z. 

 

*Note*: These results are very preliminary and are based on 100,000 pseudo-experiments. The uncertainty on the quoted Z-values is approximately ±0.1 (for the 7 sigma case).

*Checking*: 
If you reduce the number of probed masses, you should expect the global significance to increase.
To check the code, you may enable only one mass (`"jj"`) or two masses (`"jj"` and `"jb"`) in the list of masses on line 98, and only trigger 2, which has the largest statistics (see line 248 and change it to `range(2, 3)`). This will be similar to the "global" BumpHunter p-value.
You should set `--doFit false` since the statistics is good for these histograms.


|Required local Z| Global Z (only jj)      | Global Z (jj & jb)  |
|----------------|-------------------------|---------------------|  
| 4              | 2.35                    | 2.13                |
| 5              | 3.60                    | 3.51                | 

It shows the reduction of the global significance  with increase of the number of histograms.
The jj+jb case are done without overlap.

## Toy pseudoexperiments for independent distributions

This section describes the code used in the paper:

**S. V. Chekanov, E. Weik**   *On the Statistical Interpretation of Discoveries in LHC Data*, Preprint **HEP-ANL-203752**, May 16, 2026, arXiv:**2605.24441**

To reproduce the results presented in the paper, run:

```bash
python  globalDiscovery.py  --ExpectedLocalZvalue 3  --MaxEvents 20000  --interactive false --MinEntries  0.00009483  --doFit false > /dev/null 2>&1
```
This code runs tests for "3 sigma" local significances. The parameter "MinEntries" defines the normalization of the distribution (tune it as needed).  The deafult of this code assumes 280 independent distributions (on the line 127 for the list "CHANNELS"). For 28 masses, comment out the line 127. For just 1 distribution, just use CHANNELS=["jj"] after the line 127. The output of this program gives the global p-value.

## Additional Details

- The script reads the template configuration and performs the global p-value calculation based on the pseudo-experiment.
- The code also produces a **Jet+b plot**, assuming that **a fraction of the events originate from Jet+Jet processes**.

## Notes

- Review the implementation in `globalPvalue.py` to understand how the pseudo-experiment and template handling are performed.
- The number of pseudo-experiments and templates can be modified in the script if needed.

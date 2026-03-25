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
python globalPvalue.py --ExpectedLocalZvalue 5 --MaxEvents 10000  --noOverlap false
```

which runs 10,000 pseudo-experiments with 9 invariant masses in 7 independent triggers, to estimate the 
global significance for the expected Z=5 sigma (local) from BumpHunter in any distribution. 


The file with most significant bumps is stored in "figs/bumps.root". Use the  showBumps.py to plot them for debugging.


## Configuration

You can configure all input values at the beginning of the script **globalPvalue.py**. Make sure you run a sufficient number of experiments to address the 6 - 7  sigma requirement for local statistical deviation.

## Benchmark results 

This table summarizes benchmark results. The first column gives the required local significance in any of the 63 histograms, for any width. The second column shows the corresponding global significance obtained. The third column shows the case with no mass overlaps, i.e. when all 63 histograms are treated as independent.


|Required local Z| Found global Z (overlap) | Found global Z (no overlap) |
|----------------|--------------------------|-----------------------------|
| 3              | INF (p-value=1)          | INF (p-value=1)             |
| 5              | -0.08 (p-value= 0.54)    | 0.15 (p-value= 0.44)         | 
| 6              | 1.49  (p-value= 0.068)   | 1.64 (p-value=0.051)        |
| 7              | 2.61  (p-value= 0.0045)  | 2.62 (p-value=0.0042)       | 

This table was created using using "refit" option, which may suppress statistics fluctuations. But these results are likely to be more stable in the case of low statistics. 
If statistics of histogram is sufficient, use the option "--doFit false". 


The observed global Z value is somewhat smaller than that expected for a fluctuation of a single bin above the Z(local) threshold in 63 histograms with 100 bins each, due to Poisson statistics. This is because the BumpHunter only selects excesses with at least two adjacent bins fluctuating upward, which is more consistent with a physical signal having finite resolution.

The “overlap” case yields slightly more conservative p-values than the “no overlap” case. This is because the overlap was introduced with a negative correlation, which reduces fluctuations in the data points. Since the toy simulation does not model the exact event-by-event correlation, this represents the most conservative assumption. By contrast, a positive correlation would be expected to increase the significance relative to the “no overlap” (fully independent) case.


*Note*: These results are very preliminary and are based on 10,000 pseudo-experiments. The uncertainty on the quoted Z-values is approximately ±0.1 (for the 7 sigma case).

*Checking*: 
If you reduce the number of probed masses, you should expect the global significance to increase.

To check the code, you may enable only one mass (`"jj"`) or two masses (`"jj"` and `"jb"`) in the list of masses on line 98, and only trigger 2, which has the largest statistics (see line 248 and change it to `range(2, 3)`). This will be similar to the "global" BumpHunter p-value.
You should set `--doFit false` since the statistics is good for these histograms.


|Required local Z| Global Z (only jj)      | Global Z (jj & jb)  |
|----------------|---------- --------------|---------------------|  
| 5              | 1.58                    | 1.39                | 
| 6              | 2.58                    | 2.32                | 
| 7              | 3.52                    | 3.23                | 

It shows the reduction of the global significance  with increase of the number of histograms.
The jj and jb case are done without overlap.


## Additional Details

- The script reads the template configuration and performs the global p-value calculation based on the pseudo-experiment.
- The code also produces a **Jet+b plot**, assuming that **a fraction of the events originate from Jet+Jet processes**.

## Notes

- Review the implementation in `globalPvalue.py` to understand how the pseudo-experiment and template handling are performed.
- The number of pseudo-experiments and templates can be modified in the script if needed.

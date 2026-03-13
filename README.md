# Calculation of Global Significance for Jet+X Analysis Using Event Overlaps

This repository contains code to compute the global p-value (global significance) for the Jet+X analysis using multiple invariant mass distributions, including the possibility of event overlap between masses.

The code reads functional templates obtained from fits and performs multiple pseudo-experiments. In each pseudo-experiment, 63 histograms are randomly generated from the analytic functions. These histograms are organized into 9 different invariant masses across 7 independent trigger streams.

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

which runs 10000 pseudo-experiments with 9 invariant masses in 7 independent triggers, to estimate the 
global significance for the expected Z=5 sigma (local) from BumpHunter in any of the pseudo-experiment (and any distribution). 


The file with most significant bumps is stored in "figs/bumps.root". Use the  showBumps.py to plot them for debugging.


## Configuration

You can configure all input values at the beginning of the script **globalPvalue.py**. Make sure you run a sufficient number of experiments to address the 6 - 7  sigma requirement for local statistical deviation.
 

## Additional Details

- The script reads the template configuration and performs the global p-value calculation based on the pseudo-experiment.
- The code also produces a **Jet+b plot**, assuming that **a fraction of the events originate from Jet+Jet processes**.

## Notes

- Review the implementation in `globalPvalue.py` to understand how the pseudo-experiment and template handling are performed.
- The number of pseudo-experiments and templates can be modified in the script if needed.

# Calculation of Global Significance for Jet+X Analysis Using Overlapping PF Events

This repository contains code to compute the **global p-value (global significance)** for the Jet+X analysis using overlapping particle-flow (PF) events.

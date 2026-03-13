# Calculation of Global Significance for Jet+X Analysis Using Overlapping PF Events

This repository contains code to compute the **global p-value (global significance)** for the Jet+X analysis using overlapping particle-flow (PF) events.

The code reads functional templates from fits and then runs multiple pseudo-experiments. In each pseudo-experiment, 63 histograms are randomly generated from the analytic functions. It then runs pyBumpHunter to search for significant bumps above a specified significance threshold Z. Next, it counts how many pseudo-experiments contain such significant bumps. Finally, it computes the corresponding statistical significances.


## Running the Code

First, copy https://github.com/scikit-hep/pyBumpHunter to the exteranal directory, and install it.

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

The file with most significant bumps is stored in "figs/bumps.root". Use the  showBumps.py to plot them for debugging.


## Current Configuration

- The code currently generates **one pseudo-experiment**.
- A total of **63 templates** are used in the calculation.

## Additional Details

- The script reads the template configuration and performs the global p-value calculation based on the pseudo-experiment.
- The code also produces a **Jet+b plot**, assuming that **a fraction of the events originate from Jet+Jet processes**.

## Notes

- Review the implementation in `globalPvalue.py` to understand how the pseudo-experiment and template handling are performed.
- The number of pseudo-experiments and templates can be modified in the script if needed.

# Calculation of Global Significance for Jet+X Analysis Using Overlapping PF Events

This repository contains code to compute the **global p-value (global significance)** for the Jet+X analysis using overlapping particle-flow (PF) events.

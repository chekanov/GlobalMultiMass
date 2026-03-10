# Calculation of Global Significance for Jet+X Analysis Using Overlapping PF Events

This repository contains code to compute the **global p-value (global significance)** for the Jet+X analysis using overlapping particle-flow (PF) events.

## Running the Code

To run the analysis:

```bash
source setup.sh
python globalPvalue.py
```

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

#!/bin/bash

# 1. Setup Environment (Specific to ATLAS/CERN LCG views)
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

# 2. Default Values
EVENTS=100
ZVAL=5.0
CHIMAX=2.0
FIT_FLAG=""

# 3. Check for the "fit" argument
# Usage: ./run_with_fit.sh fit
if [ "$1" == "fit" ]; then
    echo "⚠️  Running WITH Chi-Squared Refitting (Production Mode - SLOW)"
    FIT_FLAG="--fit"
else
    echo "🚀 Running WITHOUT Refitting (Fast Mode)"
fi

# 4. Ensure directory exists
mkdir -p figs

# 5. Execute Python Comparison
python3 BumpHunter_Comparison.py \
    --events $EVENTS \
    --zval $ZVAL \
    --chimax $CHIMAX \
    $FIT_FLAG \
    -b | tee run_log.txt

echo "-------------------------------------------------------"
echo "✅ Run Complete. Summary printed below:"
# tail -n 12 run_log.txt

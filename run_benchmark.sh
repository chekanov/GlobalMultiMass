#!/bin/bash

# 1. Setup Environment
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

# USAGE
# ./run_test.sh (Fast, with overlap)
# ./run_test.sh fit (Slow, with overlap)
# ./run_test.sh no-overlap (Fast, no overlap)
# ./run_test.sh fit no-overlap (Slow, no overlap)

# 2. Default Values
EVENTS=100
ZVAL=5.0
CHIMAX=2.0
FIT_FLAG=""
OVERLAP_FLAG=""

# 3. Parse arguments
for arg in "$@"; do
    if [ "$arg" == "fit" ]; then
        echo "⚠️  Running WITH Chi-Squared Refitting (Production Mode - SLOW)"
        FIT_FLAG="--fit"
    elif [ "$arg" == "no-overlap" ]; then
        echo "⚠️  Running WITHOUT channel overlaps (Independent Toys)"
        OVERLAP_FLAG="--no-overlap"
    fi
done

# 4. Ensure directories exist
mkdir -p figs
mkdir -p fits

# 5. Execute Python Comparison
python3 BumpHunter_Comparison.py \
    --events $EVENTS \
    --zval $ZVAL \
    --chimax $CHIMAX \
    $FIT_FLAG \
    $OVERLAP_FLAG \
    -b | tee run_log.txt

echo "-------------------------------------------------------"
echo "✅ Run Complete. Summary printed below:"
tail -n 13 run_log.txt

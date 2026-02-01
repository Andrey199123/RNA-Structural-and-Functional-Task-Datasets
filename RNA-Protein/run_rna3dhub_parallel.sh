#!/bin/bash
# Run RNA 3D Hub processing in parallel
# Splits 13,782 RNAs into parallel jobs optimized for 12 CPU threads

TOTAL_RNAS=13782
NUM_JOBS=6  # 6 jobs with 2 threads each = 12 total threads
THREADS_PER_JOB=2
BATCH_SIZE=$((TOTAL_RNAS / NUM_JOBS))

echo "Processing $TOTAL_RNAS RNAs in $NUM_JOBS parallel jobs"
echo "Each job processes ~$BATCH_SIZE RNAs"
echo ""

# Kill any existing jobs
pkill -f "run_rna3dhub_processing.py" 2>/dev/null

# Start parallel jobs
for i in $(seq 0 $((NUM_JOBS - 1))); do
    START=$((i * BATCH_SIZE))
    if [ $i -eq $((NUM_JOBS - 1)) ]; then
        END=$TOTAL_RNAS
    else
        END=$(((i + 1) * BATCH_SIZE))
    fi
    
    echo "Starting job $((i+1))/$NUM_JOBS: RNAs $START to $END"
    
    python RNA-Protein/run_rna3dhub_processing.py \
        --start $START \
        --end $END \
        --batch-size 100 \
        --threads $THREADS_PER_JOB \
        > "RNA-Protein/data/rna3dhub_a3m_output/job_${i}.log" 2>&1 &
done

echo ""
echo "All jobs started! Monitor progress with:"
echo "  tail -f RNA-Protein/data/rna3dhub_a3m_output/job_*.log"
echo ""
echo "Check running jobs:"
echo "  ps aux | grep run_rna3dhub_processing"
echo ""
echo "Wait for all jobs to complete:"
echo "  wait"

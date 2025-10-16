#!/bin/bash
# ================================================================
#  FLO Prime Finder batch runner
#  Runs ./prime_finder for a specified number of seconds
#  and appends all stdout/stderr to results-batch.txt
# ================================================================

# ====== CONFIGURATION ======
EXEC="./prime_finder"
ARGS="--target-digits 2000 --flo-predict 1"
OUTFILE="results-batch.txt"

# ====== INPUT VALIDATION ======
if [ -z "$1" ]; then
  echo "Usage: $0 <duration_seconds>"
  exit 1
fi

DURATION=$1
START_TIME=$(date +%s)
END_TIME=$((START_TIME + DURATION))

echo "===============================================================" >> "$OUTFILE"
echo "Batch started at $(date) | Duration: ${DURATION}s" >> "$OUTFILE"
echo "Command: $EXEC $ARGS" >> "$OUTFILE"
echo "---------------------------------------------------------------" >> "$OUTFILE"

# ====== MAIN LOOP ======
RUN_COUNT=0
while [ "$(date +%s)" -lt "$END_TIME" ]; do
  RUN_COUNT=$((RUN_COUNT + 1))
  echo "[Run #$RUN_COUNT] $(date)" >> "$OUTFILE"
  $EXEC $ARGS >> "$OUTFILE" 2>&1
  echo "---------------------------------------------------------------" >> "$OUTFILE"
done

echo "Batch completed at $(date) after $RUN_COUNT runs." >> "$OUTFILE"
echo "===============================================================" >> "$OUTFILE"

# usage: ./run_batch.sh 3600 (for 1 hour of usage)
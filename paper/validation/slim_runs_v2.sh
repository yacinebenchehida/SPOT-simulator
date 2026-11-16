#!/usr/bin/env bash

###################################################################################
# Additional SLiM runs for the comparison between SPOT and SLiM, second version   #
# (Supporting Information S4)                                                     #
#                                                                                 #
#  (1) One more set of 2,000 neutral populations per neutral scenario (suffix _e, #
#      seed 5), which brings each neutral scenario to 10,000 SLiM populations for #
#      the check of the false-positive rate.                                      #
#  (2) 19 more independent runs (k = 2 to 20) of every scenario, with selection   #
#      and neutral, 2,000 populations each, written to                            #
#      paper/validation/slim_out/power/<name>_r<k>.tsv. With the original runs    #
#      (k = 1, seed 12345) they give 20 independent power estimates per scenario. #
#      Runs with selection use seed 1000 + k, neutral runs seed 2000 + k.         #
#                                                                                 #
# Run from the repository root:                                                   #
#   SLIM=/path/to/slim JOBS=2 bash paper/validation/slim_runs_v2.sh               #
# The haplodiploid_N2000 runs need about 4.2 GB of memory each: with less than    #
# about 9 GB, run them with JOBS=1 (the script skips the outputs that already     #
# exist and checks that each output is complete).                                 #
###################################################################################


###########################
# Stop at the first error #
###########################
set -euo pipefail


###################################
# Paths and number of populations #
###################################
# SLIM: path to the SLiM executable; OUT: output folder; R: populations simulated per run
SLIM=${SLIM:-slim}; OUT=paper/validation/slim_out; R=2000
# Create the output folders if needed
mkdir -p "$OUT/power"


#################################################################
# The scenarios, one per line: name system NF NM P0 S H         #
# SEL: scenarios with selection; NEU: neutral scenarios (S = 0) #
#################################################################
SEL="diploid_N1000_s0.02 diploid 500 500 0.1 0.02 0.5
diploid_N1000_s0.05_rec diploid 500 500 0.1 0.05 0
diploid_N200_s0.05 diploid 100 100 0.1 0.05 0.5
haploid_N1000_s0.02 haploid 1000 0 0.05 0.02 0.5
haplodiploid_N1000_s0.02 haplodiploid 500 500 0.1 0.02 0.5
haplodiploid_N2000_s0.05 haplodiploid 1000 1000 0.01 0.05 0.5"
NEU="diploid_N1000_neutral diploid 500 500 0.1 0 0.5
diploid_N200_neutral diploid 100 100 0.1 0 0.5
haploid_N1000_neutral haploid 1000 0 0.05 0 0.5
haplodiploid_N1000_neutral haplodiploid 500 500 0.1 0 0.5
haplodiploid_N2000_neutral haplodiploid 1000 1000 0.01 0 0.5"


#############################################################################
# Function: run SLiM for one scenario and check that the output is complete #
# Arguments: output file, name, system, NF, NM, P0, S, H, seed              #
#############################################################################
run_one() {   # out name system NF NM P0 S H seed
	# Skip the run if its output file already exists and is not empty
	[ -s "$1" ] && { echo "skip $1"; return; }
	# Remove a temporary file left by an interrupted run
	rm -f "$1.tmp"
	# Run SLiM into a temporary file, then check that it has 4 lines per population (one per sampling date)
	if "$SLIM" -s $9 -d "SYSTEM='$3'" -d NF=$4 -d NM=$5 -d P0=$6 -d S=$7 -d H=$8 -d R=$R \
		-d "GENS=c(0,20,40,90)" -d "NSAMP=c(20,20,20,40)" -d "OUT='$1.tmp'" paper/validation/slim_check.slim > /dev/null \
		&& [ "$(wc -l < "$1.tmp")" -eq $((4 * R)) ]; then
		# Complete: give the file its final name
		mv "$1.tmp" "$1"; echo "done $1 $(date +%T)"
	else
		# Incomplete or failed: remove the temporary file and report it
		rm -f "$1.tmp"; echo "FAILED $1 $(date +%T)"   # e.g. killed for lack of memory: rerun with JOBS=1
	fi
}


#################################################################
# Build the list of runs (one line each) and send it to run_one #
#################################################################
# Make the function and the variables available to the commands started by xargs
export -f run_one; export SLIM R
{
  # (1) One more neutral run per neutral scenario: suffix _e, seed 5
  echo "$NEU" | while read n rest; do echo "$OUT/${n}_e.tsv $n $rest 5"; done
  # (2) Runs 2 to 20 of every scenario
  for k in $(seq 2 20); do
    # With selection: seed 1000 + k
    echo "$SEL" | while read n rest; do echo "$OUT/power/${n}_r$k.tsv $n $rest $((1000+k))"; done
    # Neutral: seed 2000 + k
    echo "$NEU" | while read n rest; do echo "$OUT/power/${n}_r$k.tsv $n $rest $((2000+k))"; done
  done
} | xargs -P ${JOBS:-1} -L 1 bash -c 'run_one "$@"' _
# Message printed when all the runs are finished
echo ALLDONE

#!/usr/bin/env bash

####################################################################################
# SLiM side of the comparison between SPOT and SLiM (Supporting Information S4)    #
#                                                                                  #
# Each line of SCENARIOS below is one SLiM run: R = 2,000 populations of one       #
# scenario, simulated with the model paper/validation/slim_check.slim.             #
#                                                                                  #
# Requires SLiM 5 (https://messerlab.org/slim/; tested with SLiM 5.2).             #
# Run from the repository root (about 30 minutes; the largest runs need about 3 GB #
# of memory each):                                                                 #
#   bash paper/validation/slim_runs.sh                                             #
# then compare with SPOT:                                                          #
#   node paper/validation/slim_compare.js                                          #
#   (writes paper/validation/slim_results.tsv, slim_power.tsv and slim_fpr.tsv)    #
#                                                                                  #
# All runs share the design of the reference case: sampling at generations 0, 20,  #
# 40 and 90 (one generation per year: 1935, 1955, 1975, 2025), 20/20/20/40 females #
# per sampling date, perfect genotypes.                                            #
####################################################################################


###########################
# Stop at the first error #
###########################
set -euo pipefail


###################################
# Paths and number of populations #
###################################
SLIM=${SLIM:-slim}              # path to the SLiM executable (set SLIM=/path/to/slim if it is not on the PATH)
OUT=paper/validation/slim_out   # Output folder (one file per run)
R=2000                          # Populations simulated per run
mkdir -p "$OUT"                 # Create the output folder if needed


##################################################################################
# The scenarios, one per line                                                    #
# Columns: name  system  NF  NM  P0  S  H  [seed]                                #
#   name   : label used for the output file paper/validation/slim_out/<name>.tsv #
#            (and read by slim_compare.js)                                       #
#   system : diploid, haploid or haplodiploid                                    #
#   NF, NM : breeding females and males per generation (haploid: NF individuals; #
#            NM not used)                                                        #
#   P0     : starting frequency of the favoured allele                           #
#   S, H   : selection coefficient and dominance (S = 0 gives the neutral runs   #
#            used to set the threshold)                                          #
#   seed   : random seed (12345 if absent); the extra neutral runs (_b, _c, _d)  #
#            bring every neutral scenario to 8,000 populations, for a precise    #
#            check of the false-positive rate                                    #
##################################################################################
SCENARIOS="
diploid_N1000_s0.02         diploid       500  500  0.1   0.02  0.5
diploid_N1000_neutral       diploid       500  500  0.1   0     0.5
diploid_N1000_s0.05_rec     diploid       500  500  0.1   0.05  0
diploid_N200_s0.05          diploid       100  100  0.1   0.05  0.5
diploid_N200_neutral        diploid       100  100  0.1   0     0.5
haploid_N1000_s0.02         haploid       1000 0    0.05  0.02  0.5
haploid_N1000_neutral       haploid       1000 0    0.05  0     0.5
haplodiploid_N1000_s0.02    haplodiploid  500  500  0.1   0.02  0.5
haplodiploid_N1000_neutral  haplodiploid  500  500  0.1   0     0.5
haplodiploid_N2000_s0.05    haplodiploid  1000 1000 0.01  0.05  0.5
haplodiploid_N2000_neutral  haplodiploid  1000 1000 0.01  0     0.5
diploid_N1000_neutral_b     diploid       500  500  0.1   0     0.5  2
diploid_N1000_neutral_c     diploid       500  500  0.1   0     0.5  3
diploid_N1000_neutral_d     diploid       500  500  0.1   0     0.5  4
diploid_N200_neutral_b      diploid       100  100  0.1   0     0.5  2
diploid_N200_neutral_c      diploid       100  100  0.1   0     0.5  3
diploid_N200_neutral_d      diploid       100  100  0.1   0     0.5  4
haplodiploid_N1000_neutral_b haplodiploid 500  500  0.1   0     0.5  2
haplodiploid_N1000_neutral_c haplodiploid 500  500  0.1   0     0.5  3
haplodiploid_N1000_neutral_d haplodiploid 500  500  0.1   0     0.5  4
haploid_N1000_neutral_b     haploid       1000 0    0.05  0     0.5  2
haploid_N1000_neutral_c     haploid       1000 0    0.05  0     0.5  3
haploid_N1000_neutral_d     haploid       1000 0    0.05  0     0.5  4
haplodiploid_N2000_neutral_b haplodiploid 1000 1000 0.01  0     0.5  2
haplodiploid_N2000_neutral_c haplodiploid 1000 1000 0.01  0     0.5  3
haplodiploid_N2000_neutral_d haplodiploid 1000 1000 0.01  0     0.5  4
"


##############################################
# Function: run SLiM for one scenario        #
# Arguments: name system NF NM P0 S H [seed] #
##############################################
run_one() {   # arguments: name system NF NM P0 S H [seed]
	# Remove the output of a previous run (SLiM adds lines to the file)
	rm -f "$OUT/$1.tsv"
	# Run SLiM: -s gives the seed (12345 if the 8th argument is absent), each -d gives one setting of the model
	"$SLIM" -s ${8:-12345} -d "SYSTEM='$2'" -d NF=$3 -d NM=$4 -d P0=$5 -d S=$6 -d H=$7 -d R=$R \
		-d "GENS=c(0,20,40,90)" -d "NSAMP=c(20,20,20,40)" -d "OUT='$OUT/$1.tsv'" \
		paper/validation/slim_check.slim > /dev/null
	# Progress message
	echo "done $1"
}


##################################################################################
# Run all the scenarios                                                          #
# One run at a time by default; "JOBS=2 bash paper/validation/slim_runs.sh" runs #
# two at once if the memory allows                                               #
##################################################################################
# Make the function and the variables available to the commands started by xargs
export -f run_one; export SLIM OUT R
# one run at a time by default (JOBS=2 bash paper/validation/slim_runs.sh runs two at once if memory allows)
# Send each non-empty line of SCENARIOS to run_one (xargs -L 1: one line per call; -P: number of runs at once)
echo "$SCENARIOS" | grep -v '^$' | xargs -P ${JOBS:-1} -L 1 bash -c 'run_one "$@"' _

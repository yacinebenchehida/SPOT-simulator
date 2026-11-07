#!/usr/bin/env bash
# Additional SLiM runs for the SLiM comparison (v2, round numbers):
#  (1) one more set of 2,000 neutral replicates per neutral scenario (suffix _e, seed 5), bringing each neutral
#      scenario to 10,000 SLiM replicates for the false-positive check;
#  (2) 19 further independent runs (k = 2..20) of every selection and neutral scenario, 2,000 replicates each,
#      written to paper/validation/slim_out/power/<name>_r<k>.tsv; together with the original runs (k = 1, seed 12345)
#      they give 20 independent power estimates per scenario. Selection runs use seed 1000 + k, neutral runs 2000 + k.
# Run from the repository root: SLIM=/path/to/slim JOBS=2 bash paper/validation/slim_runs_v2.sh
# The haplodiploid_N2000 runs need about 4.2 GB of memory each: with less than about 9 GB, run them with JOBS=1
# (the script skips outputs that already exist and checks that each output is complete).
set -euo pipefail
SLIM=${SLIM:-slim}; OUT=paper/validation/slim_out; R=2000
mkdir -p "$OUT/power"
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
run_one() {   # out name system NF NM P0 S H seed
	[ -s "$1" ] && { echo "skip $1"; return; }
	rm -f "$1.tmp"
	if "$SLIM" -s $9 -d "SYSTEM='$3'" -d NF=$4 -d NM=$5 -d P0=$6 -d S=$7 -d H=$8 -d R=$R \
		-d "GENS=c(0,20,40,90)" -d "NSAMP=c(20,20,20,40)" -d "OUT='$1.tmp'" paper/validation/slim_check.slim > /dev/null \
		&& [ "$(wc -l < "$1.tmp")" -eq $((4 * R)) ]; then
		mv "$1.tmp" "$1"; echo "done $1 $(date +%T)"
	else
		rm -f "$1.tmp"; echo "FAILED $1 $(date +%T)"   # e.g. killed for lack of memory: rerun with JOBS=1
	fi
}
export -f run_one; export SLIM R
{
  echo "$NEU" | while read n rest; do echo "$OUT/${n}_e.tsv $n $rest 5"; done
  for k in $(seq 2 20); do
    echo "$SEL" | while read n rest; do echo "$OUT/power/${n}_r$k.tsv $n $rest $((1000+k))"; done
    echo "$NEU" | while read n rest; do echo "$OUT/power/${n}_r$k.tsv $n $rest $((2000+k))"; done
  done
} | xargs -P ${JOBS:-1} -L 1 bash -c 'run_one "$@"' _
echo ALLDONE

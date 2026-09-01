#!/bin/bash
#BSUB -J test_80cores
#BSUB -q standard
#BSUB -n 84
#BSUB -R "span[ptile=20]"
#BSUB -o test_80cores.%J.out
#BSUB -e test_80cores.%J.err
#BSUB -W 00:05

# Quick check: can the "standard" queue actually reserve > 80 slots when the
# job is allowed to spread across multiple hosts?
#
# - `-n 84`                 : request 84 slots (> 80).
# - `span[ptile=20]`        : pack up to 20 slots per host, so 84 slots spread
#                             across ~5 hosts. This does NOT force one host.
#                             (Use `span[hosts=1]` ONLY if you needed a single
#                             node -- we explicitly do not here.)
#
# Submit with:   bsub < test_lsf_80cores.sh
# Then read the .out file. If the job never starts (stays PEND), run the
# diagnostics at the bottom of this file to see why.

echo "==================== LSF allocation report ===================="
echo "Job ID              : ${LSB_JOBID}"
echo "Queue               : ${LSB_QUEUE}"
echo "Submit host         : ${LSB_SUB_HOST}"
echo "Requested slots (-n): ${LSB_DJOB_NUMPROC}"
echo

# LSB_MCPU_HOSTS is "host1 nslots1 host2 nslots2 ..."; LSB_HOSTS repeats each
# host once per slot. Parse both to report distinct hosts and total slots.
echo "LSB_MCPU_HOSTS      : ${LSB_MCPU_HOSTS}"
echo "LSB_HOSTS           : ${LSB_HOSTS}"
echo

total_slots=$(echo "${LSB_HOSTS}" | wc -w)
distinct_hosts=$(echo "${LSB_HOSTS}" | tr ' ' '\n' | sort -u | grep -c .)

echo "Total slots granted : ${total_slots}"
echo "Distinct hosts      : ${distinct_hosts}"
echo
echo "Per-host slot counts:"
echo "${LSB_HOSTS}" | tr ' ' '\n' | sort | uniq -c | awk '{printf "  %-25s %s slots\n", $2, $1}'
echo

if [ "${total_slots}" -ge 80 ]; then
    echo "RESULT: SUCCESS - reserved ${total_slots} slots (>= 80) across ${distinct_hosts} host(s)."
else
    echo "RESULT: ONLY ${total_slots} slots granted (< 80). See diagnostics below."
fi
echo "==============================================================="

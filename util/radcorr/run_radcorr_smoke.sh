#!/usr/bin/env bash
# Run short integration/timing checks without modifying tracked input files.
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
trials=${1:-100000}
created_inputs=()

cleanup() {
  local input
  for input in "${created_inputs[@]}"; do
    rm -f "$input"
  done
  rm -f "$repo_root"/worksim/radcorr_smoke_*3HeFit.bin
}
trap cleanup EXIT

printf '%-10s %8s %12s %10s %12s %12s %12s %12s\n' \
  'config' 'trials' 'rc_requests' 'failures' 'min_R' 'max_R' 'mean_R' 'runtime_s'

run_case() {
  local sf=$1 rad=$2 label=$3 tag input log time_file
  tag="radcorr_smoke_${label}"
  input="$repo_root/infiles/${tag}.inp"
  log="$repo_root/runout/${tag}.out"
  time_file="$repo_root/runout/${tag}.time"
  created_inputs+=("$input")

  awk -v trials="$trials" -v sf="$sf" -v rad="$rad" '
    NR == 4 { printf "%10d      Monte-Carlo trials\n", trials; next }
    NR == 8 { print "    -10.0       M.C. DP/P  down limit"; next }
    NR == 9 { print "     10.0       M.C. DP/P  up   limit"; next }
    /SF model flag/ { print "      1          SF model flag (0=F1F2IN21, 1=3HeFit table)"; next }
    /3He fit model/ { printf "      %d          3He fit model (1=SF1, 2=SF2, 3=SF3, 4=SF4, 5=SF5)\n", sf; next }
    /Radiative correction/ { printf "      %d          Radiative correction (0=Born, 1=radiated weight)\n", rad; next }
    { print }
  ' "$repo_root/infiles/kin3_hms.inp" > "$input"

  (
    cd "$repo_root/src"
    /usr/bin/time -f '%e' -o "$time_file" ./mc_single_arm > "$log" 2>&1 <<EOF
$tag
EOF
  )

  local requests failures minmax runtime min_r max_r mean_r
  requests=$(awk -F= '/3He radiative lookup requests/{gsub(/ /,"",$2); print $2}' "$log" | tail -n 1)
  failures=$(awk -F= '/3He radiative lookup failures/{gsub(/ /,"",$2); print $2}' "$log" | tail -n 1)
  minmax=$(awk -F= '/3He rad_weight_factor min\/max\/mean/{print $2}' "$log" | tail -n 1)
  runtime=$(tr -d '[:space:]' < "$time_file")
  if [[ "$rad" == "0" ]]; then
    requests='-'; failures='-'; min_r='-'; max_r='-'; mean_r='-'
  else
    read -r min_r max_r mean_r <<< "$minmax"
  fi
  printf '%-10s %8s %12s %10s %12s %12s %12s %12s\n' \
    "$label" "$trials" "${requests:--}" "${failures:--}" \
    "$min_r" "$max_r" "$mean_r" "$runtime"
}

run_case 1 0 SF1_Born
run_case 1 1 SF1_Rad
run_case 2 1 SF2_Rad
run_case 3 0 SF3_Born
run_case 3 1 SF3_Rad
run_case 4 1 SF4_Rad
run_case 5 1 SF5_Rad

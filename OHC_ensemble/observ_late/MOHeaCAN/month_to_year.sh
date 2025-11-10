#!/usr/bin/env bash
set -euo pipefail

IN="OHC-EEI_199301_202205_v5-0.nc"
VARS="gohc_16perc" #,gohc_16perc,gohc_84perc"
OUTNC="OHC-EEI_annual_v5-0.nc"
OUTCSV="OHC-EEI_annual_v5-0.csv"

# how many time slices?
NT=$(ncks --trd -m "$IN" | sed -nE 's/.*time.*size = ([0-9]+).*/\1/p' | head -n1)
# first year from the first timestamp (any CF calendar)
#Y0=$(ncap2 -O -v -s 'print(int(year(time(0))))' "$IN" | awk 'NF{print $NF}')
Y0=1993

# how many full 12-step years? (drop a partial tail cleanly)
NY=$(( NT / 12 ))
echo $NY
TMPD=$(mktemp -d -t ann_ohc.XXXX); trap 'rm -rf "$TMPD"' EXIT

for y in $(seq 0 $((NY-1))); do
  s=$(( y*12 ))
  e=$(( s+11 ))
  # mean over a 12-step block
  ncra -O -y avg -v "$VARS" -d time,$s,$e "$IN" "$TMPD/ann_$y.nc"
  # tag its calendar year (Y0, Y0+1, …)
  ncap2 -O -s "year=$(($Y0 + y))" "$TMPD/ann_$y.nc" "$TMPD/ann_$y.nc"
done

# single concat to annual file (no warnings)
ncrcat -O "$TMPD"/ann_*.nc "$OUTNC"

# CSV (year,gohc,gohc_16perc,gohc_84perc)
printf "year,gohc,gohc_16perc,gohc_84perc\n" > "$OUTCSV"
ncks -H -C -v year,$VARS -s "%g,%.*g,%.*g,%.*g\n" "$OUTNC" >> "$OUTCSV"

echo "Wrote: $OUTNC"
echo "Wrote: $OUTCSV"


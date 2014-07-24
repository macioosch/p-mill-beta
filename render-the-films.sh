#!/usr/bin/env bash
time ./animate-3.py output/steel-rosenkratz-mu_s-1.2-mu_r-0.5.csv plots/mu_s-1.2-mu_r-0.5.mp4
echo "mu_s-1.2-mu_r-0.5 done!"
for N in {3,5,7}0; do
    time ./animate-3.py output/steel-rosenkratz-N-${N}.csv plots/N-${N}.mp4
    echo "N-${N} done!"
done

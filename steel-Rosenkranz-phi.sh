#/usr/bin/env bash
for N in 30 50 70; do
    echo "Processing N = ${N}..."
    time ./p-mill-beta     \
        --type  3          \
        --dt    2e-6       \
        --tmax  2e-0       \
        --output_lines 1e4 \
        --e     0.70       \
        --E     2e11       \
        --G     8e10       \
        --mu_r  0.01       \
        --mu_s  0.20       \
        --rho   7.8e3      \
        --rb     5e-3      \
        --rc     5e-2      \
        --rs    15e-2      \
        --wc     62.7      \
        --ws    -20.9      \
        --N     $N         \
        > output/steel-rosenkratz-N-${N}.csv
done

# (1, 2, 2.5, 3) * 20.9 = (20.9, 41.8, 52.25, 62.7)

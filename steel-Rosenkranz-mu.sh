#/usr/bin/env bash
for mu_s in 0.2 0.35 0.5 1.2; do
    for mu_r in 0.01 0.1 0.5; do 
        echo "Processing: mu_s = ${mu_s}, mu_r = ${mu_r}..."
        time ./p-mill-beta     \
            --type  3          \
            --dt    2e-6       \
            --tmax  2e-0       \
            --output_lines 1e4 \
            --e     0.70       \
            --E     2e11       \
            --G     8e10       \
            --mu_r  $mu_r      \
            --mu_s  $mu_s      \
            --rho   7.8e3      \
            --rb     5e-3      \
            --rc     5e-2      \
            --rs    15e-2      \
            --wc     62.7      \
            --ws    -20.9      \
            --N     30         \
            > output/steel-rosenkratz-mu_s-${mu_s}-mu_r-${mu_r}.csv
    done
done

# (1, 2, 2.5, 3) * 20.9 = (20.9, 41.8, 52.25, 62.7)

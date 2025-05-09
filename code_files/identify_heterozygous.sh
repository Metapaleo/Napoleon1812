#!/bin/bash

# Usage: zcat SAMPLE_merged.rmdup.rescaled.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$5 != "."' | ./identify_heterozygous.sh > output.txt
# To see results: cut -f 13 output.txt | sort | uniq -c | awk '{print $2"\t"$1}'

while IFS= read -r line; do
    last_column=$(echo "$line" | awk -F'\t' '{print $NF}')
    values=$(echo "$last_column" | awk -F':' '{print $2}')
    value1=$(echo "$values" | cut -d',' -f1)
    value2=$(echo "$values" | cut -d',' -f2)

    # Calculate the condition using awk and store it in a variable
    condition=$(awk -v val1="$value1" -v val2="$value2" 'BEGIN {print (val1 > val2 && (val1 / (val2 + val1) > 0.8)) || (val2 > val1 && (val2 / (val2 + val1) > 0.8))}')

    # Check conditions and set category and X value
    if [ "$value1" -eq 0 ] || [ "$value2" -eq 0 ]; then
        if [ "$value1" -ne 0 ] || [ "$value2" -ne 0 ]; then
            X=1
        else
            X=0
        fi
        category="Homozygous"
    elif [ "$value1" -ne 0 ] && [ "$value2" -ne 0 ] && ((value1 + value2 >= 8)) && [ "$condition" -eq 1 ]; then
        ratio=$(awk -v val1="$value1" -v val2="$value2" 'BEGIN {print (val1 > val2) ? val1 / (val2 + val1) : val2 / (val2 + val1)}')
        category="LikelyHomozygous"
        X=$ratio
    elif [ "$value1" -ne 0 ] && [ "$value2" -ne 0 ] && ((value1 + value2 >= 8)) && [ "$condition" -eq 0 ]; then
        ratio=$(awk -v val1="$value1" -v val2="$value2" 'BEGIN {print (val1 > val2) ? val1 / (val2 + val1) : val2 / (val2 + val1)}')
        category="LikelyHeterozygous"
        X=$ratio
    elif [ "$value1" -ne 0 ] && [ "$value2" -ne 0 ] && ((value1 + value2 <= 8)); then
        ratio=$(awk -v val1="$value1" -v val2="$value2" 'BEGIN {print (val1 > val2) ? val1 / (val2 + val1) : val2 / (val2 + val1)}')
        category="LowCoverage"
        X=$ratio
    fi

    # Print the modified line
    echo -e "$line\t$value1\t$value2\t$category\t$X"
done

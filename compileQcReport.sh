#!/bin/bash
set -euo pipefail

# Christopher Medway AWMGS
# compiles a file of useful QC metrics from the multitude of PICARD metrics



# check FASTQC output
countQCFlagFails() {
    #count how many core FASTQC tests failed
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}


# make empty combined QC file
echo -e "Sample\tFastQC" > combined_QC.txt


# loop through each sample and make QC file
for sampleId in $(cat completed_samples.txt); do



        dir=./"$sampleId"/FastQC/


        # check all fastqc files for any fails
        rawSequenceQuality=PASS

        for report in $dir/"$sampleId"_*_fastqc.txt;
        do
            if [ $(countQCFlagFails $report) -gt 0 ]; then
                rawSequenceQuality=FAIL
            fi
        done

        # add to combined QC file
        echo -e "$sampleId\t$rawSequenceQuality" >> combined_QC.txt

    fi

done

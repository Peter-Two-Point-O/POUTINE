#!/bin/bash
#
# $1 := main output filename
# $2 := maxT pvalue threshold
#
# Will likely write a full program later to help users process results from Homoplasy_Counter.
#
# Output:
# 1)  .all_associations := same as full homoplasy counter output file minus all the log information
# 2)  .a1_sorted := each row represents one a1 allele that is homoplasically informative (always a subset of tot # homoplasically informative sites)
# 3)  .a2_sorted := each row represents one a2 allele that is homoplasically informative (always a subset of tot # homoplasically informative sites)
# 4)  .a1_top_hits_maxT := a1 associations filtered by maxT pvalue
# 5)  .a2_top_hits_maxT := a2 associations filtered by maxT pvalue


if [ $# -eq 2 ]; then

    # separate file for all associations (main output file minus all the log information)
    awk -F "\t" '{if (NF == 19) print $0}' $1 > $1.all_associations

    # a1 associations sorted by maxT
    awk -F "\t" '{if ($18 != "NaN") print $0}' $1.all_associations | sort -t $'\t' -gk18 > $1.a1_sorted
    
    # a2 associations sorted by maxT
    awk -F "\t" '{if ($19 != "NaN") print $0}' $1.all_associations | sort -t $'\t' -gk19 > $1.a2_sorted

    # a1 top hits by maxT <= $2
    awk -F "\t" -v MAXT_THRESHOLD=$2 '{if (FNR == 1) print $0; if ($18 != "NaN" && $18 <= MAXT_THRESHOLD) print $0}' $1.a1_sorted | sort -t $'\t' -gk18 > $1.a1_top_hits_maxT 
    
    # a2 top hits by maxT <= $2
    awk -F "\t" -v MAXT_THRESHOLD=$2 '{if (FNR == 1) print $0; if ($19 != "NaN" && $19 <= MAXT_THRESHOLD) print $0}' $1.a2_sorted | sort -t $'\t' -gk19 > $1.a2_top_hits_maxT     
else
    echo "Usage:  consume_results.sh RESULTS_FILENAME MAXT_THRESHOLD"
fi

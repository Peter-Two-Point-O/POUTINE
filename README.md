# Easy-Is-Better-Than-Better
Pre-Alpha private repo for a new homoplasy counter.  This pre-alpha release is intended for a handful of users to test as major changes continue to be incorporated.



## Installation
This program requires: (links are provided if you don't already have it on your box or want the latest version)

- java (https://www.oracle.com/java/technologies/javase-downloads.html)
- python (https://www.python.org/downloads/)
- treetime (https://treetime.readthedocs.io/en/latest/installation.html)



## Executing The Program
For this pre-alpha release, a convenience homoplasy_counter.sh file is included to run the program.  Future versions will likely be a precompiled binary with full command-line switch options like any other Unix binary (e.g. less, top, etc).

From inside the pre-alpha/ directory, simply execute: `./homoplasy_counter.sh &> output_filename.out &`

The `&>` redirects all terminal streams to the output file, while the final `&` backgrounds your process.  Both are optional.

To get a sorted list of top hits, simply run: `consume_results.sh output_filename.out 0.05`

`0.05` is the maxT p value cut-off user chooses.  This is a family-wise error rate so 5% or 1% are classically reasonable.  The script will also output the following files:

1)  `.all_associations`:  same as full homoplasy counter output file minus all the log information
2)  `.a1_sorted`:  each row represents one a1 allele in sorted fashion (a1 subset from .all_associations)
3)  `.a2_sorted`:  each row represents one a2 allele in sorted fashion (a2 subset from .all_associations)
4)  `.a1_top_hits_maxT`:  a1 associations filtered by maxT pvalue
5)  `.a2_top_hits_maxT`:  a2 associations filtered by maxT pvalue



## Reference Set
Included in this release is our *reference set* comprised of 124 mtb strains we have used for internal testing.  The homoplasy_counter.sh file starts a test run using this reference set.



## Current Major Bugs
- Concurrency bug during resampling:  This bug only occurs sometimes.  If your test run doesn't run to completion (`CLEAN EXIT` should be the last line in the output file) and you see the program hang, simply kill the process and restart it for now.



## Recent Major Code Changes
- Removed q values from the significance assessment.  Thus, all things R have been removed.  Main reason for this feature removal is because the resampling-derived FWER (maxT variant) is sufficient for users to sort and look for top hits (remember the philosophy here is "easier is better than better").  It's also more robust than anything in the FDR space in the sense that any filtering of segregating sites will alter the p value distribution that FDR-based schemes use to calculate the splines, and thus can dramatically alter the final q values.  So for example, I have currently set the program to ignore all sites with <= 3 homoplasic mutations as these will not show statistical significance (less sites mean faster run with lower memory).  The maxT p values are typically immune to changes from site filtering.
- Incorporated treetime for purposes of genotypic ancestral reconstruction using the default optimized joint probabilities method.  In a future release, when we are likely to add estimation statistics like a resampling-derived effect size + confidence intervals, we can reconsider if q values will benefit the user in any way.



## Upcoming Major Code Changes
- Incorporate vcf format.  This will also remove multi-fasta file and plink map file from standard input.
- Full command-line interface.  This will be our store-front! so will try to make this elegant and easy.
- Incorporate consume_results.sh code into main program along with other facilities to sort and consume results.
- Organize output (e.g. various results files, log file, etc)
- Logging code to facilitate prepared files for users to send us for trouble-shooting.
- Platform independence (mostly for Windows compatibility)
- Final code cleanup (e.g. comments, refactoring, etc)
- Deployment:  Likely to compile entire program (minus python component) into binaries, one for each platform.  This removes the java dependency, and users will simply run the program like any other binary.  Future work:  consider replacing python component to make this program completely self-contained; no installation, simply download the binary and run!



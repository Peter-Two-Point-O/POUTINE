# Easy-Is-Better-Than-Better
Pre-Alpha private repo for a new homoplasy counter.  This pre-alpha release is intended for a handful of users to test as major changes continue to be incorporated.  We named this initial version Easy-Is-Better-Than-Better to highlight the care we're putting into making the software easy to use, robust with no bugs, and elegant for the end-user.  The current state of research software is in poor condition as many of you know, and it's our philosophy that putting real effort into our first public release of this program will go a long ways in gaining adoption, perhaps even more so than a better algorithm (which we aim to also achieve).  



## Installation
This program requires: (links are provided if you don't already have it on your box or want the latest version)

- java (https://www.oracle.com/java/technologies/javase-downloads.html)
- python (https://www.python.org/downloads/)
- treetime (https://treetime.readthedocs.io/en/latest/installation.html)



## Executing The Program
**UPDATE:  Program is now clean & user-friendly (too many changes to detail, so give it run and see for yourself)**

From inside the pre-alpha/ directory, execute: `./poutine.sh --help`

As before for testing purposes, the reference set can be found in the reference_set/ directory.  A successful run of the reference set will end with the final line:  "CLEAN EXIT"

~~Also as before, you can `&>` the program output into a flat file, and then run consume_results.sh to see the top hits (currently working on redoing the entire output so expect a new commandline grouping of output options very soon).~~

Since we're in alpha mode, please try and break the program, including all sorts of commandline shenanigans! :sunglasses:

Code is now encouraging a new directory structure:

```
out_dir/
|___ ancestral_reconstruction_dir/
     |___ ancestral_sequences.fasta
     |___ annotated_tree.nexus
     |___ sequence_evolution_model.txt     
|___ poutine_session_current_time.log (any runtime errors will be found in this file)
|___ poutine_session_current_time.debug (temporary facility, feel free to use it as it gives a bunch of internal stats)
|___ poutine_session_current_time.out (main results file)
|___ poutine_session_current_time.out.sorted_by_a1_maxT (major allele associations across all sites sorted by maxT)
|___ poutine_session_current_time.out.sorted_by_a2_maxT (minor allele associations across all sites sorted by maxT)
```

For convenience, you can view the results output files in pretty-format using the following Unix command:
`column -ts "Ctrl-v <tab>" output_filename.out.sorted_by_a2_maxT | less -S`

To set TAB as the delimiter, macOS and Linux systems usually want `Ctrl-v <tab>` which is a Ctrl-v followed immediately by the tab key.

Since a picture is worth a thousands words:
![Image of new cli](https://github.com/Peter-Two-Point-O/Easy-Is-Better-Than-Better/blob/master/pre-alpha/images/poutine_cli_screenshot_2.png)

This is what a clean run should look like in the terminal with default settings:
![Image of new console updates](https://github.com/Peter-Two-Point-O/Easy-Is-Better-Than-Better/blob/master/pre-alpha/images/poutine_console_screenshot.png)

~~For this pre-alpha release, a convenience homoplasy_counter.sh file is included to run the program.  Future versions will likely be a precompiled binary with full command-line switch options like any other Unix binary (e.g. less, top, etc).~~

~~From inside the pre-alpha/ directory, simply execute: `./homoplasy_counter.sh &> output_filename.out &`~~

~~The `&>` redirects all terminal streams to the output file, while the final `&` backgrounds your process.  Both are optional.~~

~~To get a sorted list of top hits, simply run: `consume_results.sh output_filename.out 0.05`~~

~~`0.05` is the maxT p value cut-off user chooses.  This is a family-wise error rate so 5% or 1% are classically reasonable.  The script will also output the following files:~~

~~1)  `.all_associations`:  same as full homoplasy counter output file minus all the log information~~
~~2)  `.a1_sorted`:  each row represents one a1 allele in sorted fashion (a1 subset from .all_associations)~~
~~3)  `.a2_sorted`:  each row represents one a2 allele in sorted fashion (a2 subset from .all_associations)~~
~~4)  `.a1_top_hits_maxT`:  a1 associations filtered by maxT pvalue~~
~~5)  `.a2_top_hits_maxT`:  a2 associations filtered by maxT pvalue~~




## Reference Set
Included in this release is our *reference set* comprised of 124 mtb strains we have used for internal testing.  The homoplasy_counter.sh file starts a test run using this reference set.



## Current Major Bugs
No major bugs!

FIXED:  ~~- Concurrency bug during resampling:  This bug only occurs sometimes.  If your test run doesn't run to completion (`CLEAN EXIT` should be the last line in the output file) and you see the program hang, simply kill the process and restart it for now.~~

Sometimes the program would lock up and freeze in perpetuity.  Based upon my output showing aspects of the thread pool, it looked like a few of the threads sometimes did not execute till completion, and thus the countdown latch is never decremented.  I thought this was weird because no stacktrace was ever seen in the output.  Best guess was that we have a race condition at the heart of the bug.

I explicity attempted to catch all runtime errors inside each thread, and voila:  we get the classic NullPointerException during a run where the program locks up.  Turns out the bug is during read access of the hash table cache I used to optimize the running time (worked fine in single-threaded mode).  The cache is not thread-safe; specifically, the put() operation is not atomic in the sense that a key could be fully updated while the value is not fully updated, thus not satisfying the happens-before relation necessary for memory consistency between threads.  I'm guessing this memory inconsistency can happen either because the value actually hasn't been updated OR it has been updated but not flushed from the cpu memory (i.e. register or L1/L2 caches) into main memory for other threads to see, either way the value is null to other threads executing a get() on an identical key.

The fix is to make only part of the put() operation atomic so that all reads from the cache has a proper happens-before relation with the writes to the cache on an identical key.  No need to synchronize and lock larger sections of code, thus keeping threads more active.



## Recent Major Code Changes
- Removed q values from the significance assessment.  Thus, all things R have been removed.  Main reason for this feature removal is because the resampling-derived FWER (maxT variant) is sufficient for users to sort and look for top hits (remember the philosophy here is "easier is better than better").  It's also more robust than anything in the FDR space in the sense that any filtering of segregating sites will alter the p value distribution that FDR-based schemes use to calculate the splines, and thus can dramatically alter the final q values.  So for example, I have currently set the program to ignore all sites with <= 3 homoplasic mutations as these will not show statistical significance (less sites mean faster run with lower memory).  The maxT p values are typically immune to changes from site filtering.  In a future release, when we are likely to add estimation statistics like a resampling-derived effect size + confidence intervals, we can reconsider if q values will benefit the user in any way.
- Incorporated treetime for purposes of genotypic ancestral reconstruction using the default optimized joint probabilities method.
- The program should be fully platform independent now.  ~~Waiting for Windows users to get back to us to verify.  Tests on MacOS/Intel and Linux/Intel were successful.  NOTE:  Windows users should check out the code from the platform_independence branch.  Once we verify it works, we'll merge this branch back into master.~~
- Full command-line interface.  This will be our store-front! so will try to make this elegant and easy.
- Incorporate consume_results.sh code into main program along with other facilities to sort and pretty-format results.
- Organize all program output (e.g. various results files, log file, debugging file, proper console messages, etc).
- There are now a bounty of checks (e.g. malformed input file formats, etc) that allow the program to fail fast and meaningfully.  User is now prompted with console and log messages that provide either warnings and/or points them to further action.
- All output (whether to console or file) is now pretty formatted and human understandable, so if there is anything you want to see changed, just let us know!



## Upcoming Major Code Changes
- new commandline feature:  allow user to use a previously run ancestral reconstruction
- allow for missing phenotypes
- allow for missing genotypes
- report # monomorphic, bi/tri/quad-allelic sites instead of grouping them to all into non-biallelic category (in debug mode)
- Potentially incorporate vcf format (this will be another input option in addition to the current multi-fasta + plink map option that corresponds to the physical positions of the variants in the fasta file).  We could use feedback from folks on what your bacterial genomics workflow looks like so we can better accomodate users.




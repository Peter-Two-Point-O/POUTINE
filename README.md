# POUTINE

POUTINE is a homoplasy-counting-based method for genome-wide association studies (GWAS).  It is particularly well-suited to identify causal variants in strongly clonal populations.  For a deeper understanding of POUTINE and its use for microbial GWAS, please refer to our [manuscript](https://github.com/Peter-Two-Point-O/POUTINE/blob/master/README.md#citation).  

Stay tuned as POUTINE will continue to evolve to include [new major features](https://github.com/Peter-Two-Point-O/POUTINE/blob/master/README.md#upcoming-major-featurescode-changes) and enhancements.  POUTINE is designed to be user-friendly, and we welcome suggestions for improvements in our [Discussion section](https://github.com/Peter-Two-Point-O/POUTINE/discussions).  Collaborations are welcome!

## Updates

First blog post coming soon!

## Installation

POUTINE requires:

*   java (https://www.oracle.com/java/technologies/javase-downloads.html)
*   python (https://www.python.org/downloads/)
*   treetime (https://treetime.readthedocs.io/en/latest/installation.html)

Users are encouraged to use the latest version of the Java Virtual Machine (JVM) for both substantial speed and memory improvements compared to much older JVMs (e.g. Java 8) currently installed on many machines.

A Docker and Conda version of POUTINE will likely feature soon as well.

## Executing The Program

Execute: `./poutine.sh --help` to see all command-line options.

Setting poutine.sh into your global path will allow you to call poutine.sh from anywhere (you won't need the ./).

Since a picture is worth a thousand words:  

![Image of new cli](https://github.com/Peter-Two-Point-O/Easy-Is-Better-Than-Better/blob/master/images/poutine_cli_screenshot_3.png)

This is what a clean run should look like in the terminal with default settings:  

![Image of new console updates](https://github.com/Peter-Two-Point-O/Easy-Is-Better-Than-Better/blob/master/images/poutine_console_screenshot.png)

## Understanding The Output Files

POUTINE is now encouraging the following directory structure:

```
out_dir/
|___ ancestral_reconstruction_dir/
     |___ ancestral_sequences.fasta
     |___ ancestral_tree.newick
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

## Citation

Please cite . . . (preprint link coming shortly)

## License

See the [LICENSE](https://github.com/Peter-Two-Point-O/Easy-Is-Better-Than-Better/blob/5803828df31f41951f69a70b9cca61135562ff28/LICENSE) file for license rights and limitations (GNU GPL3).

## Upcoming Major Features/Code Changes

*   Set-testing to address allelic heterogeneity
*   Bootstrapped effect size calculations with confidence intervals
*   Ability to homoplasy count in the accessory genome
*   Ability to homoplasy count indels
*   Multi-thread homoplasy identification step

## Notes From Alpha Development Of POUTINE

### Recent Major Code Changes

*   Removed q-values from the significance assessment. Thus, all things R have been removed. The main reason for this feature removal is because the resampling-derived FWER (maxT variant) is sufficient for users to sort and look for top hits (remember the philosophy here is "easier is better than better"). It's also more robust than many methods in the FDR space because max(T) better address dependency structures between segregating sites.  In a future release, when we are likely to add estimation statistics like a resampling-derived effect size + confidence intervals, we can reconsider the progress of FDR-based methods for addressing dependence structures.
*   Incorporated treetime for purposes of genotypic ancestral reconstruction using the default optimized joint probabilities method.
*   The program should be fully platform-independent now. ~Waiting for Windows users to get back to us to verify. Tests on MacOS/Intel and Linux/Intel were successful. NOTE: Windows users should check out the code from the platform\_independence branch. Once we verify it works, we'll merge this branch back into master.~
*   Full command-line interface. This will be our store-front! so will try to make this elegant and easy.
*   Incorporate consume\_results.sh code into the main program along with other facilities to sort and pretty-format results.
*   Organize all program output (e.g. various results files, log file, debugging file, proper console messages, etc).
*   There is now a bounty of checks (e.g. malformed input file formats, etc) that allow the program to fail fast and meaningfully. The user is now prompted with console and log messages that provide either warnings and/or points them to further action.
*   All output (whether to console or file) is now pretty formatted and human understandable, so if there is anything you want to see changed, just let us know!
*   New command-line feature (-u or --use-precomputed-anc-recon): this new feature allows the user to utilize a precomputed ancestral reconstruction fasta and newick file. When the option is turned on, treetime is bypassed completely and the program understands that it is now looking at a fasta file with ancestral genotypes as well as a newick file with internal nodes labeled. In addition, I've updated the program to output the Newick tree that I parse out from treetime's nexus output file. This allows users to easily reuse this ancestral\_tree.newick file (along with its accompanying ancestral\_sequences.fasta file) for subsequent poutine sessions to bypass ancestral reconstruction and more quickly explore other GWAS settings and/or phenotypes. This new feature is of great utility for moderate to large datasets where the ancestral reconstruction phase is the rate-limiting step and there is a desire to run multiple poutine sessions based upon the same input genotypes and tree.
*   Program now identifies # monomorphic, bi/tri/quad-allelic sites, and reporting this metric in the log file (and debug file). This is particularly useful as a sanity check to see which sites are being assayed by poutine (only biallelic sites for now). I anticipate users not considering this point, and thus seeing the # non-biallelic sites can potentially alert the user to either problems in their dataset and/or # sites that are not being considered by poutine. This need came from analyzing the discovery set where a substantial number of sites turned out to be monomorphic (due likely to subclonal heterogeneity being picked up by the variant caller). A future version may include the capability of considering multiallelic sites perhaps using a multinomial test.
*   Program now works with missing phenotype data.

> Philosophy behind how the algorithm treats missing phenotype data:
> 
> Two ideas:
> 
> ```
> 1. Missing data is not the absence of data but instead the presence of information that encapsulates whether the missingness is random or non-random.   
> 2. The observed and null distributions must be comparable:  we can't have one dist as a f(case, control) and the other dist as a f(case, control, missing)
> ```

> Two ways missing phenotype data can be used with the program: Option 1 requires the user modify their input geno and pheno files. Option 2 requires nothing to be changed by the user.  
> \* Option 1: preprocess and delete all samples missing pheno data from both input genotype fasta file and input pheno file. Must be in both files!  
> \* Option 2: leave all input files untouched. The program designates missing pheno data with a non-0|1 value. This allows the program to treat missing/improperly coded pheno data as the presence of information and protect against potentially non-random bias in missingness. Samples with missing pheno data are kept in both input genotype file and input phenotype file. Doing this allows test statistics to be calculated with observed data that also includes missing phenos (eg at a specific segsite, a homoplasic mutation may be missing pheno data from the sample this genotype comes from, the observed test statistic will reflect this by not counting this homoplasic mutation towards either the tot # trials or tot # cases). During resampling, missing phenos are also permuted and test statistics also incorporate missing pheno data at homoplasic mutations in the same way as test statistics for observed data. This approach allows the observed distribution to be compareable to the resampled null distribution (they are both a function of 3 states: case, control, missing). What I do not recommend is the user keeping samples with missing pheno data in their input genotype file while removing these same samples from the input phenotype file, and vice versa. If this is done the two distributions, observed and null, are not as comparable because there could be non-random bias in missingness. Consider that if you use samples with missing pheno data in the input geno file, then sites containing observed homoplasic mutations with missing pheno data will be calculated with the missing data incorporated into the test statistic, while resampled test statistics will never see missingness because all samples with missing phenos were removed from the input pheno file! As such, the observed pvalues are a function of cases, controls, and non-random missingness in phenotypes and the resampled pvalues are only a function of cases and controls.

> Both options account for potential non-random missingness; any bias is reflected in both observed and null distributions. Option 1 does this by seeing missingness implicitly (ie all missingness is removed from all input files as a preprocessing step). Option 2 does this by explicitly modeling missingness (ie missing pheno data is incorporated in permutations and also in observed test statistics).

> Differences between the 2 options: I think option 2 is less error prone because I anticipate users removing samples with missing pheno data by simply deleting those rows from the pheno file only (and not also in the input geno file as required). I could check for this with code and prompt the user to modify input files, but then modifying fasta files by an average user could also introduce errors. With option 2, the program now protects a user if they are not even aware there are phenos with missing data/improperly coded.

*   Program now understands missing genotypes.

> The philosophy here is essentially the same as missing phenotypes (though the code is not!):

> I take a middle route where I preserve as much data as I can (i.e. segsites with some missing genotypes) while not adding to the data in any way (e.g. no imputation). In a similar fashion to how the program deals with missing phenos, any non-random structure in genotype missingness is reflected in both obs and null dists, as such they can be safely compared.

> In a future version for more explicit missing genotype control, one could build in preprocessing steps to both summarize genotype missingness by both sample and site, and also to test genotype missingness conditional on phenos at each site (i.e. are missing genotypes seen in more cases than controls?) say using a fisher's exact test. This allows the user to preprocess and remove any sites where the genotype missingness looks highly non-random. This preprocessing approach is how plink handles missing geno data.

Suite of "niceties" to make the program easier to use:

1.  Program now checks for mismatched sample names across input files. This is a common error and thus worth checking or else results may be incorrect. Specifically, the program now checks for mismatched sample names across input genotype, phenotype, and tree files. Note that this program itself does not place any restrictions on the makeup of the string that comprises the sample name. The only requirement is that sample names across files match. Program fails-fast if the sample names do not match, and the log file now reports the sample names in the phenotype file that do not match with either genotype and/or tree file.
2.  Improvements to log file, command-line help message, and terminal progress indicator.

### Current Major Bugs

**No major bugs!**

**FIXED:**  Concurrency bug during resampling: This bug only occurs sometimes. If your test run doesn't run to completion (`CLEAN EXIT` should be the last line in the output file) and you see the program hang, simply kill the process and restart it for now.

Sometimes the program would lock up and freeze in perpetuity. Based upon my output showing aspects of the thread pool, it looked like a few of the threads sometimes did not execute till completion, and thus the countdown latch is never decremented. I thought this was weird because no stack trace was ever seen in the output. Best guess was that we have a race condition at the heart of the bug.

I explicitly attempted to catch all runtime errors inside each thread, and voila: we get the classic NullPointerException during a run where the program locks up. Turns out the bug is during read access of the hash table cache I used to optimize the running time (worked fine in single-threaded mode). The cache is not thread-safe; specifically, the put() operation is not atomic in the sense that a key could be fully updated while the value is not fully updated, thus not satisfying the happens-before relation necessary for memory consistency between threads. I'm guessing this memory inconsistency can happen either because the value actually hasn't been updated OR it has been updated but not flushed from the cpu memory (i.e. register or L1/L2 caches) into main memory for other threads to see, either way the value is null to other threads executing a get() on an identical key.

The fix is to make only part of the put() operation atomic so that all reads from the cache have a proper happens-before relation with the writes to the cache on an identical key. No need to synchronize and lock larger sections of code, thus keeping threads more active.

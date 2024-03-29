/**
    POUTINE:  a homoplasy counter for genome-wide association studies.
    Copyright (C) 2020  Peter E. Chen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.gersteinlab.coevolution.core.data.*;
import org.gersteinlab.coevolution.core.io.*;
import picocli.CommandLine;
import picocli.CommandLine.*;
import picocli.CommandLine.Help.*;
import picocli.CommandLine.Model.*;

//import static org.fusesource.jansi.AnsiConsole.*;
//import static org.fusesource.jansi.Ansi.*;

import java.io.*;
import java.nio.file.Files;
import java.time.Duration;
import java.time.Instant;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import static java.nio.file.StandardOpenOption.APPEND;


/**
 * @author Peter E Chen
 * @version 1.0.0
 */
@Command(name = "poutine", version = "%nPOUTINE 1.0.0%n", mixinStandardHelpOptions = true, usageHelpWidth = 210, sortOptions = false, headerHeading = "%n", optionListHeading = "%n", footerHeading = "%n")
public class Homoplasy_Counter implements Callable<Integer> {
    @Spec
    static CommandSpec spec;

    // current commandline options:

    // TODO:  change all String filenames to File?
    @ArgGroup(exclusive = false, multiplicity = "1", validate = true)
    InputFiles inputFiles;
    static class InputFiles {
        @ArgGroup(exclusive = true, multiplicity = "1", validate = true, heading = "@|bg(213) %nInput genotype file: vcf or fasta format%n|@")
        Genotypes genotypes;
        static class Genotypes {
            @Option(names = {"-v", "--vcf"}, description = "Multi-sample vcf file (this option is currently not in use/available)")
            private String vcf_filename;

            @Option(names = {"-f", "--fasta"}, description = "Multi-fasta file (ancestral fasta file if -u in use). Variable sites only. Sites must be in the same order as in the physical position file.")
            private String msa_fasta_filename;
        }

        @ArgGroup(exclusive = false, multiplicity = "1", validate = true, heading = "@|bg(123) %nOther input files%n|@")
        OtherInputFiles otherInputFiles;
        static class OtherInputFiles {
            @Option(names = {"-p", "--phenos"}, description = "Phenotype file (case = 1, control = 0. Format is two columns tab-delimited with no header: column one is sample id, column two is phenotype)", required = true)
            private String pheno_filename;

            @Option(names = {"-t", "--tree"}, description = "Newick file (ancestral newick file if -u in use)", required = true)
            private String newick_filename;

            @Option(names = {"-m", "--map"}, description = "Physical positions file (PLINK map format or tab-delimited with no header: last column represents physical positions, all other columns not used)", required = true)
            private String map_filename;
        }
    }

    @ArgGroup(validate = false, heading = "@|bg(222) %nAlgorithm parameters%n|@")
    AlgoParams algoParams;
    static class AlgoParams {
        private static int m = 100000; // default value
        @Option(names = {"-r", "--replicates"}, paramLabel = "<# replicates>", description = "# Replicates used for resampling (default: 100000)", required = false, order = 1)
        private void validate_and_set_replicates_option(int user_value) {
            if (user_value > 0)
                m = user_value;
            else
                throw new ParameterException(spec.commandLine(), String.format("%nInvalid value '%s' for option '--replicates': " + "# replicates used for resampling must be > 0", user_value));
        }

        private static int min_hcount = 1;  // default value.  1 := use all homoplasic seg sites, 0 := use all seg sites including those without any homoplasic mutations on either allele
        @Option(names = {"-c", "--min_hcount"}, paramLabel = "<count>", description = "Minimum # homoplasic mutations required at each allele (default: 1)", required = false, order = 2)
        private void validate_and_set_min_hcount_option(int user_value) {
            if (user_value >= 0)
                min_hcount = user_value;
            else
                throw new ParameterException(spec.commandLine(), String.format("%nInvalid value '%s' for option '--min_hcounts': " + "# homoplasic mutations at a segregating site must be >= 0", user_value));
        }

        @Option(names = {"-u", "--use-precomputed-anc-recon"}, description = "Use precomputed ancestral reconstruction as input. Specify --fasta for ancestral fasta file and --tree for ancestral newick file.", required = false, order = 3)
        private static boolean use_precomputed_anc_recon = false;
    }

    @ArgGroup(validate = false, heading = "@|bg(85) %nRuntime options%n|@")
    RuntimeSettings runtimeSettings;
    static class RuntimeSettings {
        private static int num_threads = Runtime.getRuntime().availableProcessors();  // default value := max # logical processors available (# cpu cores and hyperthreading)
        @Option(names = {"-T", "--threads"}, paramLabel = "<# threads>", description = "# Threads used for parallel processing (default: max # logical processors)", required = false)
        private void validate_and_set_thread_option(int user_value) {
            if (user_value > 0)
                num_threads = user_value;
            else
                throw new ParameterException(spec.commandLine(), String.format("%nInvalid value '%s' for option '--threads': " + "# threads must be > 0", user_value));
        }
    }

    @ArgGroup(validate = false, heading = "@|bg(197) %nOutput files and options%n|@")
    OutputOptions outputOptions = new OutputOptions();
    static class OutputOptions {
        @Option(names = {"-d", "--out-dir"}, description = "Output directory for all output files (default: poutine_session_current_time)", required = false)
        private File output_dir;  // default is relative to the current dir

        @Option(names = {"-o", "--out"}, description = "Output file (relative to --out-dir) containing association results (default: out-dir/poutine_session_current_time.out)", required = false)
        private File output_file;
        private File output_file_sorted_by_a1_maxT;
        private File output_file_sorted_by_a2_maxT;
        private BufferedWriter out;
        private BufferedWriter out_sorted_by_a1_maxT;
        private BufferedWriter out_sorted_by_a2_maxT;

        @Option(names = {"-l", "--log"}, description = "Log file (relative to --out-dir) containing this session's log ouput and any runtime errors that may occur (default: out-dir/poutine_session_current_time.log)", required = false)
        private File log_file;
        private BufferedWriter log;

        @Option(names = {"-a", "--anc-recon-out-dir"}, description = "Output directory (relative to --out-dir) for all files related to ancestral reconstruction (default: out-dir/ancestral_reconstruction_current_time)", required = false)
        private File anc_recon_dir;  // default is a subdirectory of output_dirname

        @Option(names = {"-s", "--timestamp"}, description = "Add current session's start timestamp to user-specified dir and filenames (default: no timestamp)", required = false)
        private boolean timestamp = false;

        @Option(names = {"-X", "--force-overwrite"}, description = "Allow overwriting of existing files having the same user-specified filenames (default: program safely exits if existing files detected)", required = false)
        private boolean force_overwrite = false;

        @Option(names = {"-D", "--debug"}, description = "Turns debug mode on: output all debugging information to debug file out-dir/poutine_session_current_time.debug (this option is likely temporary for testing purposes)", required = false)
        private boolean DEBUG_MODE = false;
        private BufferedWriter debug;
    }

// Remaining global variables:
    private final String VERSION = "1.0.0";
//    private final String newick_filename;
//    private final String nexus_filename;  // treetime's annotated nexus file containing tree with internal nodes labeled
//    private final String ancestral_reconstruction_filename;  // treetime's fasta file containing ancestral genotypes
//    private final String pos_cross_ref_filename;  // cfml position file no longer in use

//    private final String pheno_filename;
//    private final String map_filename;
//    private final String msa_fasta_filename;  // user-specified msa used to infer input phylogeny
//    private final String user_input_newick_filename;  // user-specified phylogeny

//    private final String recomb_tracts_filename;
    private final boolean qvalues_option = false;  // false := turn off qvalue calculations as well as outputting these results

//    private final String R_dir = "/Users/blame_monster/Research/Indie/Homoplasy_Counting/R_code/";  // dir that contains R code (hard-coded for local)
    private final String R_dir = "/home/blame_monster/Research/Homoplasy_Counting/R_Code/";  // hard-coded for Kimura for testing
    private final boolean EXTANT_NODES_ONLY = true;  // turns off ancestrally reconstructed phenotypes, ie only use leaves
//    private final String anc_recon_dir = "ancestral_reconstruction";  // path is relative to user's current working dir

//    private final double p_success = 0.56;  // mtb genome-wide homoplasic mutations only (derived from diagnostic code)
//    private final double p_success = 0.38;  // mtb genome-wide all mutations 47 cases /(47 + 77 controls)
//    private final double p_success = 0.68;    // breo genome-wide homoplasic mutations only (derived from diagnostic code)
//    private final double p_success_a1 = 0.58622974;  // breo homoplasic mutations a1 only (derived from diagnostic code)
//    private final double p_success_a2 = 0.7205308;  // breo homoplasic mutations a2 only (derived from diagnostic code)
//    private final double p_success_a1 = 0.5625;  // mtb homoplasic mutations a1 only (derived from diagnostic code)
//    private final double p_success_a2 = 0.5654608;  // mtb homoplasic mutations a2 only (derived from diagnostic code)
//    private final double p_success_a1 = 0.379;  // mtb all mutations resampling space (p_success|both alleles)
//    private final double p_success_a2 = 0.379;  // mtb all mutations resampling space (p_success|both alleles)

    private enum Allele_Type {A1, A2}
    private enum Resampling_Space_Type {ALL_MUTATIONS, HOMOPLASIC_MUTATIONS_ONLY}
    private final Resampling_Space_Type R_SPACE_TYPE = Resampling_Space_Type.ALL_MUTATIONS;
    private final AlternativeHypothesis TAIL_TYPE = AlternativeHypothesis.GREATER_THAN;  // Apache Commons enum {GREATER_THAN := right-sided test (cases), TWO_SIDED, LESS_THAN}
//    private final AlternativeHypothesis TAIL_TYPE = AlternativeHypothesis.TWO_SIDED;
//    private final AlternativeHypothesis TAIL_TYPE = AlternativeHypothesis.LESS_THAN;

//    private final int m = 1000000;  // # replicates (# perms)
    private int m;  // # replicates (# perms)
//    private final int min_hcount = 7;  // 1 := use all homoplasic seg sites, 0 := use all seg sites including those without any homoplasic mutations on either allele

    private boolean DEBUG_MODE;

    private final ZonedDateTime session_start_time = ZonedDateTime.now();

    private int tot_num_sample_cases;
    private int tot_num_sample_controls;

    private HashMap<String, Integer> binom_test_calls = new HashMap();
    private int tot_num_calls_by_int_counter = 0;
    private Duration tot_duration_of_all_binom_calls_to_math_library = Duration.ZERO;


    public static void main(String[] args) {
        CommandLine cmdline = new CommandLine(new Homoplasy_Counter());
        cmdline.setUsageHelpLongOptionsMaxWidth(40);

        // registering a custom commandline exception handler to get rid of all ANSI escape sequences when outputting uncaught runtime exception; this is useful because program redirects these types of exceptions to file.
//        IParameterExceptionHandler exception_handler = cmdline.getParameterExceptionHandler();
//        cmdline.setParameterExceptionHandler((param_exception, cmdline_args) ->
//                {
//                    System.setProperty("picocli.ansi", "false");
//                    return exception_handler.handleParseException(param_exception, cmdline_args);
//                }
//        );
//        ColorScheme color_scheme = cmdline.getColorScheme();
//        System.out.println("color_scheme = " + color_scheme);

        int exitCode = cmdline.execute(args);
        System.exit(exitCode);
    }


    @Override
    public Integer call() throws Exception {

        more_cmdline_magic();
        log_cmdline_global_vars();
        String input_newick = read_and_validate_input_newick_format();

        System.out.printf(Ansi.AUTO.string("@|fg(213) %nStarting poutine session: " + session_start_time.format(DateTimeFormatter.ofPattern("YYYY-LLLL-dd EEEE HH'h':mm'm':ss's' O")) + "%n|@"));

        // WARNING:  warn user if user-specified # replicates is potentially too low to properly estimate statistical significance
        if (m < 10000) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nWarning: # replicates has been set to a potentially low value of %d. Please consider a higher value (the default is 100000 replicates).%n|@"), m);
            outputOptions.log.write(String.format("%nWarning: # replicates has been set to a potentially low value of " + m + ".  Please consider a higher value (the default is 100000 replicates).%n"));
            outputOptions.log.flush();
        }

        /*
         * ancestral recon pseudocode:
         * - kick off process to start treetime <--- ancestral_reconstruction(). this should be it's own method without any other major tasks in it
         *
         * These methods should remain in main() as it is currently:
         * - read in output file nexus          <--- nexus_to_newick()
         * - read in output file msa fasta      <--- build_seg_sites()
         * - build_tree(newick)
         */
        String ancestral_newick = null;
        String ancestral_fasta = null;
        if (!algoParams.use_precomputed_anc_recon) {
            ancestral_reconstruction();

            // parse treetime's annotated nexus file (with internal node's labelled) -> newick.  In addition, set newick to point to this ancestral newick.
            ancestral_newick = nexus_to_newick();

            // write newick to file for potential subsequent use
            write_newick(ancestral_newick);

            // set ancestral_fasta to point to treetime's ancestral_sequences.fasta file
            ancestral_fasta = outputOptions.anc_recon_dir + File.separator + "ancestral_sequences.fasta";
        }  else {
            System.out.printf(Ansi.AUTO.string("@|fg(123) %nSkipping ancestral reconstruction . . .%n|@"));

            // set ancestral_newick to point to user-specified precomputed ancestral newick file
            ancestral_newick = input_newick;

            // set ancestral_fasta to point to user-specified precomputed ancestral fasta file
            ancestral_fasta = inputFiles.genotypes.msa_fasta_filename;
        }

        System.out.printf(Ansi.AUTO.string("@|fg(85) %nConstructing internal data structures for tree, segregating sites, phenotypes . . . |@"));
        // create tree data structure
        NewickTree tree = build_tree(ancestral_newick);

        // read in ancestral genotypes
        HashMap<String, char[]> seg_sites = build_seg_sites(ancestral_fasta);

        // PLINK map file that contains physical positions in the order of the ped file used to create input files for CFML
        int[] physical_poss = get_physical_positions();

        // TODO:  refactor all ArrayLists into arrays[] since I know the sizes ahead of time?

        // check for mismatch of # seg sites between ancestral genotype file and input positions file
        check_num_seg_sites(seg_sites, physical_poss);

        HashMap<String, String> phenos = read_phenos();  // only extant phenotypes

        // check for mismatched sample names across files:  input geno, pheno, tree files.  Fail fast if sample names do not completely match.
        check_sample_names(seg_sites, phenos, tree);
        System.out.printf("structures complete.%n");

        // homoplasy counts
        ArrayList<Homoplasy_Events> all_events = count_all_homoplasy_events(tree, seg_sites, physical_poss);

        /** closer look at site 1674048
//        System.out.println("all_events.size() = " + all_events.size());
//        all_events.forEach(s -> System.out.println(s.toString()));
//        System.out.println("\n\n");
        int segsite_index = -1;
        for (int i = 0; i < all_events.size(); i++) {
            if (all_events.get(i).physical_pos == 1674048)
                segsite_index = i;
        }
        this.output_segsite_events_as_branch_colors(all_events, segsite_index);
        System.exit(-1);
        */

        assoc(all_events, phenos);

        end_session();

        return 0;
    }


    /**
     *
     * @param seg_sites
     * @param physical_poss
     */
    private void check_num_seg_sites(HashMap<String,char[]> seg_sites, int[] physical_poss) {
        int num_sites_from_anc_recon_file = seg_sites.values().stream().findFirst().get().length;
        int num_sites_from_map_file = physical_poss.length;

        if (num_sites_from_anc_recon_file != num_sites_from_map_file) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while checking # segregating sites between ancestral fasta file and physical positions file.  Please see the log file for more information on this error.%n|@"));

            try {
                outputOptions.log.write(String.format("%nError: the # segregating sites in the ancestral fasta file does not equal the # of segregating sites in the input physical positions file.%n" +
                                                      "The input fasta file must contain only variable sites and must be in the same order as the input physical positions file.%n"));
                outputOptions.log.write(String.format("%n# segregating sites detected in ancestral fasta file = " + num_sites_from_anc_recon_file +
                                                      "%n# segregating sites detected in physical positions file = " + num_sites_from_map_file + "%n"));
                outputOptions.log.flush();
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the log file. Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }

            System.exit(-1);
        }
    }



    /**
     * This method checks for mismatched sample names across input geno, pheno, and tree files.
     * Note that sample size is based upon the samples seen in the pheno file.
     * Also note that this program itself does not place any restrictions on the makeup of the
     * string that comprises the sample name.  Only requirement is that sample names across
     * files match.
     *
     * Method fails-fast if the sample names do not match.  Log file reports
     * the sample names in the pheno file that do not match with either genotype
     * and/or tree file.
     *
     * @param seg_sites
     * @param phenos
     * @param tree
     */
    private void check_sample_names(HashMap<String,char[]> seg_sites, HashMap<String, String> phenos, NewickTree tree) {
        ArrayList<String> bad_sample_names = new ArrayList<>();

        List<NewickTreeNode> leaves = tree.getLeafNodes();
        HashSet<String> leaf_names = get_leaf_names(leaves);

        Iterator<Map.Entry<String, String>> it = phenos.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<String, String> pheno = it.next();
            String pheno_sample_name = pheno.getKey();

            boolean matches_against_genofile = seg_sites.containsKey(pheno_sample_name);
            boolean matches_against_treefile = leaf_names.contains(pheno_sample_name);
            if (!matches_against_genofile || !matches_against_treefile)
                bad_sample_names.add(pheno_sample_name);
        }

        if (bad_sample_names.size() > 0) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while checking sample names across input genotype, phenotype, and tree files.  Please see the log file for more information on this error.%n|@"));

            try {
                outputOptions.log.write(String.format("%nError: the following sample names from the input phenotype file do not match with sample names in either/both the input genotype and/or tree files:%n"));

                for (String bad_sample_name : bad_sample_names) {
                    outputOptions.log.write(bad_sample_name);
                    outputOptions.log.newLine();
                }

                outputOptions.log.flush();
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the log file. Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }

            System.exit(-1);
        }
    }




//    private String read_newick() {
//
//        String newick = null;
//
//        try (BufferedReader br = new BufferedReader(new FileReader(inputFiles.otherInputFiles.newick_filename))) {
//
//        } catch (IOException e) {
//            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while reading the ancestral newick tree file. Please see the stack trace in the log file for more information on this error.%n|@"));
//            e.printStackTrace();
//            System.exit(-1);
//        }
//
//        return newick;
//    }



    /**
     * This method writes out a newick tree to a file named ancestral_tree.newick (inside anc_recon_dir) so that users have the ability to reuse this
     * newick file for subsequent poutine sessions without invoking a redundant ancestral reconstruction.  See --use-precomputed-anc-recon cmdline option.
     *
     * @param newick
     */
    private void write_newick(String newick) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outputOptions.anc_recon_dir + File.separator + "ancestral_tree.newick"))) {
            bw.write(newick);
            bw.newLine();
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing the ancestral newick tree to the ancestral reconstruction directory. Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }
    }



    /**
     *  This method checks the newick format of the input newick tree file.
     *
     *  If --use-precomputed-anc-recon is turned on, this input newick file represents a newick with ancestral nodes labelled
     */
    private String read_and_validate_input_newick_format() {

        String newick = null;
        try (BufferedReader br = new BufferedReader(new FileReader(inputFiles.otherInputFiles.newick_filename))) {
            StringBuilder sb = new StringBuilder();
            String currline;
            while ((currline = br.readLine()) != null) {
                sb.append(currline);
            }
            newick = sb.toString();
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while reading the input newick tree file.  Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        try {
            NewickTree tree = null;
            StringReader sr = new StringReader(newick);
            NewickTreeTokenizer tokenizer = new NewickTreeTokenizer(sr);
            NewickTreeReader tree_reader = new NewickTreeReader(tokenizer);
            tree = tree_reader.readTree();
        } catch (DataFormatException | IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while parsing the input newick tree (this is likely due to a malformed newick file).  Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        return newick;
    }



    /*
     * At some point, all commandline validation checks may go here to keep the code tidy.
     */
    private void more_cmdline_magic() {
        m = AlgoParams.m;  // not the most elegant (keep for now), but i want to keep the m variable the same as before (before exposing the parameter on the commandline) instead of refactoring to AlgoParams.m all over the place.
        DEBUG_MODE = outputOptions.DEBUG_MODE;

        // input files: ==========================================
        if ((inputFiles.genotypes.msa_fasta_filename != null) && (!new File(inputFiles.genotypes.msa_fasta_filename).exists()))  // user has selected --fasta option AND user-specified fasta file does not exist
            throw new ParameterException(spec.commandLine(), String.format("%nFile '%s' for option '--fasta' cannot be found.  Please make sure the correct path and file name are specified.", inputFiles.genotypes.msa_fasta_filename));

        if ((inputFiles.genotypes.vcf_filename != null) && (!new File(inputFiles.genotypes.vcf_filename).exists()))  // user has selected --vcf option AND user-specified vcf file does not exist
            throw new ParameterException(spec.commandLine(), String.format("%nFile '%s' for option '--vcf' cannot be found.  Please make sure the correct path and file name are specified.", inputFiles.genotypes.vcf_filename));

        if (!new File(inputFiles.otherInputFiles.pheno_filename).exists())
            throw new ParameterException(spec.commandLine(), String.format("%nFile '%s' for option '--phenos' cannot be found.  Please make sure the correct path and file name are specified.", inputFiles.otherInputFiles.pheno_filename));

        if (!new File(inputFiles.otherInputFiles.newick_filename).exists())
            throw new ParameterException(spec.commandLine(), String.format("%nFile '%s' for option '--tree' cannot be found.  Please make sure the correct path and file name are specified.", inputFiles.otherInputFiles.newick_filename));

        if (!new File(inputFiles.otherInputFiles.map_filename).exists())
            throw new ParameterException(spec.commandLine(), String.format("%nFile '%s' for option '--map' cannot be found.  Please make sure the correct path and file name are specified.", inputFiles.otherInputFiles.map_filename));

        // output:
        final String basename = "poutine_session_";
        final String curr_time = session_start_time.format(DateTimeFormatter.ofPattern("YYYY-MM-dd_HH'h'-mm'm'-ss's'"));

        // output dirs: ===========================================
        outputOptions.output_dir = get_validated_toplevel_path(outputOptions.output_dir, basename + curr_time, curr_time);
        outputOptions.output_dir.mkdir();  // create dir now that it has passed all checks (eg: dir does not exist yet or --force-overwrite is specified by user)

        if (!algoParams.use_precomputed_anc_recon)
            outputOptions.anc_recon_dir = get_validated_relative_path(outputOptions.output_dir, outputOptions.anc_recon_dir, outputOptions.output_dir + File.separator + "ancestral_reconstruction_" + curr_time, curr_time);

        // output files: ===========================================
        outputOptions.output_file = get_validated_relative_path(outputOptions.output_dir, outputOptions.output_file, outputOptions.output_dir + File.separator + basename + curr_time + ".out", curr_time);
        outputOptions.output_file_sorted_by_a2_maxT = new File(outputOptions.output_file.getAbsolutePath() + ".sorted_by_a2_maxT");
        outputOptions.output_file_sorted_by_a1_maxT = new File(outputOptions.output_file.getAbsolutePath() + ".sorted_by_a1_maxT");
        try {
            outputOptions.out = new BufferedWriter(new FileWriter(outputOptions.output_file));
            outputOptions.out_sorted_by_a2_maxT = new BufferedWriter(new FileWriter(outputOptions.output_file_sorted_by_a2_maxT));
            outputOptions.out_sorted_by_a1_maxT = new BufferedWriter(new FileWriter(outputOptions.output_file_sorted_by_a1_maxT));
        } catch (IOException ioe) {
            throw new ParameterException(spec.commandLine(), String.format("%nFile '%s' for option '--out' is not a valid file.  Please enter a valid filename.", outputOptions.output_file.toString()));
        }

        outputOptions.log_file = get_validated_relative_path(outputOptions.output_dir, outputOptions.log_file, outputOptions.output_dir + File.separator + basename + curr_time + ".log", curr_time);
        try {
            outputOptions.log = new BufferedWriter(new FileWriter(outputOptions.log_file));

            // set all stack traces to be sent to log file
            // properly appends to file and autoflushes to make sure in the event of a runtime error all info is written out before program execution ends
            System.setErr(new PrintStream(Files.newOutputStream(outputOptions.log_file.toPath(), APPEND), true));
        } catch (IOException e) {
            throw new ParameterException(spec.commandLine(), String.format("%nFile '%s' for option '--log' is not a valid file.  Please enter a valid filename.", outputOptions.log_file.toString()));
        }

        if (DEBUG_MODE) {
            try {
                outputOptions.debug = new BufferedWriter(new FileWriter(outputOptions.output_dir + File.separator + basename + curr_time + ".debug"));
                outputOptions.debug.write("poutine " + VERSION + " session start time: " + session_start_time.format(DateTimeFormatter.ofPattern("YYYY-LLLL-dd EEEE HH'h':mm'm':ss's' O")));
                outputOptions.debug.newLine();
                outputOptions.debug.flush();  // since in DEBUG mode, I want to capture everything in the buffer immediately in case of early program termination
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the debug file. Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }
    }



    /**
     * This method is used to set a directory path.  Top-level meaning that this allows the user to set a dir anywhere on their filesystem.
     * This method also checks the following:
     *      - is the default pathname needed?  (ie poutine_session_current_time)
     *      - is --timestamp option in use?
     *      - is --force-overwrite option in use?
     *
     * @param dir
     * @param default_pathname
     * @param curr_time
     * @return dir is a top-level output dir either specified by the user on the commandline or is set to the default:  poutine_session_current_time (relative to user's current working directory)
     */
    private File get_validated_toplevel_path(File dir, String default_pathname, String curr_time) {

        // set default if no user-specified output file
        if (dir == null) // if no user-specified output file, set to default file
            dir = new File(default_pathname);
        else if (outputOptions.timestamp) // else if user has specified output file AND specified option to add timestamp to this filename
            dir = new File(dir.getPath() + "." + curr_time);

        // check if the file or dir exists, clean exit if so (unless --force-overwrite is specified by user)
        if (dir.exists() && !outputOptions.force_overwrite)
            throw new ParameterException(spec.commandLine(), String.format("%nDirectory '%s' already exists.  Use --force-overwrite option to use existing directory.", dir.toString()));

        return dir;
    }


    /**
     * This method is used to set a relative path.  Relative meaning that the returned path is always relative to the top-level output dir.  This allows the program to enforce a directory structure to contain
     * all output files.
     * This method also checks the following:
     *      - is the default pathname needed?  (ie poutine_session_current_time)
     *      - is --timestamp option in use?
     *      - is --force-overwrite option in use?
     *
     * @param output_dir
     * @param file
     * @param default_pathname
     * @param curr_time
     * @return file is a file relative to the top-level output dir
     */
    private File get_validated_relative_path(File output_dir, File file, String default_pathname, String curr_time) {

//        if (output_dir == null)

        // set default if no user-specified output file
        if (file == null) // if no user-specified output file, set to default file
            file = new File(default_pathname);
        else if (outputOptions.timestamp) // else if user has specified output file AND specified option to add timestamp to this filename
            file = new File(output_dir + File.separator + file.getName() + "." + curr_time);
        else
            file = new File(output_dir + File.separator + file.getName());  // else simply place user's custom file inside output dir

        // check if the file or dir exists, clean exit if so (unless --force-overwrite is specified by user).
        if (file.exists() && !outputOptions.force_overwrite)
            throw new ParameterException(spec.commandLine(), String.format("%nFile '%s' already exists.  Use --force-overwrite option to overwrite existing file.", file.toString()));

        return file;
    }



    /*
     * This method kicks off the program log by logging all command line options (including default values used).
     */
    private void log_cmdline_global_vars() {

        BufferedWriter log = outputOptions.log;

        try {
            log.write("poutine " + VERSION + " session start time: " + session_start_time.format(DateTimeFormatter.ofPattern("YYYY-LLLL-dd EEEE HH'h':mm'm':ss's' O")));
            log.newLine();
            log.newLine();

            log.write("Command line options used:");
            log.newLine();

            // input files:
            // genotypes
            log.write("genotype file: " + inputFiles.genotypes.msa_fasta_filename);
            log.newLine();

            // phenos
            log.write("phenotype file: " + inputFiles.otherInputFiles.pheno_filename);
            log.newLine();

            // tree
            log.write("tree file: " + inputFiles.otherInputFiles.newick_filename);
            log.newLine();

            // plink map file
            log.write("plink map file: " + inputFiles.otherInputFiles.map_filename);
            log.newLine();

            // # replicates/permutations (m)
            log.write("# replicates: " + AlgoParams.m);
            log.newLine();

            // min_hcount
            log.write("min # homoplasic mutations required at each segregating site: " + AlgoParams.min_hcount);
            log.newLine();

            // use-precomputed-anc-recon
            log.write("Is --use-precomputed-anc-recon in use?: " + algoParams.use_precomputed_anc_recon);
            log.newLine();

            // # threads
            log.write("# threads used: " + RuntimeSettings.num_threads);
            log.newLine();

            // output dir
            log.write("ouput directory: " + outputOptions.output_dir.getAbsolutePath());
            log.newLine();

            // out files
//            log.write(String.format("output files:%n" +
            log.write(String.format(
                    "     log file: " + outputOptions.log_file.getAbsolutePath() + "%n" +
                    "     results file (sorted by physical position): " + outputOptions.output_file.getAbsolutePath() + "%n" +
                    "     results file by major allele (sorted by maxT pvalue): " + outputOptions.output_file_sorted_by_a1_maxT.getAbsolutePath() + "%n" +
                    "     results file by minor allele (sorted by maxT pvalue): " + outputOptions.output_file_sorted_by_a2_maxT.getAbsolutePath()));
            log.newLine();

            if (!algoParams.use_precomputed_anc_recon) {
                // anc recon output dir
                log.write("ancestral reconstruction output directory: " + outputOptions.anc_recon_dir.getAbsolutePath());
                log.newLine();

                // anc recon output files
                log.write(String.format(
                        "     ancestral genotypes file: " + outputOptions.anc_recon_dir.getAbsolutePath() + File.separator + "ancestral_sequences.fasta%n" +
                        "     ancestral tree file: " + outputOptions.anc_recon_dir.getAbsolutePath() + File.separator + "ancestral_tree.newick"));
                log.newLine();
            } else {
                log.write("ancestral reconstruction output directory: not in use because --use-precomputed-anc-recon option turned on.");
                log.newLine();
            }

            // timestamp
            log.write("Is --timestamp option in use?: " + outputOptions.timestamp);
            log.newLine();

            // force-overwrite
            log.write("Is --force-overwrite option in use?: " + outputOptions.force_overwrite);
            log.newLine();

            // debug mode
            log.write("DEBUG MODE: " + outputOptions.DEBUG_MODE);
            log.newLine();

            // flush buffer to capture this info in case of early program termination
            log.flush();
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the log file. Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);  // TODO:  settle on a set of exit codes and document this
        }
    }


//  NOTE:  moved everything in main() to call()
//    public static void main(String[] args) {
//
//        Homoplasy_Counter homoplaser = new Homoplasy_Counter(args);
//
//        // PLINK .map file that contains physical positions in the order of the ped file used to create input files for CFML
//        int[] physical_poss = homoplaser.get_physical_positions();
//
//        /*  cfml version
//        // phylogeny reflecting ALL seg sites
//        NewickTree tree = homoplaser.build_tree();
//        */
//
//        /*
//         * ancestral recon pseudocode:
//         * - kick off process to start treetime <--- ancestral_reconstruction(). this should be it's own method without any other major tasks in it
//         *
//         * These methods should remain in main() as it is currently:
//         * - read in output file nexus          <--- nexus_to_newick()
//         * - read in output file msa fasta      <--- build_seg_sites()
//         * - build_tree(newick)
//         */
//        homoplaser.ancestral_reconstruction();
//
//        // parse treetime's annotated nexus file (with internal node's labelled) -> newick
//        String newick = homoplaser.nexus_to_newick();
//
//        // create tree data structure
//        NewickTree tree = homoplaser.build_tree(newick);
//
///*      // process cfml version of genotype ancestral reconstruction:
//        HashMap<String , char[]> seg_sites = homoplaser.build_seg_sites();
//        int[] segsite_indexes = homoplaser.get_segsite_indexes();
//        // sanity check for physical positions (i'm moving between map/ped <-> CFML position_cross_reference/ML_sequence.fasta worlds)
//        System.out.println("physical_poss.length = " + physical_poss.length);
//        System.out.println("segsite_indexes.length = " + segsite_indexes.length);  // this should = above length
////        System.out.println("seg_sites.size() = " + seg_sites.size());  // hmmm this returns # nodes (extant + ancestral) not # of alleles at each site
//*/
//
//        // process treetime's ancestral reconstruction fasta output file:
//        // read in ancestral genotypes
//        HashMap<String , char[]> seg_sites = homoplaser.build_seg_sites();
//        // output # segregating sites (map file and anc recon fasta file)
//        System.out.println("physical_poss.length = " + physical_poss.length);
////        seg_sites.values().iterator().next().length  // # segsites from treetime's ancestral fasta (in first HashMap entry only, should be same as rest of the entries)
//        // TODO:  rename build_seg_sites() -> read_ancestral_seqs()?  rename seg_sites -> ancestral_seqs?
//
//
//        // TODO:  refactor all ArrayLists into arrays[] since I know the sizes ahead of time?
//
//        HashMap<String, String> phenos = homoplaser.read_phenos();  // only extant phenotypes
//
//        homoplaser.log_session_info();
//
//        // homoplasy counts
////        ArrayList<Homoplasy_Events> all_events = homoplaser.count_all_homoplasy_events(tree, seg_sites, segsite_indexes, physical_poss);  // cfml version
//        ArrayList<Homoplasy_Events> all_events = homoplaser.count_all_homoplasy_events(tree, seg_sites, physical_poss);  // treetime version
//
//        // EXAMINE itol tree to think about effect sizes
////        homoplaser.output_segsite_MRCAs_and_clades_as_branch_colors(all_events, 1044);
////        homoplaser.output_segsite_events_as_branch_colors(all_events, 1044);
////        System.exit(-1);
//
//  /*    // assess magnitude of recombination
//        // DEBUG:  is all_events sorted by phys pos?
////        all_events.forEach(s -> System.out.println(s.segsite_ID + "\t" + s.physical_pos));
//        homoplaser.identify_recombinant_segsites(all_events, phenos);
//        System.exit(-1); // stop here for debugging
//        homoplaser.calc_magnitude_recomb(all_events);
//
//        // REMOVE RECOMBINANT REGIONS from all_events
//        homoplaser.remove_all_recombinant_regions(all_events);
//*/
//
//        homoplaser.assoc(all_events, phenos);
//
////        homoplaser.assess_significance(all_events, phenos);  // DEPRECATED
//
//        // Summary:  calculate some stats and explore both the homoplasy events and ancestral reconstructions
////        homoplaser.characterize(all_events, tree);
//
//        // DEBUG
////        homoplaser.debug(all_events, seg_sites, segsite_indexes);
//
//        // TEST: jsc library
////        homoplaser.test_jsc();
//
//        // (deprecated) test "filesystem" version of "piping" code between Java <-> R
////        homoplaser.test_piping_fs_version();
//
//        // TEST "process-to-process" version (java process <-> Rscript process)
////        homoplaser. test_p_to_p_version();
//
//        homoplaser.end_session();
//    }


    /**
     * This method . . .
     */
    private void ancestral_reconstruction() {

        System.out.printf(Ansi.AUTO.string("@|fg(123) %nExecuting ancestral reconstruction . . . |@"));

        try {
            // treetime ancestral --aln (user-specified msa fasta) --tree (user-specified newick) --outdir ancestral_reconstruction --gtr infer
//        ProcessBuilder pb = new ProcessBuilder("treetime_no_command_available", "ancestral", "--aln", msa_fasta_filename, "--tree", user_input_newick_filename, "--outdir", anc_recon_dir, "--gtr", "infer");
            ProcessBuilder pb = new ProcessBuilder("treetime", "ancestral", "--aln", inputFiles.genotypes.msa_fasta_filename, "--tree", inputFiles.otherInputFiles.newick_filename, "--outdir", outputOptions.anc_recon_dir.getPath(), "--gtr", "infer");

            // DEBUG
//        pb = pb.inheritIO();
//        System.out.println("pb.command() = " + pb.command());
//        System.out.println("pb.environment() = " + pb.environment());

            outputOptions.log.write(String.format("%ntreetime command: " + pb.command() + "%n"));
            outputOptions.log.flush();

            Process anc_recon_process = pb.start();
            anc_recon_process.waitFor();

            if (DEBUG_MODE) {
//                System.out.print("process finished with exit value = ");
//                System.out.println(anc_recon_process.exitValue());
                outputOptions.debug.write("\nancestral reconstruction process finished with exit value = " + anc_recon_process.exitValue() + "\n");
                outputOptions.debug.flush();
            }

            // check to make sure that treetime completed without any errors, if no then fail fast to allow user to diagnose problem
            if (anc_recon_process.exitValue() != 0) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred during ancestral reconstruction (this is likely due to malformed input genotype and/or newick file formats).  You may try running treetime separately to better diagnose the treetime error.%n|@"));
                System.exit(-1);
            }
        } catch (IOException | InterruptedException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred during ancestral reconstruction. Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

//        System.out.println("exiting ancestral_reconstruction()");
        System.out.printf("reconstruction complete.%n");
    }


    /**
     * This method parses out a tree in newick format from treetime's nexus file.  Note that treetime annotates the newick format making it a
     * non-standard newick format.
     *
     * @return tree in newick format
     */
    private String nexus_to_newick() {
        // pseudocode:
        // read each line till "Begin Trees;".  next line is the newick line.
        // trim away front of newick line:  everything before first "("
        // StringBuilder(String newick line)
        // while (there are still annotations in the newick line):
        //      indexOf("[") -> int start
        //      indexOf("]") -> int end
        //      delete(start, end)

        String newick = null;

        try (BufferedReader br = new BufferedReader(new FileReader(outputOptions.anc_recon_dir + File.separator + "annotated_tree.nexus"))) {
            // moved BufferedReader into try-with-resources statement above
//            BufferedReader br = new BufferedReader(new FileReader(anc_recon_dir + File.separator + "annotated_tree.nexus"));
//            BufferedReader br = new BufferedReader(new FileReader(outputOptions.anc_recon_dir + File.separator + "annotated_tree.nexus"));

            String currline;
            while (!(currline = br.readLine()).equalsIgnoreCase("Begin Trees;")) {
            }

            StringBuilder newick_line = new StringBuilder(br.readLine());
//            System.out.println("newick_line = " + newick_line);

            // trim away the front of the newick line (typically contains: "Tree tree1=[&U]"), specifically everything before first "(" which begins newick tree
            int start_of_newick_tree = newick_line.indexOf("(");
            newick_line.delete(0, start_of_newick_tree);
//            System.out.println("after front trimming = " + newick_line);

            // remove all annotations from newick line (annotations appear as [& . . .]
            int start_annot;
            while ((start_annot = newick_line.indexOf("[")) != -1) {
                int end_annot = newick_line.indexOf("]");
                newick_line.delete(start_annot, end_annot + 1);  // end index is exclusive so +1 to delete "]"
//                System.out.println("one annot deleted = " + start_annot + " " + end_annot + ": " + newick_line);
            }

            newick = newick_line.toString();
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while processing the nexus file in the ancestral reconstruction directory. Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        // write newick tree out to file inside anc_recon_dir


        return newick;
    }


    /**
     * This method uses cfml's recombinant tract file to identify all segsites that fall within any of these regions, regardless of which sample the recomb tract is
     * detected in.  Here, I'm focusing on segregating sites, not mutations within a specific sample.
     *
     * @param all_events
     * @param phenos
     * @return
     */
    /*private void identify_recombinant_segsites(ArrayList<Homoplasy_Events> all_events, HashMap<String, String> phenos) {
        // boolean[] recombinant_segsites = new boolean[all_events.size()];

        *//*
        for (recombinant region, 1 row from cfml recombinant tract file) {
            - mark all elements from all_events that are within this region
        }
        *//*

        // mark wich segsites are covered by at least 1 recomb tract in at least 1 sample from cfml recomb tract file
        try (BufferedReader br = new BufferedReader(new FileReader(recomb_tracts_filename))) {
            String currline;
            br.readLine();  // skip header
            while ((currline = br.readLine()) != null) {  // 3 cols (tab-delimited):  sample_ID, start recomb tract, end
                String[] cols = currline.split("\t");
                String sample_ID = cols[0];
                int tract_start = Integer.parseInt(cols[1]);
                int tract_end = Integer.parseInt(cols[2]);

                // only use recomb tracts from extant nodes
                if (!phenos.containsKey(sample_ID)) {
                    continue;
                }

//                all_events.stream().filter(s -> s.physical_pos >= tract_start && s.physical_pos <= tract_end).forEach(s -> s.isRecombinant = true);
                all_events.stream().filter(s -> s.segsite_ID >= tract_start && s.segsite_ID <= tract_end).forEach(s -> s.isRecombinant = true);


                // DEBUG
//                System.out.println("Node_ID\ttract_start\ttract_end");
                if (tract_start < 526 || tract_end < 526) {
                    System.out.println(sample_ID + "\t" + tract_start + "\t" + tract_end);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        // DEBUG:  did the above work?
        List<Homoplasy_Events> recombinant_sites = all_events.stream().filter(s -> s.isRecombinant == true).collect(Collectors.toList());
        System.out.println("\n\n\ntot # recombinant segregating sites = " + recombinant_sites.size());
        recombinant_sites.forEach(r -> System.out.println(r.physical_pos));



        // return list of recombinant sites as a subset of all_events?
    }*/


    /**
     * maybe have this method be a general characterization of the cfml recomb tract file?  eg how many samples show recomb?  how many tracts?  mean tracts/sample?
     *
     * @param all_events
     */
    private void calc_magnitude_recomb(ArrayList<Homoplasy_Events> all_events) {
    }


    /**
     * I'd like to compare the maxT null dist with and without recombinant mutations.
     * Note:  this approach doesn't allow me to see how many sites are mixed sites (segsite composed of both independent and recombinant mutations), but this is quick & elegant
     * way to remove recombinant regions.  I can still calc the % segsites that are recombinant according to cfml.
     *
     * @param all_events
     */
    private void remove_all_recombinant_regions(ArrayList<Homoplasy_Events> all_events) {
    }


    private void end_session() {
        ZonedDateTime session_end_time = ZonedDateTime.now();
        Duration diff = Duration.between(session_start_time, session_end_time);
        String diff_pretty = String.format("%dh %02dm %02ds", diff.toHours(), diff.toMinutesPart(), diff.toSecondsPart());

        try {
            outputOptions.out.close();
            outputOptions.out_sorted_by_a1_maxT.close();
            outputOptions.out_sorted_by_a2_maxT.close();

            outputOptions.log.write(String.format("%nEnding poutine session: " + session_end_time.format(DateTimeFormatter.ofPattern("YYYY-LLLL-dd EEEE HH'h':mm'm':ss's' O")) + "%n"));
            outputOptions.log.write(String.format("Wall clock time elapsed = " + diff_pretty + "%n"));
            outputOptions.log.write(String.format("CLEAN EXIT%n"));
            outputOptions.log.close();

            if (DEBUG_MODE) {
                outputOptions.debug.write(String.format("%nEnding poutine session: " + session_end_time.format(DateTimeFormatter.ofPattern("YYYY-LLLL-dd EEEE HH'h':mm'm':ss's' O")) + "%n"));
                outputOptions.debug.write(String.format("Wall clock time elapsed = " + diff_pretty + "%n"));
                outputOptions.debug.write(String.format("CLEAN EXIT%n"));
                outputOptions.debug.close();
            }
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the log file. Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        System.out.printf(Ansi.AUTO.string("@|fg(213) %nEnding poutine session: " + session_end_time.format(DateTimeFormatter.ofPattern("YYYY-LLLL-dd EEEE HH'h':mm'm':ss's' O")) + "%n|@"));
        System.out.printf(Ansi.AUTO.string("@|fg(123) Wall clock time elapsed = " + diff_pretty + "%n|@"));
        System.out.printf(Ansi.AUTO.string("@|bold,blink,fg(197) CLEAN EXIT%n|@"));
    }


    /**
     * Outputs variables that helps to understand and debug the program
     */
    private void log_session_info() {
        // TEST:  jansi library
//        systemInstall();
//        System.out.printf(ansi().render("%n%n@|bold,underline,white Homoplasy counting session log -----------%n|@").toString());

        // TEST:  picocli color library
        System.out.printf(Ansi.AUTO.string("%n%n@|yellow,bold,bg(222) Homoplasy counting session log -----------|@%n"));

        // ORIGINAL OUTPUT CODE
        //        System.out.printf("%n%nHomoplasy counting session log -----------%n");

        // log session start
        System.out.println("session start time:  " + session_start_time);

        // log important global variables
        System.out.println("m = " + m);
        System.out.println("min_hcount = " + AlgoParams.min_hcount);
        System.out.println("R_SPACE_TYPE = " + R_SPACE_TYPE);
        System.out.println("TAIL_TYPE = " + TAIL_TYPE);
        System.out.println("EXTANT_NODES_ONLY = " + EXTANT_NODES_ONLY);
        System.out.println("tot_num_sample_cases = " + tot_num_sample_cases + "    tot_num_sample_controls = " + tot_num_sample_controls + "    sample size = " +
                (tot_num_sample_cases + tot_num_sample_controls));
//        ORIGINAL OUTPUT CODE
//        System.out.printf("End log -------------%n%n");

//        System.out.printf(ansi().bg(Color.valueOf("222")).render("@|white End log -------------%n%n|@").toString());
        System.out.printf(Ansi.AUTO.string("@|white,bg(213) End log -------------|@%n%n"));
//        systemUninstall();
    }



/*    //    private void test_piping_fs_version() {
    private void test_p_to_p_version() {
//        ProcessBuilder pb = new ProcessBuilder("echo", "What would Snowden do? . . .");
        ProcessBuilder pb = new ProcessBuilder("Rscript", "--vanilla", "--verbose", "piping_test.R");

        try {
            Process p = pb.start();

            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
            bw.write("10,1,1,10");
            bw.newLine();
            bw.write("60,0,10,60");
            bw.newLine();
            bw.write("50,8,5,49");
            bw.newLine();
            bw.close();

            int err_code = p.waitFor();
            System.out.println("err_code = " + err_code);

            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String subprocess_output;
            while ((subprocess_output = br.readLine()) != null) {
                System.out.println("subprocess_output = " + subprocess_output);
            }

            System.out.println("\n\nEnd Test.\n\n");
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }*/



/*
    private void test_jsc() {
//        int[][] test_site_homoplasy_counts = new int[2][2];
        int[][] test_site_homoplasy_counts = {{0,10}, {10,0}};
        ContingencyTable2x2 test_site_contingency_table = new ContingencyTable2x2(test_site_homoplasy_counts);
        FishersExactTest test_statistic = new FishersExactTest(test_site_contingency_table);

        System.out.println("test_site_contingency_table = " + test_site_contingency_table);
        System.out.println("test_statistic.getTestStatistic() = " + test_statistic.getTestStatistic());
        System.out.println("test_statistic.getOneTailedSP() = " + test_statistic.getOneTailedSP());
        System.out.println("test_statistic.getSP() = " + test_statistic.getSP());
        System.out.println("test_statistic.getOppositeTailProb() = " + test_statistic.getOppositeTailProb());

    }
*/


    // NOTE:  no longer need an explicit constructor since moving to picocli library
//    private Homoplasy_Counter(String[] args) {
////        if (args.length != 6) {
//        if (args.length != 5) {
//
//            System.out.println("\nUsage:  java Homoplasy_Counter . . . ");
//            System.exit(-1);
//        }
//
//        // TODO:  code in command-line switches
////        newick_filename = args[0];
//
////        nexus_filename = args[0];
////        ancestral_reconstruction_filename = args[1];
//
//        msa_fasta_filename = args[0];
//        user_input_newick_filename = args[1];
//
////        pos_cross_ref_filename = args[2];
//        pheno_filename = args[2];
//        map_filename = args[3];
////        recomb_tracts_filename = args[5];
//        qvalues_option = Boolean.parseBoolean(args[4]);
//    }


    /**
     * Build a tree data structure from a newick string.
     */
    private NewickTree build_tree(String newick) {

        NewickTree tree = null;

        try {
            StringReader sr = new StringReader(newick);
            NewickTreeTokenizer tokenizer = new NewickTreeTokenizer(sr);
            NewickTreeReader tree_reader = new NewickTreeReader(tokenizer);
            tree = tree_reader.readTree();

            outputOptions.log.write(String.format("%ntree file: # of leaves = " + tree.getLeafNodes().size() + "%n"));
            outputOptions.log.flush();
        } catch (IOException | DataFormatException | NullPointerException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while processing the input newick tree. Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

//        System.out.println("phylogeny:  ");
//        System.out.println("# of leaves = " + tree.getLeafNodes().size());
//        int num_children_of_rootnode = tree.getRoot().getChildrenCount();
//        System.out.println("num_children_of_rootnode = " + num_children_of_rootnode);
//        System.out.println();

        return tree;
    }


    /**
     * This method works with cfml version.
     *
     * Build a tree data structure from newick grammar
     */
    /*private NewickTree build_tree() {
        NewickTree tree = null;

        try (BufferedReader br = new BufferedReader(new FileReader(newick_filename))) {
            // read in a newick tree
            NewickTreeTokenizer tokenizer = new NewickTreeTokenizer(br);
            NewickTreeReader tree_reader = new NewickTreeReader(tokenizer);
            tree = tree_reader.readTree();
        } catch (IOException | DataFormatException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        // TODO:  refactor into log code
        System.out.println("phylogeny:  ");
        System.out.println("# of leaves = " + tree.getLeafNodes().size());
//        int num_children_of_rootnode = tree.getRoot().getChildrenCount();
//        System.out.println("num_children_of_rootnode = " + num_children_of_rootnode);
        System.out.println();

        return tree;
    }*/


    /**
     * // TODO:  rewrite this in a functional programming style using Streams and BufferedReader.lines().  How does this compare to the imperative programming version?
     * <p>
     * Process PLINK's .map file and extract physical positions for each seg site.  The position should always be the last column in the .map file (the # of cols
     * may vary according to PLINK usage)
     */
    private int[] get_physical_positions() {
        ArrayList<Integer> poss = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(inputFiles.otherInputFiles.map_filename))) {
            String currline;
            while ((currline = br.readLine()) != null) {
                String[] cols = currline.split("\t");
                int p = Integer.parseInt(cols[cols.length - 1]);
                poss.add(p);
            }

            outputOptions.log.write(String.format("%nplink map file: # segregating sites = " + poss.size() + "%n"));
            outputOptions.log.flush();
//        } catch (FileNotFoundException e) {
//            e.printStackTrace();
//            System.exit(-1);
        } catch (Exception e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while reading the map file (this is likely due to a malformed map file). Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }
//        catch (NumberFormatException nfe) {
//            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while parsing the map file (this is likely due to a malformed map file). Please see the stack trace in the log file for more information on this error.%n|@"));
//            nfe.printStackTrace();
//            System.exit(-1);
//        }

        return poss.stream().mapToInt(i -> i).toArray();
    }


    /**
     * This method reads in treetime's ancestral fasta file or a user-specified ancestral fasta file if --use-precomputed-anc-recon option is turned on.
     */
    private HashMap<String, char[]> build_seg_sites(String ancestral_fasta) {
        HashMap<String, char[]> seg_sites = new HashMap<>();  // key := node name, value := snps across all segsites
//        ArrayList<String> node_names = new ArrayList<>();
//        ArrayList<char[]> snp_sets = new ArrayList<>();

        try {
            int num_segsites_detected = -1;
//            int num_seqs = 0;

    //        Fasta_Manager multi_fasta_file = new Fasta_Manager(anc_recon_dir + File.separator + "ancestral_sequences.fasta");
//            Fasta_Manager multi_fasta_file = new Fasta_Manager(outputOptions.anc_recon_dir + File.separator + "ancestral_sequences.fasta");
            Fasta_Manager multi_fasta_file = new Fasta_Manager(ancestral_fasta);
            Fasta_Record curr_record;
            while ((curr_record = multi_fasta_file.next()) != null) {
                String node_name = curr_record.getHeader();
                char[] seq = curr_record.getSequence().toCharArray();

                // TEST missing geno data
//                if (node_name.equals("M230")) {
//                    seq[1044 - 1] = 'N';  // spike this sample with a missing geno at segsite 1044 (the top hit in ref set).  Genotype for sample M230 at site 1044 is homoplasic.
//                }
//                if (node_name.equals("KZNXDR")) {
//                    seq[1044 - 1] = 'N';  // spike this sample with a missing geno at segsite 1044 (the top hit in ref set).  Genotype for sample KZNXDR at site 1044 is homoplasic.
//                }

                seg_sites.put(node_name, seq);

                num_segsites_detected = seq.length;  // TODO:  check to make sure every single fasta seq has the same # sites
                // DEBUG
    //            System.out.println(node_name + "\t" + seq.length);
//                num_seqs++;
            }

        /*// replace this section with code from my bioinformatics library to process fasta files: =======================================
        try (BufferedReader br = new BufferedReader(new FileReader(ancestral_reconstruction_filename))) {
            String currline;
            while ((currline = br.readLine()) != null) {
                if (currline.startsWith(">")) {
                    String node_name = currline.substring(1).trim();
                    node_names.add(node_name);
                } else {
                    char[] snps = currline.toCharArray();  // snps across all segsites for one sample/internal node
                    snp_sets.add(snps);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        for (int i = 0; i < node_names.size(); i++) {
            seg_sites.put(node_names.get(i), snp_sets.get(i));
        }  // ===========================================================================================================================*/

//            outputOptions.log.write(String.format("%n# samples detected from input genotype file: " + seg_sites.size() + "%n"));
//                outputOptions.log.write(String.format("%ngenotype file: # samples = " + seg_sites.size() + ", # segregating sites in each sample = " + num_segsites_detected + "%n"));
                outputOptions.log.write(String.format("%ngenotype file: # segregating sites in each sample = " + num_segsites_detected + "%n"));
                outputOptions.log.flush();
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while processing the ancestral_sequences.fasta file in the ancestral reconstruction directory. Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        return seg_sites;
    }


    /**
     * Read in .ML_sequence.fasta file.  Each fasta entry is a set of snps where the column represents the seg site
     * TODO:  use my FastaManager or write a new one to properly process any kind of fasta file
     *
     * NOTES:
     * - seg_sites[] is 0-based
     *
     * - this CFML output file contains only unique cols (this is an artifact of CFML, why they coded it this way is beyond me).
     *   Nevertheless, the cols in this file do not map 1-to-1 to any of the other files (map, ped, fasta of seg sites, etc).
     */
    /*private HashMap<String, char[]> build_seg_sites() {
        HashMap<String , char[]> seg_sites = new HashMap<>();  // key := node name, value := snps across all segsites in cfml format (i.e. only non-redundant sites are in the fasta file)
        ArrayList<String> node_names = new ArrayList<>();
        ArrayList<char[]> snp_sets = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(ancestral_reconstruction_filename))) {
            String currline;
            while ((currline = br.readLine()) != null) {
                if (currline.startsWith(">")) {
                    String node_name = currline.substring(1).trim();
                    node_names.add(node_name);
                } else {
                    char[] snps = currline.toCharArray();
                    snp_sets.add(snps);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        for (int i = 0; i < node_names.size(); i++) {
            seg_sites.put(node_names.get(i), snp_sets.get(i));
        }

        return seg_sites;
    }*/


    /**
     * Read in .position_cross_reference.txt.
     *
     * Values in segsite_indexes are 1-based according to CFML.
     *
     * This CFML file maps the CFML fasta file of ancestral reconstructions to other files that contain segsites (e.g. ped).
     * It is this CFML file that maps 1-to-1 (col-to-col) to other segsite files.  As such, the # cols in this file equals the tot # segsites.
     */
    /*private int[] get_segsite_indexes() {
        int[] segsite_indexes = null;

        try (BufferedReader br = new BufferedReader(new FileReader(pos_cross_ref_filename))) {
            String[] indexes_as_Strings = br.readLine().split(",");
            segsite_indexes = new int[indexes_as_Strings.length];

            for (int i = 0; i < segsite_indexes.length; i++) {
                segsite_indexes[i] = Integer.parseInt(indexes_as_Strings[i]);
            }
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        // sanity check:  # indexes should equal # extant seg sites
//        System.out.println("segsite_indexes.length = " + segsite_indexes.length + " (should equal # extant seg sites)");

        return segsite_indexes;
    }*/


    /**
     * Maps all snps of this seg site to the topology.
     * <p>
     * segsite_index is 0-based.
     */
    private HashMap<String, Character> map_seg_site_to_tree(HashMap<String, char[]> seg_sites, int segsite_index) {
        HashMap<String, Character> mapped_seg_site = new HashMap<>();  // key := node name, value := allele

        for (Map.Entry<String, char[]> mapping : seg_sites.entrySet()) {
            String node_name = mapping.getKey();
            char[] snps = mapping.getValue();

            mapped_seg_site.put(node_name, snps[segsite_index]);
        }

        return mapped_seg_site;
    }


    /**
     * THIS VERSION WORKS WITH TREETIME'S ANCESTRAL RECONSTRUCTION OUTPUT.
     * <p>
     * Process 1 seg site at a time (i.e. count homoplasy events for 1 seg site at a time).
     */
    private ArrayList<Homoplasy_Events> count_all_homoplasy_events(NewickTree tree, HashMap<String, char[]> seg_sites, int[] physical_poss) {

        System.out.printf(Ansi.AUTO.string("@|fg(222) %nCounting all homoplasy events . . . |@"));

//        int num_non_biallelic_sites = 0;
        int num_monomorphic_sites = 0;  // allele counts are specific to the sample population, and thus does not include the allele from the reference genome, e.g. monomorphic is not polymorphic in your sample pop.
        int num_bi_allelic_sites = 0;
        int num_tri_allelic_sites = 0;
        int num_quad_allelic_sites = 0;
        int more_than_four_alleles_detected = 0;
        ArrayList<Homoplasy_Events> all_events = new ArrayList<>();

        List<NewickTreeNode> leaves = tree.getLeafNodes();
        HashSet<String> leaf_names = get_leaf_names(leaves);

//        for (int i = 0; i < segsite_indexes.length; i++) {  // for (each seg site)
        for (int i = 0; i < physical_poss.length; i++) {  // for (each seg site)
            int segsite_ID = i + 1;  // this ID is 1-based, e.g. segsite_ID = 10 is the 10th segsite
            int p = physical_poss[i];

//            int segsite_index = segsite_indexes[i] - 1;  // segsite_index is 0-based (while values in segsite_indexes[] from position_cross_references.txt are 1-based)

            // NOTE:  below was for cfml version.  No longer checking for monomorphic sites for now, may have to add code if this arises again.
            // workaround for monomorphic sites
//            if (segsite_index == -1) {
//                num_monomorphic_sites++;
//                continue;
//            }

//            HashMap<String,Character> mapped_seg_site = map_seg_site_to_tree(seg_sites, segsite_index);
            HashMap<String, Character> mapped_seg_site = map_seg_site_to_tree(seg_sites, i);
            // TODO:  rename method to map_ancestral_alleles_to_seg_site()?  rename mapped_seg_site to just seg_site?

            /* NOTE:  Blocked this out to expand code to identify monomorphic, tri-allelic, quad-allelic polymorphic sites.  Instead of calling is_biallelic(), # of alleles is identified and returned in new method.
            if (is_biallelic(leaves, mapped_seg_site)) {
                Homoplasy_Events segsite_homoplasy_events = identify_homoplasy_events_for_seg_site(tree, mapped_seg_site, segsite_ID, p, leaf_names);
                all_events.add(segsite_homoplasy_events);
            } else {  // tri/quad-allelic site
                num_non_biallelic_sites++;
            }*/

            Map<Character, List<Character>> alleles = get_alleles(leaves, mapped_seg_site);
//            switch (num_alleles_at_site(leaves, mapped_seg_site)) {
            switch (alleles.size()) {
                case 2:  // bi-allelic site
                    num_bi_allelic_sites++;
                    Homoplasy_Events segsite_homoplasy_events = identify_homoplasy_events_for_seg_site(tree, mapped_seg_site, segsite_ID, p, alleles, leaves, leaf_names);
                    all_events.add(segsite_homoplasy_events);

                    // TEST missing geno code at top hit site 1044 in ref set
//                    if (p == 1473246) { // if (segsite 1044)
//                        System.out.println("\n\nsegsite_homoplasy_events.physical_pos = " + segsite_homoplasy_events.physical_pos);
//                        System.out.println("segsite_homoplasy_events.segsite_ID = " + segsite_homoplasy_events.segsite_ID);
                        // these are the 9 homoplasic cases by sample name
//                        System.out.println("segsite_homoplasy_events.allele2_MRCAs.keySet() = " + segsite_homoplasy_events.allele2_MRCAs.keySet());
//                        System.exit(-1);
//                    }

                    break;
                case 1:  // monomorphic site
                    num_monomorphic_sites++;
                    break;
                case 3:  // tri-allelic
                    num_tri_allelic_sites++;
                    break;
                case 4:  // quad-allelic
                    num_quad_allelic_sites++;
                    break;
                default:
                    more_than_four_alleles_detected++;
                    break;
            }
        }

        try {
            outputOptions.log.write(String.format("%n# bi-allelic sites: " + num_bi_allelic_sites));
            outputOptions.log.write(String.format("%n# monomorphic sites: " + num_monomorphic_sites));
            outputOptions.log.write(String.format("%n# tri-allelic sites: " + num_tri_allelic_sites));
            outputOptions.log.write(String.format("%n# quad-allelic sites: " + num_quad_allelic_sites + "%n"));
//            outputOptions.log.write(String.format("%n# sites with more than four alleles detected: " + more_than_four_alleles_detected + "%n"));  // report only in debug file as a sanity check
            outputOptions.log.flush();
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the log file.  Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        if (DEBUG_MODE) {
            try {
                outputOptions.debug.write("\ncount_all_homoplasy_events():\n");
//                outputOptions.debug.write("num_non_biallelic_sites = " + num_non_biallelic_sites + "\n");
                outputOptions.debug.write("# bi-allelic sites = " + num_bi_allelic_sites + "\n");
                outputOptions.debug.write("# monomorphic sites = " + num_monomorphic_sites + "\n");
                outputOptions.debug.write("# tri-allelic sites = " + num_tri_allelic_sites + "\n");
                outputOptions.debug.write("# quad-allelic sites = " + num_quad_allelic_sites + "\n");
                outputOptions.debug.write("sites with more than four alleles detected = " + more_than_four_alleles_detected + "\n");
                outputOptions.debug.flush();
//                System.out.println("num_non_biallelic_sites = " + num_non_biallelic_sites);
//                System.out.println("num_monomorphic_sites = " + num_monomorphic_sites);
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing to the debug file. Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }

        System.out.printf("homoplasy identification complete.%n");

        return all_events;
    }



    /**
     * THIS VERSION IS THE OLDER VERSION USED TO PROCESS CFML OUTPUT FILES.  THIS HAS NOW BEEN REPLACED BY THE METHOD OF THE SAME NAME WHICH NOW
     * PROCESSES TREETIME'S OUTPUT FILES.
     * <p>
     * Process 1 seg site at a time (i.e. count homoplasy events for 1 seg site at a time).
     */
    /*private ArrayList<Homoplasy_Events> count_all_homoplasy_events(NewickTree tree, HashMap<String, char[]> seg_sites, int[] segsite_indexes, int[] physical_poss) {
        int num_non_biallelic_sites = 0;
        int num_monomorphic_sites = 0;
        ArrayList<Homoplasy_Events> all_events = new ArrayList<>();

        // get leaf names (this is a bit inelegant as it is not used until new Homoplasy_Events()).  This is the closest I can place the code to where it is needed without
        // calling this code for each segsite (I only want to get the leaf names once for all sites).
        HashSet<String> leaf_names = get_leaf_names(tree);

        // TODO:  refactor all calls to get leaves such that it relies on this call here (leaves never change so there is no point calling the function for each segsite)
        List leaves = tree.getLeafNodes();

        for (int i = 0; i < segsite_indexes.length; i++) {  // for (each extant seg site)
            int segsite_ID = i + 1;  // this ID is 1-based, e.g. segsite_ID = 10 is the 10th segsite
            int p = physical_poss[i];

            int segsite_index = segsite_indexes[i] - 1;  // segsite_index is 0-based (while values in segsite_indexes[] from position_cross_references.txt are 1-based)

            // workaround for monomorphic sites
            if (segsite_index == -1) {
                num_monomorphic_sites++;
                continue;
            }

            HashMap<String,Character> mapped_seg_site = map_seg_site_to_tree(seg_sites, segsite_index);

            // TODO:  perhaps refactor elsewhere, quick do this now to get results
            if (is_biallelic(leaves, mapped_seg_site, segsite_ID)) {
                Homoplasy_Events segsite_homoplasy_events = identify_homoplasy_events_for_seg_site(tree, mapped_seg_site, segsite_ID, p, leaf_names);
                all_events.add(segsite_homoplasy_events);
            } else {  // tri/quad-allelic site
                num_non_biallelic_sites++;
            }
        }

        System.out.println("num_non_biallelic_sites = " + num_non_biallelic_sites);
        System.out.println("num_monomorphic_sites = " + num_monomorphic_sites);

        return all_events;
    }*/



    private boolean is_biallelic(List<NewickTreeNode> leaves, HashMap<String, Character> mapped_seg_site) {
        boolean is_biallelic = false;

        ArrayList<Character> all_alleles = new ArrayList<>();
        for (NewickTreeNode leaf : leaves) {
            all_alleles.add(get_node_genotype(leaf, mapped_seg_site));
        }
        Map<Character, List<Character>> alleles = all_alleles.stream().collect(Collectors.groupingBy(a -> a));

        if (alleles.size() == 2) {
            is_biallelic = true;
        } else {
            is_biallelic = false;
        }

        // TODO:  check for non-ACGT chars


        return is_biallelic;
    }


    /**
     * This method identifies the alleles at this segregating site.
     *
     * @param leaves
     * @param mapped_seg_site
     * @return Map<Character allele, List genotypes>.  Genotypes are grouped by allele.
     */
    private Map<Character, List<Character>> get_alleles(List<NewickTreeNode> leaves, HashMap<String, Character> mapped_seg_site) {
        ArrayList<Character> genotypes = new ArrayList<>();

        for (NewickTreeNode leaf : leaves)
            genotypes.add(get_node_genotype(leaf, mapped_seg_site));

        // ignore all missing genotype data (any geno codings outside of [ACGT])
        Map<Character, List<Character>> alleles = genotypes.stream().filter(geno -> geno == 'A' || geno == 'C' || geno == 'G' || geno == 'T').collect(Collectors.groupingBy(g -> g));

        return alleles;
    }



    private int num_alleles_at_site(List<NewickTreeNode> leaves, HashMap<String, Character> mapped_seg_site) {

        ArrayList<Character> all_genotypes_at_site = new ArrayList<>();

        for (NewickTreeNode leaf : leaves)
            all_genotypes_at_site.add(get_node_genotype(leaf, mapped_seg_site));

        Map<Character, List<Character>> alleles = all_genotypes_at_site.stream().collect(Collectors.groupingBy(g -> g));

        return alleles.size();
    }



    /**
     * Get the names (ie sample ID) of each leaf.
     *
     * @param tree
     * @return
     */
    private HashSet<String> get_leaf_names(List<NewickTreeNode> leaves) {
        HashSet<String> leaf_names = new HashSet<>();

//        List<NewickTreeNode> leaves = tree.getLeafNodes();
        leaves.forEach(leaf -> leaf_names.add(leaf.getId()));

        return leaf_names;
    }


    /**
     * Traverse entire tree and identify homoplasy events for this particular seg site.
     * Traversals start from each leaf all the way to the root, depth first.
     * <p>
     * already_visited_nodes is an optimization to insure that each node is visited only once.
     * Currently, optimization not in use; keep code intuitive for now.
     * <p>
     * This method identifies and stores all MRCAs for each allele.  Only when # MRCAs for an allele >= 2 is there a homoplasy event by definition.
     * // TODO:  this function could use parallization!
     */
    private Homoplasy_Events identify_homoplasy_events_for_seg_site(NewickTree tree, HashMap<String, Character> mapped_seg_site, int segsite_ID, int p,
                                                                    Map<Character, List<Character>> alleles, List<NewickTreeNode> leaves, HashSet<String> leaf_names) {
        HashMap<String, HashSet<String>> allele1_MRCAs = new HashMap<>();  // key := node name of MRCA, value := members of this MRCA: HashSet<node names>
        HashMap<String, HashSet<String>> allele2_MRCAs = new HashMap<>();  // same as above
//        HashMap<String,String> visited_nodes = new HashMap<>();  // key =: node name, value := node name of MRCA set

        // first, identify all the MRCAs for both alleles:
        NewickTreeNode root = tree.getRoot();
//        List leaves = tree.getLeafNodes();

        // deprecated:  this is where the a1 allele is assigned arbitrarily, let's add some code to figure out which allele is the major allele, then assign that to a1
        // 2 lines commented out so that allele1 is now always the major allele:
        //        NewickTreeNode first_leaf = (NewickTreeNode) leaves.get(0);  // arbitrarly pick the first leaf allele; likely the major allele
        //        char allele1 = get_node_genotype(first_leaf, mapped_seg_site);
//        char major_allele = identify_major_allele(leaves, mapped_seg_site, segsite_ID);
        HashMap<String, Character> major_minor_alleles = identify_major_minor_alleles(alleles);  // at this point, alleles must be biallelic
//        char allele1 = major_allele;
        char allele1 = major_minor_alleles.get("a1");

        for (int i = 0; i < leaves.size(); i++) {
//            NewickTreeNode leaf = (NewickTreeNode) leaves.get(i);
            NewickTreeNode leaf = leaves.get(i);
            ArrayList<NewickTreeNode> curr_path = new ArrayList<>();  // nodes visited starting from this leaf
            NewickTreeNode prev_node = leaf;

//            curr_path.add(leaf);  // No longer using HashSet clade data structure
            // leaf -> root
            for (NewickTreeNode curr_node = leaf.getParent(); !curr_node.equals(root); curr_node = curr_node.getParent()) {  // TODO:  rewrite equals() to curr_node != root ?
//                curr_path.add(curr_node);

                if (is_MRCA_found(curr_node, prev_node, mapped_seg_site)) {
                    // prev_node is the MRCA
                    add_MRCA(prev_node, mapped_seg_site, allele1, allele1_MRCAs, allele2_MRCAs, curr_path);
                }
//                else {
//                    curr_path.add(curr_node);  // only add non-MRCA nodes.  No longer using Hashset clade data structure
//                }

                prev_node = curr_node;
            }

            // at this point:  root was visited
            add_MRCA(root, mapped_seg_site, allele1, allele1_MRCAs, allele2_MRCAs, curr_path);
        }

        // at this point:  curr seg site has been fully counted.  The 2 allele HashMaps currently store all identified MRCAs.  Check each allele HashMap and remove all non-homoplasic mutations.
        remove_non_homoplasic_mutations(allele1_MRCAs);
        remove_non_homoplasic_mutations(allele2_MRCAs);

//        return new Homoplasy_Events(allele1_MRCAs, allele2_MRCAs, segsite_ID, p, alleles, leaves);
        return new Homoplasy_Events(allele1_MRCAs, allele2_MRCAs, segsite_ID, p, major_minor_alleles, leaf_names);
    }



    /**
     * This method identifies the major and minor alleles at this segregating site.
     * This method is only for biallelic sites.
     *
     * @param leaves
     * @param mapped_seg_site
     * @param segsite_id
     * @param alleles
     * @return HashMap<String, Character> alleles, key := "a1" (major) or "a2" (minor), value := allele [ACGT]
     */
    private HashMap<String, Character> identify_major_minor_alleles(Map<Character, List<Character>> alleles) {

        char major_allele = 0;
        char minor_allele = 0;
//        int major_allele_count = 0;
//        int minor_allele_count = 0;
        Iterator<Map.Entry<Character, List<Character>>> it = alleles.entrySet().iterator();  // biallelic sites have exactly 2 Map entries, one for each allele
        Map.Entry<Character, List<Character>> arbitrary_a1 = it.next();
        Map.Entry<Character, List<Character>> arbitrary_a2 = it.next();
        int arbitrary_a1_size = arbitrary_a1.getValue().size();
        int arbitrary_a2_size = arbitrary_a2.getValue().size();
        HashMap<String, Character> major_minor_alleles = new HashMap<>();

        if (arbitrary_a1_size >= arbitrary_a2_size) {  // for equal frequency alleles (a1 size == a2 size), arbitrary_a1 is picked as the major allele
            major_allele = arbitrary_a1.getKey();
//            major_allele_count = arbitrary_a1_size;

            minor_allele = arbitrary_a2.getKey();
//            minor_allele_count = arbitrary_a2_size;
        } else {
            major_allele = arbitrary_a2.getKey();
//            major_allele_count = arbitrary_a2_size;

            minor_allele = arbitrary_a1.getKey();
//            minor_allele_count = arbitrary_a1_size;
        }

        major_minor_alleles.put("a1", major_allele);
        major_minor_alleles.put("a2", minor_allele);

        return major_minor_alleles;
    }



    /**
     * DEPRECATED:  replaced with identify_major_minor_allele()
     *
     * Returns the major allele at this site.
     *
     * @param leaves
     * @param mapped_seg_site
     * @param segsite_ID
     * @return
     */
    private char identify_major_allele(List<NewickTreeNode> leaves, HashMap<String, Character> mapped_seg_site, int segsite_ID) {

        // 2 sanity checks:  1)  there must be exactly l # of leaves, where l := # samples   2) there must be exactly 2 alleles at each site
//        if (leaves.size() != )

        ArrayList<Character> all_alleles = new ArrayList<>();

        for (NewickTreeNode leaf : leaves) {
            all_alleles.add(get_node_genotype(leaf, mapped_seg_site));
        }

        Map<Character, List<Character>> alleles = all_alleles.stream().collect(Collectors.groupingBy(a -> a));

        if (alleles.size() != 2) {
            System.out.println("\nError:  there are more than 2 alleles at this segregating site");
            System.out.println("segsite_ID = " + segsite_ID);
            System.out.println("alleles = " + alleles);
            // counts for each allele (should total # of samples actually used)
            alleles.values().forEach(a -> System.out.println(a.size()));
            System.exit(-1);
        }

        char major_allele = 0;
        char minor_allele = 0;
        int major_allele_count = 0;
        int minor_allele_count = 0;
        Iterator<Map.Entry<Character, List<Character>>> it = alleles.entrySet().iterator();
        Map.Entry<Character, List<Character>> arbitrary_a1 = it.next();
        Map.Entry<Character, List<Character>> arbitrary_a2 = it.next();

        if (arbitrary_a1.getValue().size() > arbitrary_a2.getValue().size()) {
            major_allele = arbitrary_a1.getKey();
            major_allele_count = arbitrary_a1.getValue().size();

            minor_allele = arbitrary_a2.getKey();
            minor_allele_count = arbitrary_a2.getValue().size();
        } else {
            major_allele = arbitrary_a2.getKey();
            major_allele_count = arbitrary_a2.getValue().size();

            minor_allele = arbitrary_a1.getKey();
            minor_allele_count = arbitrary_a1.getValue().size();
        }

        // DIAGNOSTICS:
//        NewickTreeNode first_leaf = (NewickTreeNode) leaves.get(0);  // arbitrarly pick the first leaf allele; likely the major allele
//        char a1 = get_node_genotype(first_leaf, mapped_seg_site);  // this is what a1 is currently in my code

//        if (a1 != major_allele) {
//            System.out.println("picked wrong major allele");
//        }

//        System.out.println("segsite_ID = " + segsite_ID);
//        System.out.println("a1 = " + a1);
//        System.out.println("major_allele = " + major_allele + "\tmajor_allele_count = " + major_allele_count);
//        System.out.println("minor_allele = " + minor_allele + "\tminor_allele_count = " + minor_allele_count);
//        System.out.println("leaves.size() = " + leaves.size());
//        System.out.println();

        return major_allele;
    }


    /**
     * Return the genotype mapped to this node.
     */
    private char get_node_genotype(NewickTreeNode node, HashMap<String, Character> mapped_seg_site) {
        // DEBUG
//        if (!mapped_seg_site.containsKey(node.getId())) {  // INTERNAL NODE NAME USED HERE
//            System.out.println("missing sample from input segsites fasta file (but present in leaves of input tree) = " + node.getId());
//        }

        return mapped_seg_site.get(node.getId());
    }


    /**
     * Returns true if the MRCA has been found.  Otherwise false.
     * Note that if the prev_node (the MRCA of interest) is a missing genotype value (any value not [ACGT]) then method returns false.
     */
    private boolean is_MRCA_found(NewickTreeNode curr_node, NewickTreeNode prev_node, HashMap<String, Character> mapped_seg_site) {
        boolean is_MRCA_found = false;

        char curr_geno = get_node_genotype(curr_node, mapped_seg_site);
        char prev_geno = get_node_genotype(prev_node, mapped_seg_site);

        // never consider missing genos in homoplasy identification
        if (prev_geno != 'A' && prev_geno != 'C' && prev_geno != 'G' && prev_geno != 'T') {  // if this node corresponds to a missing genotype
            is_MRCA_found = false;
        } else if (curr_geno != prev_geno) {
            is_MRCA_found = true;
        } // else these two genotypes are non-missing and equal (i.e. prev_node is not homoplasic)

//        return curr_geno != prev_geno;
        return is_MRCA_found;
    }



    private void add_MRCA(NewickTreeNode MRCA_node, HashMap<String, Character> mapped_seg_site, char allele1, HashMap<String, HashSet<String>> allele1_MRCAs,
                          HashMap<String, HashSet<String>> allele2_MRCAs, ArrayList<NewickTreeNode> curr_path) {

        char MRCA_geno = get_node_genotype(MRCA_node, mapped_seg_site);

        if (MRCA_geno == allele1) {
            update_MRCAs(allele1_MRCAs, MRCA_node, curr_path);
        } else {  // a2
            update_MRCAs(allele2_MRCAs, MRCA_node, curr_path);
        }
    }



    private void update_MRCAs(HashMap<String, HashSet<String>> MRCAs, NewickTreeNode MRCA_node, ArrayList<NewickTreeNode> node_path) {
        String MRCA_node_name = MRCA_node.getId();

        /* HashSet clade is only used to draw out subclades of MRCA nodes, removing its use for now as this is bloat without providing utility beyond the drawing
        if (!MRCAs.containsKey(MRCA_node_name)) {  // new MRCA to be added
            HashSet<String> clade = new HashSet<>();
            add_to_clade(clade, node_path);
            MRCAs.put(MRCA_node_name, clade);
        } else {  // add clade members to this MRCA
            HashSet<String> clade = MRCAs.get(MRCA_node_name);

            // add new nodes from node_path to this MRCA's clade.  Adding node names only.
            add_to_clade(clade, node_path);
        }*/

        if (!MRCAs.containsKey(MRCA_node_name)) {
            MRCAs.put(MRCA_node_name, null);
        }
    }


    /**
     * Add all nodes in ArrayList node_path to HashSet clade.
     */
    private void add_to_clade(HashSet<String> clade, ArrayList<NewickTreeNode> node_path) {
        for (int i = 0; i < node_path.size(); i++) {
            clade.add(node_path.get(i).getId());  // if adding nodes themselves, I can use HashSet.addAll(node_path) without looping
        }
    }


    /**
     * DEPRECATED:  combining this method with identify_major_minor_alleles().
     *
     * This method simply figures out the two alleles at this seg site.
     * a1 := allele 1
     * a2 := allele 2
     * <p>
     * Also, allele1 must be the same allele used for allele1_MRCAs (i.e. so we don't have allele flipping)
     */
    private HashMap<String, Character> get_alleles_for_site(char allele1, HashMap<String, Character> mapped_seg_site) {
        HashMap<String, Character> alleles = new HashMap<>();

/*
        // IMPERATIVE version
        for (Character snp : mapped_seg_site.values()) {
            if (snp == allele1) {  // can I use == operator here?
                alleles.put("a1", snp);
            } else {
                alleles.put("a2", snp);
            }
        }
*/

        // FUNCTIONAL version 1
        // is there a way to do this without explicit casting?  Why does the type not carry through the Stream?
        // also likely a better way to do this without creating a new object (toArray()), i.e. maybe use a lambda to insert proper allele into HashMap?
//        Character[] alleles_unordered = mapped_seg_site.values().stream().distinct().toArray();
//        if (alleles_unordered[0] == allele1) {
//            alleles.put("a1", allele1);
//        } else {
//            alleles.put("a2", alleles_unordered[0]);
//        }
//
//         sanity check (length must equal 2 alleles)
//        if (alleles_unordered.length != 2) {
//            System.out.println("alleles_unordered.length = " + alleles_unordered.length);
//        }


        // FUNCTIONAL version 2
        alleles.put("a1", allele1);
        alleles.put("a2", (Character) mapped_seg_site.values().stream().distinct().filter(x -> allele1 != x).toArray()[0]);

        // TODO:  verify code is correct by printing out mapped_seg_site and allele1_MRCAs.keys()

        return alleles;
    }


    /**
     * By definition, a homoplasy event must consist of >= 2 homoplasic mutations.  Thus, if (# MRCAs in a HashMap < 2) simply remove all MRCA nodes because this allele
     * shows 0 homoplasy events.
     */
    private void remove_non_homoplasic_mutations(HashMap<String, HashSet<String>> MRCAs) {
        if (MRCAs.size() < 2) {
//            for (Map.Entry<String, HashSet<String>>mapping : MRCAs.entrySet()) {
//                MRCAs.remove(mapping.getKey());
//            }

            MRCAs.clear();
        }
    }


    /**
     * This method reads in the phenotype file.
     *
     * control := 0
     * case := 1
     * All other pheno values are considered missing data.  User may use as many missing labels as they wish.
     */
    private HashMap<String, String> read_phenos() {
        HashMap<String, String> phenos = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(inputFiles.otherInputFiles.pheno_filename))) {
            String currline;
            while ((currline = br.readLine()) != null) {
                String[] pheno_cols = currline.split("\t");
                String strain_id = pheno_cols[0];
                String pheno = pheno_cols[1];
                phenos.put(strain_id, pheno);
            }
        } catch (Exception e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while reading the phenotype file (this is likely due to a malformed phenotype file). Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        // Sanity check:
//        System.out.println("phenos.size() = " + phenos.size());
//        System.out.println("phenos = " + phenos);

        // ORIGINAL 2 LINES THAT DO NOT WORK W FEATURE: MISSING PHENO DATA
        // following only works if 0 := control and 1 := case
        // int tot_num_cases = phenos.entrySet().stream().map(es -> Integer.parseInt(es.getValue())).mapToInt(i -> i).sum();
        // int tot_num_controls = phenos.size() - tot_num_cases;

        // filter out missing phenotypes (eg "NA") before summation of ints, ie keep all pheno data matching a "0" or "1" (control or case)
        List<String> available_phenos = phenos.values().stream().filter(s -> s.equals("0") || s.equals("1")).collect(Collectors.toList());
        int tot_num_cases = available_phenos.stream().map(s -> Integer.parseInt(s)).mapToInt(i -> i).sum();
        int tot_num_controls = available_phenos.size() - tot_num_cases;
        int tot_num_missing_phenos = phenos.size() - available_phenos.size();
//        System.out.println("tot_num_cases = " + tot_num_cases + "\ttot_num_controls = " + tot_num_controls);

        // set global variables
        tot_num_sample_cases = tot_num_cases;
        tot_num_sample_controls = tot_num_controls;

        try {
            outputOptions.log.write(String.format("%nphenotype file: cases = " + tot_num_cases + ", controls = " + tot_num_controls + ", samples missing phenotype = " + tot_num_missing_phenos + ", total sample size (cases + controls) = " + (tot_num_cases + tot_num_controls) + "%n"));
            outputOptions.log.flush();  // damn I think the %n aren't triggering flushes
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing to the log file.  Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        return phenos;
    }


    /**
     * Currently implement this method such that I can easily turn on/off a test statistic (e.g. fisher's, phyC, etc).  I plan on experimenting/comparing diff
     * association test statistics
     *
     * @param all_events each element of ArrayList represents 1 segregating site's homoplasy events for both alleles
     * @param phenos
     */
    private void assoc(ArrayList<Homoplasy_Events> all_events, HashMap<String, String> phenos) {

//        assoc_test_stat_fishers_exact(all_events, phenos);

//        assoc_test_stat_phyC(all_events, phenos);

        assoc_test_stat_binomial_test(all_events, phenos);
    }


    /**
     * test stat #1:  fisher's exact
     * <p>
     * This association test statistic only operates on those sites where there are homoplasic mutations occurring in both alleles so that a 2x2 contingency table can be
     * constructed.
     *
     * @param all_events
     * @param phenos
     */
    private void assoc_test_stat_fishers_exact(ArrayList<Homoplasy_Events> all_events, HashMap<String, String> phenos) {

        // use only "homoplasically segregating" sites
        ArrayList<Homoplasy_Events> homoplasically_segregating_sites = get_only_homoplasically_segregating_sites(all_events);
        ArrayList<Fishers_Exact_Statistic> test_statistics = fishers_exact(homoplasically_segregating_sites, phenos);

        // qvalues
        QSet qset = qvalues(test_statistics);
        output_qvalues(homoplasically_segregating_sites, test_statistics, qset);
    }


    /**
     * test stat #2:  phyC
     * <p>
     * This association test stat more closely resembles the original phyC test stat.  Specifically, the original test stat is used for point-estimates while the dn test
     * stat is used for comparing across sites to derived a maxT FWER.
     *
     * @param all_events
     * @param phenos
     */
    private void assoc_test_stat_phyC(ArrayList<Homoplasy_Events> all_events, HashMap<String, String> phenos) {

        // use only "homoplasically informative sites"
//        int min_hcount = 6;
        ArrayList<Homoplasy_Events> homoplasically_informative_sites = subset_by_min_hcount(all_events, AlgoParams.min_hcount);
        phyC_Test_Statistic[] resampled_test_statistics = phyC(homoplasically_informative_sites, phenos);

        // DEBUG
        // output_point_estimates(homoplasically_informative_sites, resampled_test_statistics);

        // DIAGNOSTIC
        diagnostic_observed_case_control_ratio_homoplasic_mutations_only(homoplasically_informative_sites);

        // qvalues
        QSet[] qset_all_alleles = qvalues_phyC(homoplasically_informative_sites, resampled_test_statistics);

        // output significance assessments (pointwise, qvalues, maxTs)
        output_significance_assessments(homoplasically_informative_sites, resampled_test_statistics, qset_all_alleles);
    }


    /**
     * test stat #3:  binomial test
     *
     * @param all_events
     * @param phenos
     */
    private void assoc_test_stat_binomial_test(ArrayList<Homoplasy_Events> all_events, HashMap<String, String> phenos) {

        // use only "homoplasically informative sites"
//        int min_hcount = 0;
        ArrayList<Homoplasy_Events> homoplasically_informative_sites = subset_by_min_hcount(all_events, AlgoParams.min_hcount);

        Binomial_Test_Stat[] resampled_test_statistics = null;
        switch (R_SPACE_TYPE) {
            case ALL_MUTATIONS:
                // RESAMPLING SPACE:  ALL MUTATIONS
//                resampled_test_statistics = resample_all_mutations_binom_test(homoplasically_informative_sites, phenos);
//                resampled_test_statistics = resample_all_mutations_binom_test_combined_nulldists(homoplasically_informative_sites, phenos);
//                resampled_test_statistics = resample_all_mutations_binom_test_combined_nulldists_memoization(homoplasically_informative_sites, phenos);
                resampled_test_statistics = resample_all_mutations_binom_test_combined_nulldists_memoization_concurrent(homoplasically_informative_sites, phenos);
                break;

            case HOMOPLASIC_MUTATIONS_ONLY:
                // RESAMPLING SPACE:  HOMOPLASIC MUTATIONS ONLY
                resampled_test_statistics = resample_homoplasic_mutations_only_binom_test(homoplasically_informative_sites, phenos);
                break;
        }

        // DIAGNOSTIC
//        diagnostic_output_resampling_derived_pvalues_only(homoplasically_informative_sites, resampled_test_statistics);

        // DIAGNOSTIC
//        diagnostic_observed_case_control_ratio_homoplasic_mutations_only(homoplasically_informative_sites);

        if (qvalues_option) {
            // qvalues
            QSet[] qset_all_alleles = qvalues_binom(homoplasically_informative_sites, resampled_test_statistics);

            // output significance assessments (pointwise, maxTs, qvalues)
            output_significance_assessments_binom(homoplasically_informative_sites, resampled_test_statistics, qset_all_alleles);
        } else {
            // output significance assessments (pointwise, maxTs, qvalues)
            output_significance_assessments_binom(homoplasically_informative_sites, resampled_test_statistics);
        }
    }


    private void output_significance_assessments_binom(ArrayList<Homoplasy_Events> homoplasically_informative_sites, Binomial_Test_Stat[] resampled_test_statistics) {

        System.out.printf(Ansi.AUTO.string("@|fg(222) %nWriting results to file . . . |@"));

        // print out 1 result row per seg site:  segsite ID, pos, a1, a2, a1 homoplasy count, a2 homoplasy count, pvalue, or, ci, qvalue, local fdr
        // all_events and test_statistics should both be in the same order (i.e. ordered by segsite_ID)
        // sanity check (both sizes must be equal)
//        System.out.println("(# homoplasically informative sites) homoplasically_informative_sites.size() = " + homoplasically_informative_sites.size());
//        System.out.println("resampled_test_statistics.length = " + resampled_test_statistics.length);

//        System.out.println(Homoplasy_Events.COL_NAMES + "\t" + Binomial_Test_Stat.COL_NAMES);  // header row
//
//        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
//            System.out.println(homoplasically_informative_sites.get(i) + "\t" + resampled_test_statistics[i]);
//        }

        try {
            // .out file
            outputOptions.out.write(Homoplasy_Events.COL_NAMES + "\t" + Binomial_Test_Stat.COL_NAMES);
            outputOptions.out.newLine();

            for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
                outputOptions.out.write(homoplasically_informative_sites.get(i) + "\t" + resampled_test_statistics[i]);
                outputOptions.out.newLine();
            }

            // a1 and a2 out files sorted by maxT
            outputOptions.out_sorted_by_a2_maxT.write(Homoplasy_Events.COL_NAMES + "\t" + Binomial_Test_Stat.COL_NAMES);
            outputOptions.out_sorted_by_a2_maxT.newLine();
            outputOptions.out_sorted_by_a1_maxT.write(Homoplasy_Events.COL_NAMES + "\t" + Binomial_Test_Stat.COL_NAMES);
            outputOptions.out_sorted_by_a1_maxT.newLine();

            // create composite objects to hold both halves of the row
            Row[] rows = new Row[homoplasically_informative_sites.size()];
            for (int i = 0; i < homoplasically_informative_sites.size(); i++)
                rows[i] = new Row(homoplasically_informative_sites.get(i), resampled_test_statistics[i]);

            // sort by a2 maxT
            // sort using lambda to create Comparator Functional Interface :)
            Arrays.parallelSort(rows, (a, b) -> {
                return Double.compare(a.test_stat_for_one_site.familywise_pvalue_a2, b.test_stat_for_one_site.familywise_pvalue_a2);
            });

            // this is elegant as well!  Loving the static method Comparator.comparing().
//            Arrays.parallelSort(rows, Comparator.comparing(Row::get_maxT));

            // write out sorted a2 maxT
            for (Row row : rows) {
                if (!Double.isNaN(row.test_stat_for_one_site.familywise_pvalue_a2)) {  // don't print out rows where a2 maxT is NaN, i.e. only a1 maxT has a value (some sites have both a1 and a2 maxT if both are homoplasically informative alleles)
                    outputOptions.out_sorted_by_a2_maxT.write(row.homoplasy_events_for_one_site.toString() + "\t" + row.test_stat_for_one_site.toString());
                    outputOptions.out_sorted_by_a2_maxT.newLine();
                }
            }

            // sort by a1 maxT
            Arrays.parallelSort(rows, (a, b) -> {
                return Double.compare(a.test_stat_for_one_site.familywise_pvalue_a1, b.test_stat_for_one_site.familywise_pvalue_a1);
            });

            // write out sorted a2 maxT
            for (Row row : rows) {
                if (!Double.isNaN(row.test_stat_for_one_site.familywise_pvalue_a1)) {  // don't print out rows where a1 maxT is NaN, i.e. only a2 maxT has a value (some sites have both a1 and a2 maxT if both are homoplasically informative alleles)
                    outputOptions.out_sorted_by_a1_maxT.write(row.homoplasy_events_for_one_site.toString() + "\t" + row.test_stat_for_one_site.toString());
                    outputOptions.out_sorted_by_a1_maxT.newLine();
                }
            }
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing results out to to file.  Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        System.out.printf("writing complete.%n");
    }


    /**
     * When one needs to see the final results without calling any of the qvalue code (only resampling-derived point-wise and family-wise estimates).
     *
     * @param homoplasically_informative_sites
     * @param resampled_test_statistics
     */
    private void diagnostic_output_resampling_derived_pvalues_only(ArrayList<Homoplasy_Events> homoplasically_informative_sites,
                                                                   Binomial_Test_Stat[] resampled_test_statistics) {

        System.out.println(Homoplasy_Events.COL_NAMES + "\t" + Binomial_Test_Stat.COL_NAMES);  // header row

        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            System.out.println(homoplasically_informative_sites.get(i) + "\t" + resampled_test_statistics[i]);
        }
    }


    private void spike_site(Homoplasy_Events spiked_site) {
        System.out.println("\nDIAGNOSTIC:  spke_site");

        //        Homoplasy_Events spiked_site = homoplasically_informative_sites.stream().filter(s -> s.segsite_ID == segsite_ID).collect(Collectors.toList()).get(0);

//        spiked_site.a2_count_extant_only = 12;  // include all internal nodes
        spiked_site.obs_counts[2] = 12;  // make the new internal nodes all cases in a2
    }


    /**
     * The base association test statistic is the binomial test (currently one-tailed test).  Resampling-derived point-wise and family-wise estimates are generated in this
     * function.
     *
     * @param homoplasically_informative_sites
     * @param phenos
     * @return
     */
    private Binomial_Test_Stat[] resample_homoplasic_mutations_only_binom_test(ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos) {

        BinomialTest binom_test = new BinomialTest();

        Binomial_Test_Stat[] resampled_test_statistics = new Binomial_Test_Stat[homoplasically_informative_sites.size()];

        // observed counts
        ArrayList<int[]> obs_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, phenos);

        // use observed counts to build resampling space
        Resampling_Space r_space = new Resampling_Space(obs_homoplasy_counts);

        // set observed values
        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {

            // set observed counts for each Homoplasy_Events obj
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            site.obs_counts = obs_homoplasy_counts.get(s);

            // DIAGNOSTIC:  spike site:  to test internal nodes as hypothetically all cases
//            if (site.segsite_ID == 581) {
//                spike_site(site);
//            }

            // calc and set observed binom test pvalues for each Binomial_Test_Stat obj
            resampled_test_statistics[s] = new Binomial_Test_Stat(site.a1_in_use, site.a2_in_use);

            if (site.a1_in_use) {
                int tot_num_homoplasic_mutations_a1 = site.obs_counts[0] + site.obs_counts[1];
                int num_cases_a1 = site.obs_counts[0];

                resampled_test_statistics[s].obs_binom_pvalue_a1 = binom_test.binomialTest(tot_num_homoplasic_mutations_a1, num_cases_a1, r_space.p_success_a1,
                        TAIL_TYPE);
            }

            if (site.a2_in_use) {
                int tot_num_homoplasic_mutations_a2 = site.obs_counts[2] + site.obs_counts[3];
                int num_cases_a2 = site.obs_counts[2];

                resampled_test_statistics[s].obs_binom_pvalue_a2 = binom_test.binomialTest(tot_num_homoplasic_mutations_a2, num_cases_a2, r_space.p_success_a2,
                        TAIL_TYPE);
            }
        }

        // DEBUG:  sanity check on EXTANT_NODES_ONLY code
//        check_extant_nodes_only_code(homoplasically_informative_sites);

        // construct both point-wise and family-wise null dists (resampling begins)
//        int m = 1000000;  // # resamples/replicates/perms
        // DEBUG
//        System.out.println("m = " + m);

        // construct a separate familywise null dist for a1 (major) and a2 (minor)
        double[] maxT_nulldist_a1 = new double[m];
        double[] maxT_nulldist_a2 = new double[m];

        // DIAGNOSTIC:  point-wise null dist by segsite_ID
//        ArrayList<int[]> pointwise_nulldist = new ArrayList<>();  // both alleles
//        int segsite_index = debug_get_index_by_segsite(75910, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(20829, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(1142, homoplasically_informative_sites);  // figure out index of this segsite

        // DIAGNOSTIC:  resampled case & control counts (all sites).  To turn on, comment out a total of 3 lines below (same comment line)
//        ArrayList<ArrayList<int[]>> all_resampled_homoplasy_counts = new ArrayList<>();  // each element represents 1 replicate across all informative sites
        Resampling_Space_Diagnostics r_space_diagnostics = new Resampling_Space_Diagnostics();

        ZonedDateTime start_of_resampling = ZonedDateTime.now();
        for (int curr_rep = 0; curr_rep < m; curr_rep++) {
            // 1 complete resampling := 3 methods:
            r_space.permute();
            ArrayList<int[]> resampled_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, phenos, r_space);
            ArrayList[] replicate = calc_binomial_test_stats(homoplasically_informative_sites, resampled_homoplasy_counts, resampled_test_statistics, binom_test, r_space);

            // DIAGNOSTIC:  point-wise null dist by segsite_ID
//            pointwise_nulldist.add(resampled_homoplasy_counts.get(segsite_index));

            // identify the replicate maxT for each allele
            ArrayList<Double> replicate_a1 = replicate[0];
            ArrayList<Double> replicate_a2 = replicate[1];
            double replicate_maxT_a1 = replicate_a1.stream().mapToDouble(d -> d).min().getAsDouble();  // there should be no Double.NaNs in this collection
            double replicate_maxT_a2 = replicate_a2.stream().mapToDouble(d -> d).min().getAsDouble();

            // allele-specific null dists (all maxT null dists should always be size m regardless of "maf" filter)
            maxT_nulldist_a1[curr_rep] = replicate_maxT_a1;
            maxT_nulldist_a2[curr_rep] = replicate_maxT_a2;

            // DIAGNOSTIC:  resampled case & control counts (all sites)
//            if (curr_rep < 100000) {  // store only the first 100,000 replicates due to memory constraints
//                all_resampled_homoplasy_counts.add(resampled_homoplasy_counts);
//            }
            r_space_diagnostics.tally_current_replicate(resampled_homoplasy_counts);

            // progress
            if ((curr_rep % 50000) == 0) {
                System.out.println("replicate: " + curr_rep + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()));
            }
        }

        calc_pointwise_estimate_binom(homoplasically_informative_sites, resampled_test_statistics, m);
        calc_familywise_estimates_binom(homoplasically_informative_sites, resampled_test_statistics, m, maxT_nulldist_a1, maxT_nulldist_a2);

        // DIAGNOSTIC
//        diagnostic_pointwise_resamplings(segsite_index, homoplasically_informative_sites, resampled_test_statistics, pointwise_nulldist);

        // DIAGNOSTIC
//        diagnostic_resampled_case_control_ratio_by_site(homoplasically_informative_sites, segsite_index, pointwise_nulldist);

        // DIAGNOSTIC:  resampled case & control counts (all sites)
//        diagnostic_resampled_pheno_counts(all_resampled_homoplasy_counts);
        r_space_diagnostics.output_resampled_pheno_counts();

        // DIAGNOSTIC:
//        double tail_cutoff = 0.007257167447129453;  // inclusive
//        diagnostic_identify_sites_contributing_to_maxT_tail(homoplasically_informative_sites, resampled_test_statistics, maxT_nulldist_a1, maxT_nulldist_a2, tail_cutoff);

        return resampled_test_statistics;
    }


    /**
     * Identifies the segsites that contribute to the right tail of a1 and/or a2 maxT null dists.
     *
     * @param homoplasically_informative_sites
     * @param resampled_test_statistics
     * @param maxT_nulldist_a1
     * @param maxT_nulldist_a2
     * @param tail_cutoff                      // function find all segsites <= this threshold (inclusive)
     */
    private void diagnostic_identify_sites_contributing_to_maxT_tail(ArrayList<Homoplasy_Events> homoplasically_informative_sites,
                                                                     Binomial_Test_Stat[] resampled_test_statistics,
                                                                     double[] maxT_nulldist_a1, double[] maxT_nulldist_a2,
                                                                     double tail_cutoff) {

        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {

        }
    }


    private ArrayList[] calc_binomial_test_stats(ArrayList<Homoplasy_Events> homoplasically_informative_sites, ArrayList<int[]> resampled_homoplasy_counts,
                                                 Binomial_Test_Stat[] resampled_test_statistics, BinomialTest binom_test, Resampling_Space r_space) {

        ArrayList<Double> binom_test_stats_a1 = new ArrayList<>();  // only stores values for alleles that are in use, thus size of ArrayList will be <= # sites
        ArrayList<Double> binom_test_stats_a2 = new ArrayList<>();

        for (int s = 0; s < homoplasically_informative_sites.size(); s++) { // for each site
            Homoplasy_Events curr_site = homoplasically_informative_sites.get(s);

            // DEBUG:  checking maxT for only site 1044
//            if (curr_site.segsite_ID != 1044) {
//                continue;
//            }
            // since a1 maxT null dist cannot be empty, give it a dummy value
//            binom_test_stats_a1.add(1d);


            Binomial_Test_Stat curr_test_stat = resampled_test_statistics[s];

            int[] curr_site_resampled_counts = resampled_homoplasy_counts.get(s);

            if (curr_site.a1_in_use) {
                // calc binom test for resampled counts
                int tot_num_trials = curr_site_resampled_counts[0] + curr_site_resampled_counts[1];  // tot # homoplasic mutations at this site for a1
                int num_cases = curr_site_resampled_counts[0];
                double resampled_binom_test_pvalue = binom_test.binomialTest(tot_num_trials, num_cases, r_space.p_success_a1, TAIL_TYPE);
                binom_test_stats_a1.add(resampled_binom_test_pvalue);

                // if resampled test stat equal or more extreme than obs test stat
                if (resampled_binom_test_pvalue <= curr_test_stat.obs_binom_pvalue_a1) {
                    curr_test_stat.r_a1++;
                }
            }

            if (curr_site.a2_in_use) {
                // calc binom test for resampled counts
                int tot_num_trials = curr_site_resampled_counts[2] + curr_site_resampled_counts[3];  // tot # homoplasic mutations at this site for a1
                int num_cases = curr_site_resampled_counts[2];
                double resampled_binom_test_pvalue = binom_test.binomialTest(tot_num_trials, num_cases, r_space.p_success_a2, TAIL_TYPE);
                binom_test_stats_a2.add(resampled_binom_test_pvalue);

                // if resampled test stat equal or more extreme than obs test stat
                if (resampled_binom_test_pvalue <= curr_test_stat.obs_binom_pvalue_a2) {
                    curr_test_stat.r_a2++;
                }
            }
        }

        // TODO BINOM:  what sites contribute to the right-tail end of maxT null dists?
        //  maybe create new data structures to hold the segsite_ID corresponding to each value in resampled binomial test pvalues, binom_test_stats_a2

        return new ArrayList[]{binom_test_stats_a1, binom_test_stats_a2};
    }


    private ArrayList<int[]> count_homoplasies(ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos, Resampling_Space r_space) {

        ArrayList<int[]> homoplasy_counts_all_sites = new ArrayList<>();

        for (Homoplasy_Events site : homoplasically_informative_sites) {
            int[] site_counts = new int[4];

            // homoplasy counts per allele (extension to tri and quad-allelic sites is straight-forward)
            int[] counts_a1 = {0, 0};  // a1 [0] := case, [1] := control
            int[] counts_a2 = {0, 0};  // a2

            if (site.a1_in_use) {
                counts_a1 = counts_by_allele(Allele_Type.A1, site.allele1_MRCAs, phenos, r_space);
            }

            if (site.a2_in_use) {
                counts_a2 = counts_by_allele(Allele_Type.A2, site.allele2_MRCAs, phenos, r_space);
            }

            // int[4] counts := a1p1, a1p2, a2p1, a2p2   p1 := case, p2 := control (same as count_phenos())
            site_counts[0] = counts_a1[0];
            site_counts[1] = counts_a1[1];
            site_counts[2] = counts_a2[0];
            site_counts[3] = counts_a2[1];

            homoplasy_counts_all_sites.add(site_counts);
        }

        return homoplasy_counts_all_sites;
    }


    /**
     * Homoplasy counts for a1 || a2 at a specific segregating site.
     * <p>
     * pheno label convention:  int[4] counts := a1p1, a1p2, a2p1, a2p2   p1 := case p2 := control (same as count_phenos() used in fisher's exact)
     *
     * @param homoplasic_mutations site.allele1_MRCAs || site.allele2_MRCAs
     * @param phenos               HashMap of all available (either all or extant only) phenos
     * @return homoplasy counts for a1 or a2.  int[2]:  [0] := p1 := case, [1] := p2 := control
     */
    private int[] counts_by_allele(Allele_Type allele_type, HashMap<String, HashSet<String>> homoplasic_mutations, HashMap<String, String> phenos, Resampling_Space r_space) {

        int[] counts_for_allele = {0, 0};

        for (Map.Entry<String, HashSet<String>> homoplasic_mutation : homoplasic_mutations.entrySet()) {
            String node_name = homoplasic_mutation.getKey();

            // if no pheno is available for this node (i.e. this is an internal node), simply continue without updating any counters
            if (phenos.containsKey(node_name)) {
//                String pheno = phenos.get(node_name);
                int pheno = r_space.next_pheno(allele_type);

                if (pheno == 1) {  // case
                    counts_for_allele[0]++;
                } else {  // control:  pheno = 0
                    counts_for_allele[1]++;
                }
            }
        }

        return counts_for_allele;
    }


    /**
     * THIS VERSION USES THE OLD RESAMPLING SPACE (IE:  ALL MUTATIONS, INCLUDING NON-HOMOPLASIC MUTATIONS)
     * <p>
     * The base association test statistic is the binomial test.  Resampling-derived point-wise and family-wise estimates are generated in this
     * function.
     *
     * @param homoplasically_informative_sites
     * @param phenos
     * @return
     */
    private Binomial_Test_Stat[] resample_all_mutations_binom_test(ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos) {

        BinomialTest binom_test = new BinomialTest();

        // calc expected background p_success values:  currently just use the p_success|both alleles value for a1 and a2 as a close approximation
        double background_p_success_given_both_alleles = calc_background_p_success();
        double p_success_a1 = background_p_success_given_both_alleles;
        double p_success_a2 = background_p_success_given_both_alleles;
        System.out.println("p_success_a1 (used in all a1 binomial tests) = " + p_success_a1);
        System.out.println("p_success_a2 (used in all a2 binomial tests) = " + p_success_a2);

        Binomial_Test_Stat[] resampled_test_statistics = new Binomial_Test_Stat[homoplasically_informative_sites.size()];

        // observed counts
        ArrayList<int[]> obs_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, phenos);

        // set observed values
        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {

            // set observed counts for each Homoplasy_Events obj
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            site.obs_counts = obs_homoplasy_counts.get(s);

            // calc and set observed binom test pvalues for each Binomial_Test_Stat obj
            resampled_test_statistics[s] = new Binomial_Test_Stat(site.a1_in_use, site.a2_in_use);

            if (site.a1_in_use) {
                int tot_num_homoplasic_mutations_a1 = site.obs_counts[0] + site.obs_counts[1];
                int num_cases_a1 = site.obs_counts[0];
                resampled_test_statistics[s].obs_binom_pvalue_a1 = binom_test.binomialTest(tot_num_homoplasic_mutations_a1, num_cases_a1, p_success_a1,
                        TAIL_TYPE);
            }

            if (site.a2_in_use) {
                int tot_num_homoplasic_mutations_a2 = site.obs_counts[2] + site.obs_counts[3];
                int num_cases_a2 = site.obs_counts[2];
                resampled_test_statistics[s].obs_binom_pvalue_a2 = binom_test.binomialTest(tot_num_homoplasic_mutations_a2, num_cases_a2, p_success_a2,
                        TAIL_TYPE);
            }
        }


        // DIAGNOSTIC
        diagnostic_observed_case_control_ratio_homoplasic_mutations_only(homoplasically_informative_sites);
        // TEMPORARY:  just need the obs p_success values, don't resample as that takes too long
//        System.exit(-1);

        // DEBUG:  sanity check on EXTANT_NODES_ONLY code
//        check_extant_nodes_only_code(homoplasically_informative_sites);

        // construct both point-wise and family-wise null dists (resampling begins)
//        int m = 1000000;  // # resamples/replicates/perms
        // DEBUG
//        System.out.println("m = " + m);

        // construct a separate familywise null dist for a1 (major) and a2 (minor)
        double[] maxT_nulldist_a1 = new double[m];
        double[] maxT_nulldist_a2 = new double[m];

        // DIAGNOSTIC:  point-wise null dist by segsite_ID
//        ArrayList<int[]> pointwise_nulldist = new ArrayList<>();  // both alleles
//        int segsite_index = debug_get_index_by_segsite(75910, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(20829, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(1142, homoplasically_informative_sites);  // figure out index of this segsite

        // DIAGNOSTIC:  resampled case & control counts
//        ArrayList<ArrayList<int[]>> all_resampled_homoplasy_counts = new ArrayList<>();  // each element represents 1 replicate across all informative sites
        Resampling_Space_Diagnostics r_space_diagnostics = new Resampling_Space_Diagnostics();

        ZonedDateTime start_of_resampling = ZonedDateTime.now();
        for (int curr_rep = 0; curr_rep < m; curr_rep++) {
            // 1 complete resampling := 3 methods:
            HashMap<String, String> permuted_phenos = permute_phenos(phenos);
            ArrayList<int[]> resampled_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, permuted_phenos);

            // DIAGNOSTIC:  output obs p_success for this replicate
//            r_space_diagnostics.output_obs_p_success(resampled_homoplasy_counts);

///*  TEMPORARY:  blocked out to quickly calc obs p_success for each repicate
            ArrayList[] replicate = calc_binomial_test_stats(homoplasically_informative_sites, resampled_homoplasy_counts, resampled_test_statistics, binom_test,
                    p_success_a1, p_success_a2);

            // DIAGNOSTIC:  point-wise null dist by segsite_ID
//            pointwise_nulldist.add(resampled_homoplasy_counts.get(segsite_index));

            // identify the replicate maxT for each allele
            ArrayList<Double> replicate_a1 = replicate[0];
            ArrayList<Double> replicate_a2 = replicate[1];
            double replicate_maxT_a1 = replicate_a1.stream().mapToDouble(d -> d).min().getAsDouble();  // there should be no Double.NaNs in this collection
            double replicate_maxT_a2 = replicate_a2.stream().mapToDouble(d -> d).min().getAsDouble();

            // allele-specific null dists (all maxT null dists should always be size m regardless of "maf" filter)
            maxT_nulldist_a1[curr_rep] = replicate_maxT_a1;
            maxT_nulldist_a2[curr_rep] = replicate_maxT_a2;

            // DIAGNOSTIC:  resampled case & control counts (all sites)
//            if (curr_rep < 100000) {  // store only the first 100,000 replicates due to memory constraints
//                all_resampled_homoplasy_counts.add(resampled_homoplasy_counts);
//            }
            r_space_diagnostics.tally_current_replicate(resampled_homoplasy_counts);
//*/

            // progress
            if ((curr_rep % 50000) == 0) {
                System.out.println("replicate: " + curr_rep + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()));
            }
        }

        //  TEMPORARY:  just output obs p_success for each replicate, avoid all other code
//        System.exit(-1);

        calc_pointwise_estimate_binom(homoplasically_informative_sites, resampled_test_statistics, m);
        calc_familywise_estimates_binom(homoplasically_informative_sites, resampled_test_statistics, m, maxT_nulldist_a1, maxT_nulldist_a2);

        // DIAGNOSTIC
//        diagnostic_pointwise_resamplings(segsite_index, homoplasically_informative_sites, resampled_test_statistics, pointwise_nulldist);

        // DIAGNOSTIC
//        diagnostic_resampled_case_control_ratio_by_site(homoplasically_informative_sites, segsite_index, pointwise_nulldist);

        // DIAGNOSTIC:  show resampled case & control counts
//        diagnostic_resampled_pheno_counts(all_resampled_homoplasy_counts);
        r_space_diagnostics.output_resampled_pheno_counts();

        return resampled_test_statistics;
    }


    /**
     * This method is the same as resample_all_mutations_binom_test() except that the family-wise null dists for a1 and a2 are combined.
     * <p>
     * THIS VERSION USES THE OLD RESAMPLING SPACE (IE:  ALL MUTATIONS, INCLUDING NON-HOMOPLASIC MUTATIONS)
     * <p>
     * The base association test statistic is the binomial test.  Resampling-derived point-wise and family-wise estimates are generated in this
     * function.
     *
     * @param homoplasically_informative_sites
     * @param phenos
     * @return
     */
    private Binomial_Test_Stat[] resample_all_mutations_binom_test_combined_nulldists(ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos) {

        BinomialTest binom_test = new BinomialTest();

        // calc expected background p_success values:  currently just use the p_success|both alleles value for a1 and a2 as a close approximation
        double background_p_success_given_both_alleles = calc_background_p_success();
        double p_success_a1 = background_p_success_given_both_alleles;
        double p_success_a2 = background_p_success_given_both_alleles;
        System.out.println("p_success_a1 (used in all a1 binomial tests) = " + p_success_a1);
        System.out.println("p_success_a2 (used in all a2 binomial tests) = " + p_success_a2);

        Binomial_Test_Stat[] resampled_test_statistics = new Binomial_Test_Stat[homoplasically_informative_sites.size()];

        // observed counts
        ArrayList<int[]> obs_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, phenos);

        // set observed values
        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {

            // set observed counts for each Homoplasy_Events obj
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            site.obs_counts = obs_homoplasy_counts.get(s);

            // calc and set observed binom test pvalues for each Binomial_Test_Stat obj
            resampled_test_statistics[s] = new Binomial_Test_Stat(site.a1_in_use, site.a2_in_use);

            if (site.a1_in_use) {
                int tot_num_homoplasic_mutations_a1 = site.obs_counts[0] + site.obs_counts[1];
                int num_cases_a1 = site.obs_counts[0];
                resampled_test_statistics[s].obs_binom_pvalue_a1 = binom_test.binomialTest(tot_num_homoplasic_mutations_a1, num_cases_a1, p_success_a1,
                        TAIL_TYPE);
            }

            if (site.a2_in_use) {
                int tot_num_homoplasic_mutations_a2 = site.obs_counts[2] + site.obs_counts[3];
                int num_cases_a2 = site.obs_counts[2];
                resampled_test_statistics[s].obs_binom_pvalue_a2 = binom_test.binomialTest(tot_num_homoplasic_mutations_a2, num_cases_a2, p_success_a2,
                        TAIL_TYPE);
            }
        }


        // DIAGNOSTIC
        diagnostic_observed_case_control_ratio_homoplasic_mutations_only(homoplasically_informative_sites);
        // TEMPORARY:  just need the obs p_success values, don't resample as that takes too long
//        System.exit(-1);

        // DEBUG:  sanity check on EXTANT_NODES_ONLY code
//        check_extant_nodes_only_code(homoplasically_informative_sites);

        // construct both point-wise and family-wise null dists (resampling begins)
//        int m = 1000000;  // # resamples/replicates/perms
        // DEBUG
//        System.out.println("m = " + m);

        // construct a separate familywise null dist for a1 (major) and a2 (minor)
//        double[] maxT_nulldist_a1 = new double[m];
//        double[] maxT_nulldist_a2 = new double[m];
        double[] maxT_nulldist_a1_a2_combined = new double[m];  // combined a1 and a2 family-wise null dist

        // DIAGNOSTIC:  point-wise null dist by segsite_ID
//        ArrayList<int[]> pointwise_nulldist = new ArrayList<>();  // both alleles
//        int segsite_index = debug_get_index_by_segsite(75910, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(20829, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(1142, homoplasically_informative_sites);  // figure out index of this segsite

        // DIAGNOSTIC:  resampled case & control counts
//        ArrayList<ArrayList<int[]>> all_resampled_homoplasy_counts = new ArrayList<>();  // each element represents 1 replicate across all informative sites
        Resampling_Space_Diagnostics r_space_diagnostics = new Resampling_Space_Diagnostics();

        ZonedDateTime start_of_resampling = ZonedDateTime.now();
        for (int curr_rep = 0; curr_rep < m; curr_rep++) {  // PARALLEL OPTIMIZATION: each replicate is a thread
            // 1 complete resampling := 3 methods:
            HashMap<String, String> permuted_phenos = permute_phenos(phenos);
            ArrayList<int[]> resampled_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, permuted_phenos);

            // DIAGNOSTIC:  output obs p_success for this replicate
//            r_space_diagnostics.output_obs_p_success(resampled_homoplasy_counts);

///*  TEMPORARY:  blocked out to quickly calc obs p_success for each repicate
            ArrayList[] replicate = calc_binomial_test_stats(homoplasically_informative_sites, resampled_homoplasy_counts, resampled_test_statistics, binom_test,
                    p_success_a1, p_success_a2);

            // DIAGNOSTIC:  point-wise null dist by segsite_ID
//            pointwise_nulldist.add(resampled_homoplasy_counts.get(segsite_index));

            // identify the replicate maxT for each allele
//            ArrayList<Double> replicate_a1 = replicate[0];
//            ArrayList<Double> replicate_a2 = replicate[1];
//            double replicate_maxT_a1 = replicate_a1.stream().mapToDouble(d -> d).min().getAsDouble();  // there should be no Double.NaNs in this collection
//            double replicate_maxT_a2 = replicate_a2.stream().mapToDouble(d -> d).min().getAsDouble();
//
//            // allele-specific null dists (all maxT null dists should always be size m regardless of "maf" filter)
//            maxT_nulldist_a1[curr_rep] = replicate_maxT_a1;
//            maxT_nulldist_a2[curr_rep] = replicate_maxT_a2;

            // COMBINED a1 and a2 FAMILY-WISE NULLDIST VERSION:
            // first combine resampled pvalues from a1 and a2 into one full repicate of resampled pvalues
            // NOTE:  I do it like this to keep more of the previous code the same, i.e. separate a1 and a2 null dists were coded first
            replicate[0].addAll(replicate[1]);
            ArrayList<Double> replicate_a1_and_a2 = replicate[0];
            // identify the replicate maxT from pooled a1 and a2 resampled pvalues
            double replicate_maxT = replicate_a1_and_a2.stream().mapToDouble(d -> d).min().getAsDouble();  // there should be no Double.NaNs in this collection
            maxT_nulldist_a1_a2_combined[curr_rep] = replicate_maxT;

            // DIAGNOSTIC:  resampled case & control counts (all sites)
//            if (curr_rep < 100000) {  // store only the first 100,000 replicates due to memory constraints
//                all_resampled_homoplasy_counts.add(resampled_homoplasy_counts);
//            }
            r_space_diagnostics.tally_current_replicate(resampled_homoplasy_counts);
//*/

            // progress
            if ((curr_rep % 50000) == 0) {
                System.out.println("replicate: " + curr_rep + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()));
            }
        }

        // OPTIMIZATION
        summarize_redundant_binom_calls();
//        System.exit(-1);

        //  TEMPORARY:  just output obs p_success for each replicate, avoid all other code
//        System.exit(-1);

        calc_pointwise_estimate_binom(homoplasically_informative_sites, resampled_test_statistics, m);
//        calc_familywise_estimates_binom(homoplasically_informative_sites, resampled_test_statistics, m, maxT_nulldist_a1, maxT_nulldist_a2);
        calc_familywise_estimates_binom_combined_nulldists(homoplasically_informative_sites, resampled_test_statistics, m, maxT_nulldist_a1_a2_combined);

        // DIAGNOSTIC
//        diagnostic_pointwise_resamplings(segsite_index, homoplasically_informative_sites, resampled_test_statistics, pointwise_nulldist);

        // DIAGNOSTIC
//        diagnostic_resampled_case_control_ratio_by_site(homoplasically_informative_sites, segsite_index, pointwise_nulldist);

        // DIAGNOSTIC:  show resampled case & control counts
//        diagnostic_resampled_pheno_counts(all_resampled_homoplasy_counts);
        r_space_diagnostics.output_resampled_pheno_counts();

        return resampled_test_statistics;
    }


    private void summarize_redundant_binom_calls() {
        // TODO:  print out total # calls, sum of all values for sanity check of total # calls, proportion unique calls, let program run till end
        // TODO:  calc time spent making ALL binom calls
        // TODO:  Tablesaw for histogram plot
//        System.out.println("\n\n");
        System.out.printf("%n%nSummary of redundant binom calls:%n");

        // total # calls calculated in 2 ways as a sanity check:
        // By HashMap:
        int tot_num_calls_storedinHashMap = 0;
        for (Integer calls : binom_test_calls.values()) {
            tot_num_calls_storedinHashMap += calls;
        }
        System.out.println("tot_num_calls_storedinHashMap = " + tot_num_calls_storedinHashMap);

        // alt way of doing += using streams
        System.out.println("tot # calls stored in HashMap using streams: " + binom_test_calls.values().stream().mapToInt(Integer::intValue).sum());

        // By int counters:
        System.out.println("tot_num_calls_by_int_counter = " + tot_num_calls_by_int_counter);

        // sanity check
        assert (tot_num_calls_storedinHashMap == tot_num_calls_by_int_counter) : "total # binom calls sanity check failed";

        // proportion of calls that are unique (i.e. same binom params of: tot # homoplasic mutations/trials, tot # cases
        double tot_num_unique_calls = binom_test_calls.size();
        System.out.println("tot_num_unique_calls = " + tot_num_unique_calls);
        double proportion_unique_calls = tot_num_unique_calls / (double) tot_num_calls_by_int_counter;
        System.out.println("proportion_unique_calls = " + proportion_unique_calls);
        System.out.println("tot_duration_of_all_binom_calls_to_math_library = " + tot_duration_of_all_binom_calls_to_math_library);
        System.out.println("unique binom calls:  " + binom_test_calls.toString());
        System.out.printf("End of summary ----------%n%n");
    }


    /**
     * This method is the same as resample_all_mutations_binom_test_combined_nulldists() except it's now optimized with memoization to guard against
     * redundant binom test calls to the Apache math library.  Hmmm, not sure if this is the method I want to branch off; where does it make sense to
     * instatiate one HashMap to be used by every binom call?
     * <p>
     * This method is the same as resample_all_mutations_binom_test() except that the family-wise null dists for a1 and a2 are combined.
     * <p>
     * THIS VERSION USES THE OLD RESAMPLING SPACE (IE:  ALL MUTATIONS, INCLUDING NON-HOMOPLASIC MUTATIONS)
     * <p>
     * The base association test statistic is the binomial test.  Resampling-derived point-wise and family-wise estimates are generated in this
     * function.
     *
     * @param homoplasically_informative_sites
     * @param phenos
     * @return
     */
    private Binomial_Test_Stat[] resample_all_mutations_binom_test_combined_nulldists_memoization(ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos) {

        BinomialTest binom_test = new BinomialTest();

        // calc expected background p_success values:  currently just use the p_success|both alleles value for a1 and a2 as a close approximation
        double background_p_success_given_both_alleles = calc_background_p_success();
        double p_success_a1 = background_p_success_given_both_alleles;
        double p_success_a2 = background_p_success_given_both_alleles;
        System.out.println("p_success_a1 (used in all a1 binomial tests) = " + p_success_a1);
        System.out.println("p_success_a2 (used in all a2 binomial tests) = " + p_success_a2);

        Binomial_Test_Stat[] resampled_test_statistics = new Binomial_Test_Stat[homoplasically_informative_sites.size()];

        // observed counts
        ArrayList<int[]> obs_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, phenos);

        // set observed values
        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {

            // set observed counts for each Homoplasy_Events obj
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            site.obs_counts = obs_homoplasy_counts.get(s);

            // calc and set observed binom test pvalues for each Binomial_Test_Stat obj
            resampled_test_statistics[s] = new Binomial_Test_Stat(site.a1_in_use, site.a2_in_use);

            if (site.a1_in_use) {
                int tot_num_homoplasic_mutations_a1 = site.obs_counts[0] + site.obs_counts[1];
                int num_cases_a1 = site.obs_counts[0];
                resampled_test_statistics[s].obs_binom_pvalue_a1 = binom_test.binomialTest(tot_num_homoplasic_mutations_a1, num_cases_a1, p_success_a1,
                        TAIL_TYPE);
            }

            if (site.a2_in_use) {
                int tot_num_homoplasic_mutations_a2 = site.obs_counts[2] + site.obs_counts[3];
                int num_cases_a2 = site.obs_counts[2];
                resampled_test_statistics[s].obs_binom_pvalue_a2 = binom_test.binomialTest(tot_num_homoplasic_mutations_a2, num_cases_a2, p_success_a2,
                        TAIL_TYPE);
            }
        }

        // DIAGNOSTIC
        diagnostic_observed_case_control_ratio_homoplasic_mutations_only(homoplasically_informative_sites);
        // TEMPORARY:  just need the obs p_success values, don't resample as that takes too long
//        System.exit(-1);

        // DEBUG:  sanity check on EXTANT_NODES_ONLY code
//        check_extant_nodes_only_code(homoplasically_informative_sites);

        // construct both point-wise and family-wise null dists (resampling begins)
//        int m = 1000000;  // # resamples/replicates/perms
        // DEBUG
//        System.out.println("m = " + m);

        // construct a separate familywise null dist for a1 (major) and a2 (minor)
//        double[] maxT_nulldist_a1 = new double[m];
//        double[] maxT_nulldist_a2 = new double[m];
        double[] maxT_nulldist_a1_a2_combined = new double[m];  // combined a1 and a2 family-wise null dist

        // DIAGNOSTIC:  point-wise null dist by segsite_ID
//        ArrayList<int[]> pointwise_nulldist = new ArrayList<>();  // both alleles
//        int segsite_index = debug_get_index_by_segsite(75910, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(20829, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(1142, homoplasically_informative_sites);  // figure out index of this segsite

        // START OF RESAMPLING ---------------------------------
        // DIAGNOSTIC:  resampled case & control counts
//        ArrayList<ArrayList<int[]>> all_resampled_homoplasy_counts = new ArrayList<>();  // each element represents 1 replicate across all informative sites
        Resampling_Space_Diagnostics r_space_diagnostics = new Resampling_Space_Diagnostics();
        HashMap<String, Double> resampled_pvalue_cache = new HashMap<>();  // cache to replace redundant binom calls
        ZonedDateTime start_of_resampling = ZonedDateTime.now();
        for (int curr_rep = 0; curr_rep < m; curr_rep++) {  // PARALLEL OPTIMIZATION: each replicate is a thread
            // 1 complete resampling := 3 methods:
            HashMap<String, String> permuted_phenos = permute_phenos(phenos);
            ArrayList<int[]> resampled_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, permuted_phenos);

            // DIAGNOSTIC:  output obs p_success for this replicate
//            r_space_diagnostics.output_obs_p_success(resampled_homoplasy_counts);

///*  TEMPORARY:  blocked out to quickly calc obs p_success for each repicate
            ArrayList[] replicate = calc_binomial_test_stats_memoization(homoplasically_informative_sites, resampled_homoplasy_counts, resampled_test_statistics, binom_test,
                    p_success_a1, p_success_a2, resampled_pvalue_cache);

            // DIAGNOSTIC:  point-wise null dist by segsite_ID
//            pointwise_nulldist.add(resampled_homoplasy_counts.get(segsite_index));

            // identify the replicate maxT for each allele
//            ArrayList<Double> replicate_a1 = replicate[0];
//            ArrayList<Double> replicate_a2 = replicate[1];
//            double replicate_maxT_a1 = replicate_a1.stream().mapToDouble(d -> d).min().getAsDouble();  // there should be no Double.NaNs in this collection
//            double replicate_maxT_a2 = replicate_a2.stream().mapToDouble(d -> d).min().getAsDouble();
//
//            // allele-specific null dists (all maxT null dists should always be size m regardless of "maf" filter)
//            maxT_nulldist_a1[curr_rep] = replicate_maxT_a1;
//            maxT_nulldist_a2[curr_rep] = replicate_maxT_a2;

            // COMBINED a1 and a2 FAMILY-WISE NULLDIST VERSION:
            // first combine resampled pvalues from a1 and a2 into one full repicate of resampled pvalues
            // NOTE:  I do it like this to keep more of the previous code the same, i.e. separate a1 and a2 null dists were coded first
            replicate[0].addAll(replicate[1]);
            ArrayList<Double> replicate_a1_and_a2 = replicate[0];
            // identify the replicate maxT from pooled a1 and a2 resampled pvalues
            double replicate_maxT = replicate_a1_and_a2.stream().mapToDouble(d -> d).min().getAsDouble();  // there should be no Double.NaNs in this collection
            maxT_nulldist_a1_a2_combined[curr_rep] = replicate_maxT;

            // DIAGNOSTIC:  resampled case & control counts (all sites)
//            if (curr_rep < 100000) {  // store only the first 100,000 replicates due to memory constraints
//                all_resampled_homoplasy_counts.add(resampled_homoplasy_counts);
//            }
            r_space_diagnostics.tally_current_replicate(resampled_homoplasy_counts);
//*/

            // progress
            if ((curr_rep % 50000) == 0) {
                System.out.println("replicate: " + curr_rep + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()));
            } else if (curr_rep == m - 1) { // last replicate
                System.out.println("replicate: " + curr_rep + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()) +
                        " <- end of resampling");
            }
        }

        // REDUNDANT BINOM CALLS
//        summarize_redundant_binom_calls();
//        System.exit(-1);

        //  TEMPORARY:  just output obs p_success for each replicate, avoid all other code
//        System.exit(-1);

        calc_pointwise_estimate_binom(homoplasically_informative_sites, resampled_test_statistics, m);
//        calc_familywise_estimates_binom(homoplasically_informative_sites, resampled_test_statistics, m, maxT_nulldist_a1, maxT_nulldist_a2);
        calc_familywise_estimates_binom_combined_nulldists(homoplasically_informative_sites, resampled_test_statistics, m, maxT_nulldist_a1_a2_combined);

        // DIAGNOSTIC
//        diagnostic_pointwise_resamplings(segsite_index, homoplasically_informative_sites, resampled_test_statistics, pointwise_nulldist);

        // DIAGNOSTIC
//        diagnostic_resampled_case_control_ratio_by_site(homoplasically_informative_sites, segsite_index, pointwise_nulldist);

        // DIAGNOSTIC:  show resampled case & control counts
//        diagnostic_resampled_pheno_counts(all_resampled_homoplasy_counts);
        r_space_diagnostics.output_resampled_pheno_counts();

        return resampled_test_statistics;
    }


    /**
     * This method is the same as resample_all_mutations_binom_test_combined_nulldists() except it's now optimized with memoization to guard against
     * redundant binom test calls to the Apache math library.  Hmmm, not sure if this is the method I want to branch off; where does it make sense to
     * instatiate one HashMap to be used by every binom call?
     * <p>
     * This method is the same as resample_all_mutations_binom_test() except that the family-wise null dists for a1 and a2 are combined.
     * <p>
     * THIS VERSION USES THE OLD RESAMPLING SPACE (IE:  ALL MUTATIONS, INCLUDING NON-HOMOPLASIC MUTATIONS)
     * <p>
     * The base association test statistic is the binomial test.  Resampling-derived point-wise and family-wise estimates are generated in this
     * function.
     *
     * @param homoplasically_informative_sites
     * @param phenos
     * @return
     */
    private Binomial_Test_Stat[] resample_all_mutations_binom_test_combined_nulldists_memoization_concurrent(ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos) {

        BinomialTest binom_test = new BinomialTest();

        // DEBUG
//        double zero_counts = binom_test.binomialTest(0, 0, 0.5, TAIL_TYPE);
//        double one_case = binom_test.binomialTest(1, 1, 0.5, TAIL_TYPE);
//        double one_control = binom_test.binomialTest(1, 0, 0.5, TAIL_TYPE);
//        System.out.println("zero_counts = " + zero_counts);
//        System.out.println("one_case = " + one_case);
//        System.out.println("one_control = " + one_control);
//        System.exit(-1);


        // calc expected background p_success values:  currently just use the p_success|both alleles value for a1 and a2 as a close approximation
        double background_p_success_given_both_alleles = calc_background_p_success();
        double p_success_a1 = background_p_success_given_both_alleles;
        double p_success_a2 = background_p_success_given_both_alleles;

        Binomial_Test_Stat[] resampled_test_statistics = new Binomial_Test_Stat[homoplasically_informative_sites.size()];

        // observed counts
        ArrayList<int[]> obs_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, phenos);

        // set observed values
        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {

            // set observed counts for each Homoplasy_Events obj
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            site.obs_counts = obs_homoplasy_counts.get(s);

            // calc and set observed binom test pvalues for each Binomial_Test_Stat obj
            resampled_test_statistics[s] = new Binomial_Test_Stat(site.a1_in_use, site.a2_in_use);

            if (site.a1_in_use) {
                int tot_num_homoplasic_mutations_a1 = site.obs_counts[0] + site.obs_counts[1];
                int num_cases_a1 = site.obs_counts[0];
                resampled_test_statistics[s].obs_binom_pvalue_a1 = binom_test.binomialTest(tot_num_homoplasic_mutations_a1, num_cases_a1, p_success_a1, TAIL_TYPE);
            }

            if (site.a2_in_use) {
                int tot_num_homoplasic_mutations_a2 = site.obs_counts[2] + site.obs_counts[3];
                int num_cases_a2 = site.obs_counts[2];
                resampled_test_statistics[s].obs_binom_pvalue_a2 = binom_test.binomialTest(tot_num_homoplasic_mutations_a2, num_cases_a2, p_success_a2, TAIL_TYPE);
            }
        }

        if (DEBUG_MODE) {
            diagnostic_observed_case_control_ratio_homoplasic_mutations_only(homoplasically_informative_sites);
        }

//        concurrecy_test();
//        System.exit(-1);

        // START OF RESAMPLING -----------------------------------
        System.out.printf(Ansi.AUTO.string("@|fg(85) %nStarting resampling . . . %n|@"));
        ExecutorService thread_pool = Executors.newFixedThreadPool(RuntimeSettings.num_threads);
        CountDownLatch count_down_latch = new CountDownLatch(m);
        System.out.printf("New thread pool created.  Status: " + get_thread_pool_status(thread_pool) + "%n");

        if (DEBUG_MODE) {
            try {
                outputOptions.debug.write("\nresample_all_mutations_binom_test_combined_nulldists_memoization_concurrent():\n");
                int max_threads = Runtime.getRuntime().availableProcessors();
                outputOptions.debug.write("max_threads = " + max_threads + "\n"); // max logical processors
                outputOptions.debug.write("RuntimeSettings.num_threads = " + RuntimeSettings.num_threads + "\n");  // testing user-specified # threads and default value

                // did thread pool get set to user-specified num_thread or default value?
                ThreadPoolExecutor thread_pool_executor = (ThreadPoolExecutor) thread_pool;  // needed an explicit typecast to get to the concrete implementation of ExecutorService so i can access the method below!
                outputOptions.debug.write("\nBefore submitting replicates to thread pool:\n");
                outputOptions.debug.write("thread_pool_executor.getMaximumPoolSize() = " + thread_pool_executor.getMaximumPoolSize() + "\n");
                outputOptions.debug.write("thread_pool status: " + thread_pool.toString() + "\n");
                outputOptions.debug.write("count_down_latch status: " + count_down_latch.toString() + "\n");
                outputOptions.debug.flush();
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing to the debug file.  Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }

        double[] maxT_nulldist_a1_a2_combined = new double[m];  // combined a1 and a2 family-wise null dist
        Resampling_Space_Diagnostics r_space_diagnostics = new Resampling_Space_Diagnostics();
        ConcurrentHashMap<String, Double> resampled_pvalue_cache = new ConcurrentHashMap<>();  // thread-safe cache to replace redundant binom calls, curcially put() is atomic while get() is not synchronized
        ZonedDateTime start_of_resampling = ZonedDateTime.now();
        for (int curr_replicate = 0; curr_replicate < m; curr_replicate++) {  // PARALLEL OPTIMIZATION: each replicate is a thread
            thread_pool.submit(new Replicate(curr_replicate, count_down_latch,
                                             homoplasically_informative_sites, phenos,
                                             maxT_nulldist_a1_a2_combined, r_space_diagnostics, resampled_pvalue_cache, start_of_resampling,
                                             binom_test, p_success_a1, p_success_a2, resampled_test_statistics));

            /*// 1 complete resampling := 3 methods:
            HashMap<String, String> permuted_phenos = permute_phenos(phenos);
            ArrayList<int[]> resampled_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, permuted_phenos);
            ArrayList[] replicate = calc_binomial_test_stats_memoization(homoplasically_informative_sites, resampled_homoplasy_counts, resampled_test_statistics, binom_test,
                                                                         p_success_a1, p_success_a2, resampled_pvalue_cache);
            // COMBINED a1 and a2 FAMILY-WISE NULLDIST VERSION:
            // first combine resampled pvalues from a1 and a2 into one full repicate of resampled pvalues
            // NOTE:  I do it like this to keep more of the previous code the same, i.e. separate a1 and a2 null dists were coded first
            replicate[0].addAll(replicate[1]);
            ArrayList<Double> replicate_a1_and_a2 = replicate[0];
            // identify the replicate maxT from pooled a1 and a2 resampled pvalues
            double replicate_maxT = replicate_a1_and_a2.stream().mapToDouble(d -> d).min().getAsDouble();  // there should be no Double.NaNs in this collection
            maxT_nulldist_a1_a2_combined[curr_rep] = replicate_maxT;

            r_space_diagnostics.tally_current_replicate(resampled_homoplasy_counts);

            // progress
            if ((curr_rep % 50000) == 0) {
                System.out.println("replicate: " + curr_rep + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()));
            } else if (curr_rep == m - 1) { // last replicate
                System.out.println("replicate: " + curr_rep + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()) +
                        " <- end of resampling");
            }*/
        }  // END OF RESAMPLING ----------------------------------
        try {
//            System.out.println("count_down_latch.await() BEFORE call");
            count_down_latch.await();
//            System.out.println("count_down_latch.await() AFTER call");  // are all the replicates complete when this call is reached?
        } catch (InterruptedException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while working with the count down latch for multi-threading.  Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }


        if (DEBUG_MODE) {
            try {
                outputOptions.debug.write("\nCount down complete:\n");
                outputOptions.debug.write("count_down_latch.getCount() = " + count_down_latch.getCount() + "\n");
                outputOptions.debug.write("thread_pool status: " + thread_pool.toString() + "\n");
                outputOptions.debug.flush();
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing to the debug file.  Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }

        thread_pool.shutdown();

        if (DEBUG_MODE) {
            try {
                outputOptions.debug.write("thread_pool is shutting down, current status: " + thread_pool.toString() + "\n");
                outputOptions.debug.flush();
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing to the debug file.  Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }

        System.out.printf("All replicates processed, thread pool is shutting down.  Status: " + get_thread_pool_status(thread_pool) + "%n");
        System.out.printf("Resampling complete.%n");

        // REDUNDANT BINOM CALLS
//        summarize_redundant_binom_calls();
//        System.exit(-1);
        calc_pointwise_estimate_binom(homoplasically_informative_sites, resampled_test_statistics, m);
        calc_familywise_estimates_binom_combined_nulldists(homoplasically_informative_sites, resampled_test_statistics, m, maxT_nulldist_a1_a2_combined);
        // DIAGNOSTIC
//        diagnostic_pointwise_resamplings(segsite_index, homoplasically_informative_sites, resampled_test_statistics, pointwise_nulldist);

        // DIAGNOSTIC
//        diagnostic_resampled_case_control_ratio_by_site(homoplasically_informative_sites, segsite_index, pointwise_nulldist);

        // DIAGNOSTIC:  show resampled case & control counts
//        diagnostic_resampled_pheno_counts(all_resampled_homoplasy_counts);
        if (DEBUG_MODE) {
            r_space_diagnostics.output_resampled_pheno_counts();
        }

        return resampled_test_statistics;
    }



    private String get_thread_pool_status(ExecutorService thread_pool) {
        String raw = thread_pool.toString();
        return raw.substring(raw.indexOf("["));
    }


    private class Replicate implements Runnable {

        private int replicate_num;
        private CountDownLatch count_down_latch;
        private ArrayList<Homoplasy_Events> homoplasically_informative_sites;
        private HashMap<String, String> phenos;
        private double[] maxT_nulldist_a1_a2_combined;
        private Resampling_Space_Diagnostics r_space_diagnostics;
        private ConcurrentHashMap<String, Double> resampled_pvalue_cache;
        private ZonedDateTime start_of_resampling;
        private BinomialTest binom_test;
        private double p_success_a1;
        private double p_success_a2;
        private Binomial_Test_Stat[] resampled_test_statistics;

        public Replicate(int replicate_num, CountDownLatch count_down_latch,
                         ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos,
                         double[] maxT_nulldist_a1_a2_combined, Resampling_Space_Diagnostics r_space_diagnostics, ConcurrentHashMap<String, Double> resampled_pvalue_cache, ZonedDateTime start_of_resampling,
                         BinomialTest binom_test, double p_success_a1, double p_success_a2, Binomial_Test_Stat[] resampled_test_statistics) {
            this.replicate_num = replicate_num;
            this.count_down_latch = count_down_latch;
            this.homoplasically_informative_sites = homoplasically_informative_sites;
            this.phenos = phenos;
            this.maxT_nulldist_a1_a2_combined = maxT_nulldist_a1_a2_combined;
            this.r_space_diagnostics = r_space_diagnostics;
            this.resampled_pvalue_cache = resampled_pvalue_cache;
            this.start_of_resampling = start_of_resampling;
            this.binom_test = binom_test;
            this.p_success_a1 = p_success_a1;
            this.p_success_a2 = p_success_a2;
            this.resampled_test_statistics = resampled_test_statistics;
        }

        @Override
        public void run() {
            try {
                HashMap<String, String> permuted_phenos = permute_phenos(phenos);
                ArrayList<int[]> resampled_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, permuted_phenos);
                ArrayList[] replicate = calc_binomial_test_stats_memoization_concurrent(homoplasically_informative_sites, resampled_homoplasy_counts, resampled_test_statistics, binom_test,
                        p_success_a1, p_success_a2, resampled_pvalue_cache);
                // COMBINED a1 and a2 FAMILY-WISE NULLDIST VERSION:
                // first combine resampled pvalues from a1 and a2 into one full repicate of resampled pvalues
                // NOTE:  I do it like this to keep more of the previous code the same, i.e. separate a1 and a2 null dists were coded first
                replicate[0].addAll(replicate[1]);
                ArrayList<Double> replicate_a1_and_a2 = replicate[0];
                // identify the replicate maxT from pooled a1 and a2 resampled pvalues
                double replicate_maxT = replicate_a1_and_a2.stream().mapToDouble(d -> d).min().getAsDouble();  // there should be no Double.NaNs in this collection
                maxT_nulldist_a1_a2_combined[replicate_num] = replicate_maxT;

                r_space_diagnostics.tally_current_replicate(resampled_homoplasy_counts);

                // progress indicator for resampling:
                if (DEBUG_MODE) {
                    if ((replicate_num % 10000) == 0)
                        System.out.println("replicate: " + replicate_num + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()) + "  count_down_latch.getCount() = " + count_down_latch.getCount());
                    else if (replicate_num == m - 1) // last replicate
                        System.out.println("replicate: " + replicate_num + "     time elapsed since start of resampling: " + Duration.between(start_of_resampling, ZonedDateTime.now()) + "  count_down_latch.getCount() = " + count_down_latch.getCount() + " <- end of resampling");
                } else {  // else not in debug mode, thus user does not need as much info about the thread pool and count down latch. Also make replicate count 1-based (fun with modulo math!)
                    if (((replicate_num + 1) % 10000) == 0) {
                        Duration diff = Duration.between(start_of_resampling, ZonedDateTime.now());
                        String diff_pretty = String.format("%dh %02dm %02ds", diff.toHours(), diff.toMinutesPart(), diff.toSecondsPart());
                        System.out.printf("replicate: " + (replicate_num + 1) + "     (time elapsed since start of resampling: " + diff_pretty + ")%n");
                    } else if ((replicate_num + 1) == m) { // last replicate
                        Duration diff = Duration.between(start_of_resampling, ZonedDateTime.now());
                        String diff_pretty = String.format("%dh %02dm %02ds", diff.toHours(), diff.toMinutesPart(), diff.toSecondsPart());
                        System.out.printf("replicate: " + (replicate_num + 1) + "     (time elapsed since start of resampling: " + diff_pretty + ")%n");
                    }
                }

//            System.out.println("replicate_num = " + replicate_num);
                count_down_latch.countDown();
            } catch (Exception e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred in a thread during resampling.  Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }
    }


    /**
     * This method tests the concurrency architecture I've designed that houses the resampling code.  Namely, I'm using a fixed thread pool
     * with a CountDownLatch synchronizer.
     */
    private void concurrecy_test() {
        int tot_num_replicates = 50;
        int max_threads = Runtime.getRuntime().availableProcessors();
        System.out.printf("%n%nmax_threads = " + max_threads + "%n");
        ExecutorService thread_pool = Executors.newFixedThreadPool(max_threads);
        CountDownLatch count_down_latch = new CountDownLatch(tot_num_replicates);

        System.out.println("thread_pool status: " + thread_pool.toString());
        System.out.println("count_down_latch status: " + count_down_latch.toString());

        for (int curr_replicate = 0; curr_replicate < tot_num_replicates; curr_replicate++) {
            thread_pool.submit(new Replicate_Test(curr_replicate, count_down_latch));
        }

        try {
            System.out.println("count_down_latch.await() BEFORE call");
            count_down_latch.await();
            System.out.println("count_down_latch.await() AFTER call");  // are all the replicate complete when this call is reached?
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        System.out.println("count down complete, count_down_latch.getCount() = " + count_down_latch.getCount());
        System.out.println("thread_pool status: " + thread_pool.toString());

        // probably should shutdown thread pool
        thread_pool.shutdown();
        System.out.println("thread_pool is shut down, current status: " + thread_pool.toString());
    }


    /**
     * Used to test concurrency code. Each Replicate_Test represents one thread.
     */
    private class Replicate_Test implements Runnable {

        private int replicate_num;
        private CountDownLatch count_down_latch;

        public Replicate_Test(int replicate_num, CountDownLatch count_down_latch) {
            this.replicate_num = replicate_num;
            this.count_down_latch = count_down_latch;
        }

        @Override
        public void run() {
            System.out.println("replicate_num = " + replicate_num);
            count_down_latch.countDown();
        }
    }


    /**
     * @return p_success | both alleles (based upon overall sample case:control)
     */
    private double calc_background_p_success() {

        double background_p_success_given_both_alleles = (double) tot_num_sample_cases / (double) (tot_num_sample_cases + tot_num_sample_controls);

        if (DEBUG_MODE) {
            try {
                outputOptions.debug.write("\ncalc_background_p_success():\n");
                outputOptions.debug.write("background_p_success_given_both_alleles = " + background_p_success_given_both_alleles + "\n");
                outputOptions.debug.flush();
//                System.out.println("p_success_a2 (used in all a2 binomial tests) = " + p_success_a2);
//                System.out.println("p_success_a1 (used in all a1 binomial tests) = " + p_success_a1);
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the debug file.  Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }

        return background_p_success_given_both_alleles;
    }


    /**
     * Resampled case & control counts over only homoplasic mutations and over all homoplasically informative sites.
     *
     * @param all_resampled_homoplasy_counts
     */
    private void diagnostic_resampled_pheno_counts(ArrayList<ArrayList<int[]>> all_resampled_homoplasy_counts) {
        System.out.printf("%n%nDIAGNOSTIC:  Resampled case and control counts over all replicates%n");

        int tot_resampled_cases_a1 = 0;
        int tot_resampled_controls_a1 = 0;
        int tot_resampled_cases_a2 = 0;
        int tot_resampled_controls_a2 = 0;

        for (ArrayList<int[]> replicate : all_resampled_homoplasy_counts) {
            for (int[] site : replicate) {
                tot_resampled_cases_a1 += site[0];
                tot_resampled_controls_a1 += site[1];
                tot_resampled_cases_a2 += site[2];
                tot_resampled_controls_a2 += site[3];
            }
        }

        System.out.println("tot_resampled_cases_a1 = " + tot_resampled_cases_a1);
        System.out.println("tot_resampled_controls_a1 = " + tot_resampled_controls_a1);
        System.out.println("tot_resampled_cases_a2 = " + tot_resampled_cases_a2);
        System.out.println("tot_resampled_controls_a2 = " + tot_resampled_controls_a2);
        System.out.println("tot_resampled_a1_counts = " + (tot_resampled_cases_a1 + tot_resampled_controls_a1));
        System.out.println("tot_resampled_a2_counts = " + (tot_resampled_cases_a2 + tot_resampled_controls_a2));

        float resampled_case_control_ratio_a1 = (float) tot_resampled_cases_a1 / (float) tot_resampled_controls_a1;
        float resampled_case_control_ratio_a2 = (float) tot_resampled_cases_a2 / (float) tot_resampled_controls_a2;
        System.out.println("resampled_case_control_ratio_a1 = " + resampled_case_control_ratio_a1);
        System.out.println("resampled_case_control_ratio_a2 = " + resampled_case_control_ratio_a2);

        float p_success_both_alleles = (float) (tot_resampled_cases_a1 + tot_resampled_cases_a2) / (float) (tot_resampled_cases_a1 + tot_resampled_cases_a2 + tot_resampled_controls_a1 + tot_resampled_controls_a2);
        float p_success_a1 = (float) tot_resampled_cases_a1 / (float) (tot_resampled_cases_a1 + tot_resampled_controls_a1);
        float p_success_a2 = (float) tot_resampled_cases_a2 / (float) (tot_resampled_cases_a2 + tot_resampled_controls_a2);
        System.out.println("Resampled p(success|both alleles) := # cases / (# cases + # controls) = " + p_success_both_alleles);
        System.out.println("Resampled p(success|a1 only) = " + p_success_a1);
        System.out.println("Resampled p(success|a2 only) = " + p_success_a2);
        System.out.println();
    }


    /**
     * Observed case & control counts over only homoplasic mutations.
     *
     * @param homoplasically_informative_sites
     */
    private void diagnostic_observed_case_control_ratio_homoplasic_mutations_only(ArrayList<Homoplasy_Events> homoplasically_informative_sites) {
        try {
            outputOptions.debug.write("\ndiagnostic_observed_case_control_ratio_homoplasic_mutations_only():\n");
            outputOptions.debug.write("Observed case and control counts over only homoplasic mutations:\n");

            int tot_obs_cases_a1 = 0;
            int tot_obs_controls_a1 = 0;
            int tot_obs_cases_a2 = 0;
            int tot_obs_controls_a2 = 0;

            // NOTE:  noticed I'm not using a1.is_in_use boolean, thus this loop adds up ALL cases and controls regardless if the allele passes the min_hcount threshold
            for (Homoplasy_Events site : homoplasically_informative_sites) {
                tot_obs_cases_a1 += site.obs_counts[0];
                tot_obs_controls_a1 += site.obs_counts[1];

                tot_obs_cases_a2 += site.obs_counts[2];
                tot_obs_controls_a2 += site.obs_counts[3];
            }

//        System.out.println("p_success_a1 = " + p_success_a1);
//        System.out.println("p_success_a2 = " + p_success_a2);

            outputOptions.debug.write("tot_obs_cases_a1 = " + tot_obs_cases_a1 + "\n");
            outputOptions.debug.write("tot_obs_controls_a1 = " + tot_obs_controls_a1 + "\n");
            outputOptions.debug.write("tot_obs_cases_a2 = " + tot_obs_cases_a2 + "\n");
            outputOptions.debug.write("tot_obs_controls_a2 = " + tot_obs_controls_a2 + "\n");
            outputOptions.debug.write("tot_observed_a1_counts = " + (tot_obs_cases_a1 + tot_obs_controls_a1) + "\n");
            outputOptions.debug.write("tot_observed_a2_counts = " + (tot_obs_cases_a2 + tot_obs_controls_a2) + "\n");

            float obs_homoplasic_mutations_case_control_ratio_a1 = (float) tot_obs_cases_a1 / (float) tot_obs_controls_a1;
            float obs_homoplasic_mutations_case_control_ratio_a2 = (float) tot_obs_cases_a2 / (float) tot_obs_controls_a2;
            outputOptions.debug.write("obs_homoplasic_mutations_case_control_ratio_a1 = " + obs_homoplasic_mutations_case_control_ratio_a1 + "\n");
            outputOptions.debug.write("obs_homoplasic_mutations_case_control_ratio_a2 = " + obs_homoplasic_mutations_case_control_ratio_a2 + "\n");

            float p_success_both_alleles = (float) (tot_obs_cases_a1 + tot_obs_cases_a2) / (float) (tot_obs_cases_a1 + tot_obs_cases_a2 + tot_obs_controls_a1 + tot_obs_controls_a2);
            float p_success_a1 = (float) tot_obs_cases_a1 / (float) (tot_obs_cases_a1 + tot_obs_controls_a1);
            float p_success_a2 = (float) tot_obs_cases_a2 / (float) (tot_obs_cases_a2 + tot_obs_controls_a2);
            outputOptions.debug.write("Observed p(success|both alleles) := # cases / (# cases + # controls) = " + p_success_both_alleles + "\n");
            outputOptions.debug.write("Observed p(success|a1 only) = " + p_success_a1 + "\n");
            outputOptions.debug.write("Observed p(success|a2 only) = " + p_success_a2 + "\n");
            outputOptions.debug.flush();
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the debug file.  Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }
    }


    /**
     * This method will likely only be used for debugging purposes since the end-user is interested in the multi-hypo corrected significance values.
     *
     * @param homoplasically_informative_sites
     * @param resampled_test_statistics
     */
    private void output_point_estimates(ArrayList<Homoplasy_Events> homoplasically_informative_sites, phyC_Test_Statistic[] resampled_test_statistics) {

        System.out.println(Homoplasy_Events.COL_NAMES + "\t" + phyC_Test_Statistic.COL_NAMES);  // header row

        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            System.out.println(homoplasically_informative_sites.get(i) + "\t" + resampled_test_statistics[i]);
        }
    }


    /**
     * Both alleles' qvalues are outputted across all sites.
     *
     * @param homoplasically_informative_sites
     * @param resampled_test_statistics
     * @param qset_all_alleles
     */
    private void output_significance_assessments(ArrayList<Homoplasy_Events> homoplasically_informative_sites, phyC_Test_Statistic[] resampled_test_statistics,
                                                 QSet[] qset_all_alleles) {

        // print out 1 result row per seg site:  segsite ID, pos, a1, a2, a1 homoplasy count, a2 homoplasy count, pvalue, or, ci, qvalue, local fdr
        // all_events and test_statistics should both be in the same order (i.e. ordered by segsite_ID)
        // sanity check (both sizes must be equal)
//        System.out.println("(# homoplasically informative sites) homoplasically_informative_sites.size() = " + homoplasically_informative_sites.size());
        System.out.println("resampled_test_statistics.length = " + resampled_test_statistics.length);

        ArrayList<String> qvalues_a1 = qset_all_alleles[0].qvalues;
        ArrayList<String> qvalues_a2 = qset_all_alleles[1].qvalues;
        ArrayList<String> lfdrs_a1 = qset_all_alleles[0].local_fdrs;
        ArrayList<String> lfdrs_a2 = qset_all_alleles[1].local_fdrs;

        System.out.println(Homoplasy_Events.COL_NAMES + "\t" + phyC_Test_Statistic.COL_NAMES + "\tqvalue_a1\tqvalue_a2\tlfdr_a1\tlfdr_a2");  // header row

        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            System.out.println(homoplasically_informative_sites.get(i) + "\t" +
                    resampled_test_statistics[i] + "\t" +
                    qvalues_a1.get(i) + "\t" + qvalues_a2.get(i) + "\t" +
                    lfdrs_a1.get(i) + "\t" + lfdrs_a2.get(i));
        }
    }


    private void output_significance_assessments_binom(ArrayList<Homoplasy_Events> homoplasically_informative_sites, Binomial_Test_Stat[] resampled_test_statistics,
                                                       QSet[] qset_all_alleles) {

        // print out 1 result row per seg site:  segsite ID, pos, a1, a2, a1 homoplasy count, a2 homoplasy count, pvalue, or, ci, qvalue, local fdr
        // all_events and test_statistics should both be in the same order (i.e. ordered by segsite_ID)
        // sanity check (both sizes must be equal)
//        System.out.println("(# homoplasically informative sites) homoplasically_informative_sites.size() = " + homoplasically_informative_sites.size());
        System.out.println("resampled_test_statistics.length = " + resampled_test_statistics.length);

        ArrayList<String> qvalues_a1 = qset_all_alleles[0].qvalues;
        ArrayList<String> qvalues_a2 = qset_all_alleles[1].qvalues;
        ArrayList<String> lfdrs_a1 = qset_all_alleles[0].local_fdrs;
        ArrayList<String> lfdrs_a2 = qset_all_alleles[1].local_fdrs;

        System.out.println(Homoplasy_Events.COL_NAMES + "\t" + Binomial_Test_Stat.COL_NAMES + "\tqvalue_a1\tqvalue_a2\tlfdr_a1\tlfdr_a2");  // header row

        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            System.out.println(homoplasically_informative_sites.get(i) + "\t" +
                    resampled_test_statistics[i] + "\t" +
                    qvalues_a1.get(i) + "\t" + qvalues_a2.get(i) + "\t" +
                    lfdrs_a1.get(i) + "\t" + lfdrs_a2.get(i));
        }
    }


    /**
     * I could use one qvalues() and pass in pvalues directly (so phyC would call qvalues() twice, one for each allele).  But, output_qvalues() would still require a new function
     * simply because of the 2 alleles.  So for now, I've opted to simply write a new function for qvalues() for phyC.  Will think about design again after debugging 1.0.
     *
     * @param homoplasically_informative_sites
     * @param resampled_test_statistics
     * @return
     */
    private QSet[] qvalues_phyC(ArrayList<Homoplasy_Events> homoplasically_informative_sites, phyC_Test_Statistic[] resampled_test_statistics) {

        // Keeping site-centric code:
        // get pvalues for all alleles:
        double[] pvalues_a1 = new double[resampled_test_statistics.length];
        for (int i = 0; i < pvalues_a1.length; i++) {
            pvalues_a1[i] = resampled_test_statistics[i].pointwise_pvalue_a1;
        }

        double[] pvalues_a2 = new double[resampled_test_statistics.length];
        for (int i = 0; i < pvalues_a2.length; i++) {
            pvalues_a2[i] = resampled_test_statistics[i].pointwise_pvalue_a2;
        }

        QSet qset_a1 = qvalues_R(pvalues_a1);
        QSet qset_a2 = qvalues_R(pvalues_a2);

/*
        // DEBUG:  is R code processing NaN values as NAs?
        long num_NANs_in_pvalues_a1 = Arrays.stream(pvalues_a1).filter(p -> Double.isNaN(p)).count();
        long num_NANs_in_pvalues_a2 = Arrays.stream(pvalues_a2).filter(p -> Double.isNaN(p)).count();
        System.out.println("num_NANs_in_pvalues_a1 = " + num_NANs_in_pvalues_a1);
        System.out.println("num_NANs_in_pvalues_a2 = " + num_NANs_in_pvalues_a2);
        System.out.println("pvalues_a1.length = " + pvalues_a1.length);
        System.out.println("qset_a1.qvalues.size() = " + qset_a1.qvalues.size());
        System.out.println("qset_a1.local_fdrs.size() = " + qset_a1.local_fdrs.size());
        System.out.println("pvalues_a2.length = " + pvalues_a2.length);
        System.out.println("qset_a2.qvalues.size() = " + qset_a2.qvalues.size());
        System.out.println("qset_a2.local_fdrs.size() = " + qset_a2.local_fdrs.size());
        System.out.println(String.valueOf(Double.NaN));
*/

/*
        // ALLELE-CENTRIC
        ArrayList<Double> pvalues_a1 = new ArrayList<>();
        ArrayList<Double> pvalues_a2 = new ArrayList<>();
        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(i);

            if (site.a1_in_use) {
                pvalues_a1.add(resampled_test_statistics[i].pointwise_pvalue_a1);
            }

            if (site.a2_in_use) {
                pvalues_a2.add(resampled_test_statistics[i].pointwise_pvalue_a2);
            }
        }

        // old site-centric code used double[] for pvalues, so here I'll continue with that (streams and auto-unboxing to the rescue) so I don't modify qvalues_R().
        // However, coding qvalues_R(ArrayList pvalues) should be a minor modification.
        QSet qset_a1 = qvalues_R(pvalues_a1.stream().mapToDouble(d -> d).toArray());
        QSet qset_a2 = qvalues_R(pvalues_a2.stream().mapToDouble(d -> d).toArray());
*/
/*
        // DEBUG:  sizes of pvalues should equal # informative alleles
        System.out.println("pvalues_a1.size() = " + pvalues_a1.size());
        System.out.println("pvalues_a2.size() = " + pvalues_a2.size());
        System.out.println("tot # pvalues = " + (pvalues_a1.size() + pvalues_a2.size()));
*//*


        // DEBUG:  missing values problem
        System.out.println("pvalues_a1.size() = " + pvalues_a1.size());
        System.out.println("qset_a1.qvalues.size() = " + qset_a1.qvalues.size());
        System.out.println("qset_a1.local_fdrs.size() = " + qset_a1.local_fdrs.size());
        System.out.println("pvalues_a2.size() = " + pvalues_a2.size());
        System.out.println("qset_a2.qvalues.size() = " + qset_a2.qvalues.size());
        System.out.println("qset_a2.local_fdrs.size() = " + qset_a2.local_fdrs.size());
*/

        return new QSet[]{qset_a1, qset_a2};

    }


    private QSet[] qvalues_binom(ArrayList<Homoplasy_Events> homoplasically_informative_sites, Binomial_Test_Stat[] resampled_test_statistics) {

        // Keeping site-centric code:
        // get pvalues for all alleles:
        double[] pvalues_a1 = new double[resampled_test_statistics.length];
        for (int i = 0; i < pvalues_a1.length; i++) {
            pvalues_a1[i] = resampled_test_statistics[i].pointwise_pvalue_a1;
        }

        double[] pvalues_a2 = new double[resampled_test_statistics.length];
        for (int i = 0; i < pvalues_a2.length; i++) {
            pvalues_a2[i] = resampled_test_statistics[i].pointwise_pvalue_a2;
        }

        QSet qset_a1 = qvalues_R(pvalues_a1);
        QSet qset_a2 = qvalues_R(pvalues_a2);

/*
        // DEBUG:  is R code processing NaN values as NAs?
        long num_NANs_in_pvalues_a1 = Arrays.stream(pvalues_a1).filter(p -> Double.isNaN(p)).count();
        long num_NANs_in_pvalues_a2 = Arrays.stream(pvalues_a2).filter(p -> Double.isNaN(p)).count();
        System.out.println("num_NANs_in_pvalues_a1 = " + num_NANs_in_pvalues_a1);
        System.out.println("num_NANs_in_pvalues_a2 = " + num_NANs_in_pvalues_a2);
        System.out.println("pvalues_a1.length = " + pvalues_a1.length);
        System.out.println("qset_a1.qvalues.size() = " + qset_a1.qvalues.size());
        System.out.println("qset_a1.local_fdrs.size() = " + qset_a1.local_fdrs.size());
        System.out.println("pvalues_a2.length = " + pvalues_a2.length);
        System.out.println("qset_a2.qvalues.size() = " + qset_a2.qvalues.size());
        System.out.println("qset_a2.local_fdrs.size() = " + qset_a2.local_fdrs.size());
        System.out.println(String.valueOf(Double.NaN));
*/

/*
        // ALLELE-CENTRIC
        ArrayList<Double> pvalues_a1 = new ArrayList<>();
        ArrayList<Double> pvalues_a2 = new ArrayList<>();
        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(i);

            if (site.a1_in_use) {
                pvalues_a1.add(resampled_test_statistics[i].pointwise_pvalue_a1);
            }

            if (site.a2_in_use) {
                pvalues_a2.add(resampled_test_statistics[i].pointwise_pvalue_a2);
            }
        }

        // old site-centric code used double[] for pvalues, so here I'll continue with that (streams and auto-unboxing to the rescue) so I don't modify qvalues_R().
        // However, coding qvalues_R(ArrayList pvalues) should be a minor modification.
        QSet qset_a1 = qvalues_R(pvalues_a1.stream().mapToDouble(d -> d).toArray());
        QSet qset_a2 = qvalues_R(pvalues_a2.stream().mapToDouble(d -> d).toArray());
*/
/*
        // DEBUG:  sizes of pvalues should equal # informative alleles
        System.out.println("pvalues_a1.size() = " + pvalues_a1.size());
        System.out.println("pvalues_a2.size() = " + pvalues_a2.size());
        System.out.println("tot # pvalues = " + (pvalues_a1.size() + pvalues_a2.size()));
*//*


        // DEBUG:  missing values problem
        System.out.println("pvalues_a1.size() = " + pvalues_a1.size());
        System.out.println("qset_a1.qvalues.size() = " + qset_a1.qvalues.size());
        System.out.println("qset_a1.local_fdrs.size() = " + qset_a1.local_fdrs.size());
        System.out.println("pvalues_a2.size() = " + pvalues_a2.size());
        System.out.println("qset_a2.qvalues.size() = " + qset_a2.qvalues.size());
        System.out.println("qset_a2.local_fdrs.size() = " + qset_a2.local_fdrs.size());
*/

        return new QSet[]{qset_a1, qset_a2};

    }


    /**
     * // TODO:  refactor this method so that all test stats use this method
     * <p>
     * Counts the geno/pheno pairings for both alleles.
     * Currently, only considering case/control phenotypes (2 phenos).
     *
     * @param homoplasically_informative_sites
     * @param phenos
     * @return
     */
    private ArrayList<int[]> count_homoplasies(ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos) {
        ArrayList<int[]> homoplasy_counts_all_sites = new ArrayList<>();

        for (Homoplasy_Events site : homoplasically_informative_sites) {
            int[] site_counts = new int[4];

            // homoplasy counts per allele (extension to tri and quad-allelic sites is straight-forward)
            int[] counts_a1 = {0, 0};  // a1 [0] := case, [1] := control
            int[] counts_a2 = {0, 0};  // a2

            if (site.a1_in_use) {
                counts_a1 = counts_by_allele(site.allele1_MRCAs, phenos);
            }

            if (site.a2_in_use) {
                counts_a2 = counts_by_allele(site.allele2_MRCAs, phenos);
            }

            // int[4] counts := a1p1, a1p2, a2p1, a2p2   p1 := case, p2 := control (same as count_phenos())
            site_counts[0] = counts_a1[0];
            site_counts[1] = counts_a1[1];
            site_counts[2] = counts_a2[0];
            site_counts[3] = counts_a2[1];

            homoplasy_counts_all_sites.add(site_counts);
        }

        return homoplasy_counts_all_sites;
    }


    /**
     * Homoplasy counts for a specific allele at a specific segregating site.
     * <p>
     * pheno label convention:  int[4] counts := a1p1, a1p2, a2p1, a2p2   p1 := case p2 := control (same as count_phenos() used in fisher's exact)
     *
     * @param homoplasic_mutations
     * @param phenos
     * @return homoplasy counts for this allele.  int[2]:  [0] := p1 := case, [1] := p2 := control
     */
    private int[] counts_by_allele(HashMap<String, HashSet<String>> homoplasic_mutations, HashMap<String, String> phenos) {

        int[] counts_for_allele = {0, 0};

        for (Map.Entry<String, HashSet<String>> homoplasic_mutation : homoplasic_mutations.entrySet()) {  // only non-missing genotypes are in this set of homoplasic mutations
            String node_name = homoplasic_mutation.getKey();

            // if no pheno is available for this node (i.e. this is an internal node or sample is missing pheno data), simply continue without updating any counters of the contingency table
            if (phenos.containsKey(node_name)) {
                String pheno = phenos.get(node_name);

               /* if (pheno.equals("0")) {
                    counts_for_allele[0]++;
                } else if (pheno.equals("1")) {
                    counts_for_allele[1]++;
                }*/

                if (pheno.equals("1")) {  // case
                    counts_for_allele[0]++;
                } else if (pheno.equals("0")) {  // control
                    counts_for_allele[1]++;
                } // all other pheno values are considered missing pheno data
            }
        }

        return counts_for_allele;
    }


    /**
     * Subset sites by min_hcount to use only "homoplasically informative" sites := sites where (a1 h count >= min_hcount || a2 h count >= min_hcount).  Allele homoplasic
     * counts either use only extant nodes or include homoplasic mutations from internal nodes (ie phenotypic ancestral reconstruction is in play).
     * At least one of the alleles at each segregating site must possess a min homoplasic mutation count > 0, otherwise it is statistically pointless to test this site;
     * sites where both alleles harbor 0 homoplasic mutations do not contribute to either the background p_success value nor a meaningful association test.
     * A min_hcount = 1 filters out these non-informative sites.  Keep in mind, many sites show 0 homoplasic mutations for one or both alleles, so this function merely
     * identifies informative sites.
     * <p>
     * To be clear:
     * Setting min_hcount = 1 includes all sites where there is at least one allele harboring at least 1 homoplasic mutation.
     * This is true for both extant nodes only and internal nodes modes.
     * Keep in mind that when using extant nodes only, alleles can harbor only 1 homoplasic mutation (ie the other homoplasic mutation, required to make this site
     * homoplasic, is on an internal branch).
     * <p>
     * Also, this function was useful in exploring higher thresholds like a maf filter.
     *
     * @param all_events all sites, all homoplasy events
     * @param min_hcount
     */
    private ArrayList<Homoplasy_Events> subset_by_min_hcount(ArrayList<Homoplasy_Events> all_events, int min_hcount) {
        ArrayList<Homoplasy_Events> homoplasically_informative_sites = new ArrayList<>();

        // DEBUG:  informative alleles
        int a1_informative_count = 0;
        int a2_informative_count = 0;
        int num_sites_both_alleles_informative = 0;

        for (Homoplasy_Events site : all_events) {
            int a1_count;
            int a2_count;

            if (EXTANT_NODES_ONLY) {
                a1_count = site.a1_count_extant_only;
                a2_count = site.a2_count_extant_only;
            } else {
                a1_count = site.a1_count;
                a2_count = site.a2_count;
            }

/*          old site-centric code
            if (a1_count >= min_hcount || a2_count >= min_hcount) {
                homoplasically_informative_sites.add(site);
            }
*/

            if (a1_count >= min_hcount) {
                site.a1_in_use = true;

                // DEBUG
                a1_informative_count++;
            }

            if (a2_count >= min_hcount) {
                site.a2_in_use = true;

                // DEBUG
                a2_informative_count++;
            }

            /* TESTING: block commented out to test maxT null dist behavior if sites with a tot # homoplasic mutations >= 9 is removed except for site 1044 [9,0] top hit
            if (site.segsite_ID == 1044) {  // top a2 hit [9,0] I want to keep
                if (a1_count >= min_hcount) {
                    site.a1_in_use = true;

                    // DEBUG
                    a1_informative_count++;
                }

                if (a2_count >= min_hcount) {
                    site.a2_in_use = true;

                    // DEBUG
                    a2_informative_count++;
                }
            } else {  // else filter out all other sites where tot # homoplasic extant mutations is >= 9
                if (a1_count >= min_hcount && a1_count < 9) {
                    site.a1_in_use = true;

                    // DEBUG
                    a1_informative_count++;
                }

                if (a2_count >= min_hcount && a2_count < 9) {
                    site.a2_in_use = true;

                    // DEBUG
                    a2_informative_count++;
                }
            }*/  // ********************* TESTING maxT null dist behavior *********************


            if (site.a1_in_use || site.a2_in_use) {
                homoplasically_informative_sites.add(site);
            }

            // DEBUG:  duplicates
            if (site.a1_in_use && site.a2_in_use) {
                num_sites_both_alleles_informative++;
            }
        }

        // if user has set the min_hcount too high and thus filtered out all segsites, then exit gracefully.
        if (homoplasically_informative_sites.size() == 0) {
            System.out.printf("%nmin_hcount = %s is set too high and has filtered out all segregating sites.  Try setting min_hcount to a lower value.%n%n", AlgoParams.min_hcount);
            System.exit(-1);
        }

        if (DEBUG_MODE) {
            try {
                outputOptions.debug.write("\nsubset_by_min_hcount():\n");
                outputOptions.debug.write("min_hcount = " + min_hcount + "\n");
                outputOptions.debug.write("all_events.size() = " + all_events.size() + "\n");
                outputOptions.debug.write("homoplasically_informative_sites.size() = " + homoplasically_informative_sites.size() + "\n");
                outputOptions.debug.write("a1_informative_count = " + a1_informative_count  + "\n");
                outputOptions.debug.write("a2_informative_count = " + a2_informative_count + "\n");
                outputOptions.debug.write("tot # informative alleles = " + (a1_informative_count + a2_informative_count) + "\n");
                outputOptions.debug.write("num_sites_both_alleles_informative = " + num_sites_both_alleles_informative + "\n");
                //        homoplasically_informative_sites.forEach(i -> System.out.println(i));

                // DEBUG:  duplicate segsite entries
                List<Integer> segsite_IDs = homoplasically_informative_sites.stream().map(s -> s.segsite_ID).collect(Collectors.toList());
                List<Integer> duplicate_segsites = segsite_IDs.stream().filter(ID -> Collections.frequency(segsite_IDs, ID) > 1).collect(Collectors.toList());
                outputOptions.debug.write("duplicate_segsites = " + duplicate_segsites + "\n");
                outputOptions.debug.flush();
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the debug file.  Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }

        try {
            outputOptions.log.write(String.format("%n# segregating sites with at least " + AlgoParams.min_hcount + " homoplasic mutations: " + homoplasically_informative_sites.size() + "%n"));
            outputOptions.log.flush();
        } catch (IOException e) {
            System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the log file.  Please see the stack trace in the log file for more information on this error.%n|@"));
            e.printStackTrace();
            System.exit(-1);
        }

        return homoplasically_informative_sites;
    }


    /**
     * This method produces resampling-derived point-estimates (the original phyC test stat does not exist outside of resampling) as well as resampling-derived family-wise
     * estimates.
     *
     * @param homoplasically_informative_sites
     * @param phenos
     * @return
     */
    private phyC_Test_Statistic[] phyC(ArrayList<Homoplasy_Events> homoplasically_informative_sites, HashMap<String, String> phenos) {
        phyC_Test_Statistic[] resampled_test_statistics = new phyC_Test_Statistic[homoplasically_informative_sites.size()];

        // observed counts
        ArrayList<int[]> obs_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, phenos);

        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            site.obs_counts = obs_homoplasy_counts.get(s);  // set obs counts for each Homoplasy_Events obj
            resampled_test_statistics[s] = new phyC_Test_Statistic(site.a1_in_use, site.a2_in_use);  // initialize
        }

        // DEBUG:  sanity check on EXTANT_NODES_ONLY code
//        check_extant_nodes_only_code(homoplasically_informative_sites);

        // construct both point-wise and family-wise null dists (resampling begins)
//        int m = 100;  // # resamples/replicates/perms
        // DEBUG
//        System.out.println("m = " + m);

        // construct 1 familywise null dist per allele
        float[] maxT_nulldist_a1 = new float[m];
        float[] maxT_nulldist_a2 = new float[m];

        // DEBUG:  point-wise null dist by segsite_ID
//        ArrayList<int[]> pointwise_nulldist = new ArrayList<>();  // both alleles
//        int segsite_index = debug_get_index_by_segsite(75910, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(20829, homoplasically_informative_sites);  // figure out index of this segsite
//        int segsite_index = debug_get_index_by_segsite(1142, homoplasically_informative_sites);  // figure out index of this segsite

        for (int curr_rep = 0; curr_rep < m; curr_rep++) {
            // 1 complete resampling := 3 methods:
            HashMap<String, String> permuted_phenos = permute_phenos(phenos);
            ArrayList<int[]> resampled_homoplasy_counts = count_homoplasies(homoplasically_informative_sites, permuted_phenos);
            phyC_test_stat(homoplasically_informative_sites, obs_homoplasy_counts, resampled_homoplasy_counts, resampled_test_statistics);

            // DEBUG:  point-wise null dist by segsite_ID
//            pointwise_nulldist.add(resampled_homoplasy_counts.get(segsite_index));

            // progress
//            if ((curr_rep % 100) == 0) {
//                System.out.println("replicate: " + curr_rep);
//            }

            // get the replicate maxT for each allele
            float[] replicate_maxT = maxT(homoplasically_informative_sites, resampled_homoplasy_counts);

            // allele-specific null dists (all maxT null dists should always be size m regardless of "maf" filter)
            maxT_nulldist_a1[curr_rep] = replicate_maxT[0];
            maxT_nulldist_a2[curr_rep] = replicate_maxT[1];
        }

        calc_pointwise_estimate(homoplasically_informative_sites, resampled_test_statistics, m);
        calc_familywise_estimates(homoplasically_informative_sites, resampled_test_statistics, m, maxT_nulldist_a1, maxT_nulldist_a2, obs_homoplasy_counts);

        // DIAGNOSTIC
//        diagnostic_pointwise_resamplings(segsite_index, homoplasically_informative_sites, resampled_test_statistics, pointwise_nulldist);

        // DIAGNOSTIC
//        diagnostic_resampled_case_control_ratio_by_site(homoplasically_informative_sites, segsite_index, pointwise_nulldist);

        return resampled_test_statistics;
    }


    /**
     * THIS VERSION USES THE OLD RESAMPLING SPACE (IE:  ALL MUTATIONS, INCLUDING NON-HOMOPLASIC MUTATIONS)
     * <p>
     * New association test stat:  binomial test
     * <p>
     * For all sites in this current replicate:  r++ for the allele if . . . resampled binomial test pvalue is <= observed binomial test pvalue, else r is unchanged.
     * *
     * homoplasy counts := int[4] := a1p1, a1p2, a2p1, a2p2 := {a1/case, a1/control}, {a2/case, a2/control}
     *
     * @param homoplasically_informative_sites
     * @param resampled_homoplasy_counts
     * @param resampled_test_statistics
     * @param binom_test
     * @param p_success_a1
     * @param p_success_a2
     * @return ArrayList[0] := binom test pvalues for a1, ArrayList[1] := binom test pvalues for a2
     */
    private ArrayList[] calc_binomial_test_stats(ArrayList<Homoplasy_Events> homoplasically_informative_sites, ArrayList<int[]> resampled_homoplasy_counts,
                                                 Binomial_Test_Stat[] resampled_test_statistics, BinomialTest binom_test, double p_success_a1, double p_success_a2) {

        ArrayList<Double> binom_test_stats_a1 = new ArrayList<>();  // only stores values for alleles that are in use, thus size of ArrayList will be <= # sites
        ArrayList<Double> binom_test_stats_a2 = new ArrayList<>();

        for (int s = 0; s < homoplasically_informative_sites.size(); s++) { // for each site
            Homoplasy_Events curr_site = homoplasically_informative_sites.get(s);
            Binomial_Test_Stat curr_test_stat = resampled_test_statistics[s];

            int[] curr_site_resampled_counts = resampled_homoplasy_counts.get(s);

            if (curr_site.a1_in_use) {
                // calc binom test for resampled counts
                int tot_num_trials = curr_site_resampled_counts[0] + curr_site_resampled_counts[1];  // tot # homoplasic mutations at this site for a1
                int num_cases = curr_site_resampled_counts[0];

                // REDUNDANT BINOM CALLS ============
                String k = tot_num_trials + "-" + num_cases;
                if (binom_test_calls.containsKey(k)) {  // TODO: try tertiary operator ?:
                    binom_test_calls.put(k, binom_test_calls.get(k) + 1);
                } else {
                    binom_test_calls.put(k, 1);
                }

                tot_num_calls_by_int_counter++;
                // REDUNDANT BINOM CALLS =============

                // Measure wall time for one binomialTest() call
                Instant start_time = Instant.now();
                double resampled_binom_test_pvalue = binom_test.binomialTest(tot_num_trials, num_cases, p_success_a1, TAIL_TYPE);
                tot_duration_of_all_binom_calls_to_math_library = Duration.between(start_time, Instant.now()).plus(tot_duration_of_all_binom_calls_to_math_library);
                binom_test_stats_a1.add(resampled_binom_test_pvalue);

                // if resampled test stat equal or more extreme than obs test stat
                if (resampled_binom_test_pvalue <= curr_test_stat.obs_binom_pvalue_a1) {
                    curr_test_stat.r_a1++;
                }
            }

            if (curr_site.a2_in_use) {
                // calc binom test for resampled counts
                int tot_num_trials = curr_site_resampled_counts[2] + curr_site_resampled_counts[3];  // tot # homoplasic mutations at this site for a1
                int num_cases = curr_site_resampled_counts[2];

                // REDUNDANT BINOM CALLS ==========
                String k = tot_num_trials + "-" + num_cases;
                if (binom_test_calls.containsKey(k)) {
                    binom_test_calls.put(k, binom_test_calls.get(k) + 1);
                } else {
                    binom_test_calls.put(k, 1);
                }

                tot_num_calls_by_int_counter++;
                // REDUNDANT BINOM CALLS ============

                // Measure wall time for one binomialTest() call
                Instant start_time = Instant.now();
                double resampled_binom_test_pvalue = binom_test.binomialTest(tot_num_trials, num_cases, p_success_a2, TAIL_TYPE);
                tot_duration_of_all_binom_calls_to_math_library = Duration.between(start_time, Instant.now()).plus(tot_duration_of_all_binom_calls_to_math_library);
                binom_test_stats_a2.add(resampled_binom_test_pvalue);

                // if resampled test stat equal or more extreme than obs test stat
                if (resampled_binom_test_pvalue <= curr_test_stat.obs_binom_pvalue_a2) {
                    curr_test_stat.r_a2++;
                }
            }
        }

        return new ArrayList[]{binom_test_stats_a1, binom_test_stats_a2};
    }


    /**
     * THIS VERSION USES THE OLD RESAMPLING SPACE (IE:  ALL MUTATIONS, INCLUDING NON-HOMOPLASIC MUTATIONS).
     * THIS VERSION NOW ALSO IS OPTIMIZED USING MEMOIZATION of resampled pavlues.
     * <p>
     * New association test stat:  binomial test
     * <p>
     * For all sites in this current replicate:  r++ for the allele if . . . resampled binomial test pvalue is <= observed binomial test pvalue, else r is unchanged.
     * *
     * homoplasy counts := int[4] := a1p1, a1p2, a2p1, a2p2 := {a1/case, a1/control}, {a2/case, a2/control}
     *
     * @param homoplasically_informative_sites
     * @param resampled_homoplasy_counts
     * @param resampled_test_statistics
     * @param binom_test
     * @param p_success_a1
     * @param p_success_a2
     * @param resampled_pvalues_cache
     * @return ArrayList[0] := binom test pvalues for a1, ArrayList[1] := binom test pvalues for a2
     */
    private ArrayList[] calc_binomial_test_stats_memoization(ArrayList<Homoplasy_Events> homoplasically_informative_sites, ArrayList<int[]> resampled_homoplasy_counts,
                                                             Binomial_Test_Stat[] resampled_test_statistics, BinomialTest binom_test, double p_success_a1, double p_success_a2,
                                                             HashMap<String, Double> resampled_pvalues_cache) {

        ArrayList<Double> binom_test_stats_a1 = new ArrayList<>();  // only stores values for alleles that are in use, thus size of ArrayList will be <= # sites
        ArrayList<Double> binom_test_stats_a2 = new ArrayList<>();

        for (int s = 0; s < homoplasically_informative_sites.size(); s++) { // for each site
            Homoplasy_Events curr_site = homoplasically_informative_sites.get(s);
            Binomial_Test_Stat curr_test_stat = resampled_test_statistics[s];

            int[] curr_site_resampled_counts = resampled_homoplasy_counts.get(s);

            if (curr_site.a1_in_use) {
                // calc binom test for resampled counts
                int tot_num_trials = curr_site_resampled_counts[0] + curr_site_resampled_counts[1];  // tot # homoplasic mutations at this site for a1
                int num_cases = curr_site_resampled_counts[0];

                /*// REDUNDANT CALLS ============
                String k = tot_num_trials + "-" + num_cases;
                if (binom_test_calls.containsKey(k)) {
                    binom_test_calls.put(k, binom_test_calls.get(k) + 1);
                } else {
                    binom_test_calls.put(k, 1);
                }
                tot_num_calls_by_int_counter++;
                // REDUNDANT CALLS =============*/

                // String experiment
                String k = new StringBuilder().append(tot_num_trials).append("_").append(num_cases).toString();
//                String k = tot_num_trials + "-" + num_cases;  // the 2 params that make a binom call unique
                double resampled_binom_test_pvalue = 0;
                if (resampled_pvalues_cache.containsKey(k)) {  // retrieve pvalue from cache or make call to Apache match library
                    resampled_binom_test_pvalue = resampled_pvalues_cache.get(k);
                } else {
//                    Instant start_time = Instant.now();
                    resampled_binom_test_pvalue = binom_test.binomialTest(tot_num_trials, num_cases, p_success_a1, TAIL_TYPE);
//                    tot_duration_of_all_binom_calls_to_math_library = Duration.between(start_time, Instant.now()).plus(tot_duration_of_all_binom_calls_to_math_library);
                    resampled_pvalues_cache.put(k, resampled_binom_test_pvalue);
                }
                binom_test_stats_a1.add(resampled_binom_test_pvalue);

                // if resampled test stat equal or more extreme than obs test stat
                if (resampled_binom_test_pvalue <= curr_test_stat.obs_binom_pvalue_a1) {
                    curr_test_stat.r_a1++;
                }
            }

            if (curr_site.a2_in_use) {
                // calc binom test for resampled counts
                int tot_num_trials = curr_site_resampled_counts[2] + curr_site_resampled_counts[3];  // tot # homoplasic mutations at this site for a1
                int num_cases = curr_site_resampled_counts[2];

                /*// REDUNDANT CALLS ==========
                String k = tot_num_trials + "-" + num_cases;
                if (binom_test_calls.containsKey(k)) {
                    binom_test_calls.put(k, binom_test_calls.get(k) + 1);
                } else {
                    binom_test_calls.put(k, 1);
                }
                tot_num_calls_by_int_counter++;
                // REDUNDANT CALLS ============*/

                String k = new StringBuilder().append(tot_num_trials).append("_").append(num_cases).toString();
//                String k = tot_num_trials + "-" + num_cases;
                double resampled_binom_test_pvalue = 0;
                if (resampled_pvalues_cache.containsKey(k)) {
                    resampled_binom_test_pvalue = resampled_pvalues_cache.get(k);
                } else {
//                    Instant start_time = Instant.now();
                    resampled_binom_test_pvalue = binom_test.binomialTest(tot_num_trials, num_cases, p_success_a2, TAIL_TYPE);
//                    tot_duration_of_all_binom_calls_to_math_library = Duration.between(start_time, Instant.now()).plus(tot_duration_of_all_binom_calls_to_math_library);
                    resampled_pvalues_cache.put(k, resampled_binom_test_pvalue);
                }
                binom_test_stats_a2.add(resampled_binom_test_pvalue);

                // if resampled test stat equal or more extreme than obs test stat
                if (resampled_binom_test_pvalue <= curr_test_stat.obs_binom_pvalue_a2) {
                    curr_test_stat.r_a2++;
                }
            }
        }

        return new ArrayList[]{binom_test_stats_a1, binom_test_stats_a2};
    }


    /**
     * THIS VERSION USES THE OLD RESAMPLING SPACE (IE:  ALL MUTATIONS, INCLUDING NON-HOMOPLASIC MUTATIONS).
     * THIS VERSION NOW ALSO IS OPTIMIZED USING MEMOIZATION of resampled pavlues.
     * <p>
     * New association test stat:  binomial test
     * <p>
     * For all sites in this current replicate:  r++ for the allele if . . . resampled binomial test pvalue is <= observed binomial test pvalue, else r is unchanged.
     * *
     * homoplasy counts := int[4] := a1p1, a1p2, a2p1, a2p2 := {a1/case, a1/control}, {a2/case, a2/control}
     *
     * @param homoplasically_informative_sites
     * @param resampled_homoplasy_counts
     * @param resampled_test_statistics
     * @param binom_test
     * @param p_success_a1
     * @param p_success_a2
     * @param resampled_pvalues_cache
     * @return ArrayList[0] := binom test pvalues for a1, ArrayList[1] := binom test pvalues for a2
     */
    private ArrayList[] calc_binomial_test_stats_memoization_concurrent(ArrayList<Homoplasy_Events> homoplasically_informative_sites, ArrayList<int[]> resampled_homoplasy_counts,
                                                                        Binomial_Test_Stat[] resampled_test_statistics, BinomialTest binom_test, double p_success_a1, double p_success_a2,
                                                                        ConcurrentHashMap<String, Double> resampled_pvalues_cache) {

        ArrayList<Double> binom_test_stats_a1 = new ArrayList<>();  // only stores values for alleles that are in use, thus size of ArrayList will be <= # sites
        ArrayList<Double> binom_test_stats_a2 = new ArrayList<>();

        for (int s = 0; s < homoplasically_informative_sites.size(); s++) { // for each site
            Homoplasy_Events curr_site = homoplasically_informative_sites.get(s);
            Binomial_Test_Stat curr_test_stat = resampled_test_statistics[s];

            int[] curr_site_resampled_counts = resampled_homoplasy_counts.get(s);

            if (curr_site.a1_in_use) {
                // calc binom test for resampled counts
                int tot_num_trials = curr_site_resampled_counts[0] + curr_site_resampled_counts[1];  // tot # homoplasic mutations at this site for a1
                int num_cases = curr_site_resampled_counts[0];

                /*// REDUNDANT CALLS ============
                String k = tot_num_trials + "-" + num_cases;
                if (binom_test_calls.containsKey(k)) {
                    binom_test_calls.put(k, binom_test_calls.get(k) + 1);
                } else {
                    binom_test_calls.put(k, 1);
                }
                tot_num_calls_by_int_counter++;
                // REDUNDANT CALLS =============*/

                // String experiment
                String k = new StringBuilder().append(tot_num_trials).append("_").append(num_cases).toString();
//                String k = tot_num_trials + "-" + num_cases;  // the 2 params that make a binom call unique
                double resampled_binom_test_pvalue = 0;
                if (resampled_pvalues_cache.containsKey(k)) {  // retrieve pvalue from cache or make call to Apache match library
                    resampled_binom_test_pvalue = resampled_pvalues_cache.get(k);
                } else {
//                    Instant start_time = Instant.now();
                    resampled_binom_test_pvalue = binom_test.binomialTest(tot_num_trials, num_cases, p_success_a1, TAIL_TYPE);
//                    tot_duration_of_all_binom_calls_to_math_library = Duration.between(start_time, Instant.now()).plus(tot_duration_of_all_binom_calls_to_math_library);
                    resampled_pvalues_cache.put(k, resampled_binom_test_pvalue);
                }
                binom_test_stats_a1.add(resampled_binom_test_pvalue);

                // if resampled test stat equal or more extreme than obs test stat
                if (resampled_binom_test_pvalue <= curr_test_stat.obs_binom_pvalue_a1) {
                    curr_test_stat.r_a1++;
                }
            }

            if (curr_site.a2_in_use) {
                // calc binom test for resampled counts
                int tot_num_trials = curr_site_resampled_counts[2] + curr_site_resampled_counts[3];  // tot # homoplasic mutations at this site for a1
                int num_cases = curr_site_resampled_counts[2];

                /*// REDUNDANT CALLS ==========
                String k = tot_num_trials + "-" + num_cases;
                if (binom_test_calls.containsKey(k)) {
                    binom_test_calls.put(k, binom_test_calls.get(k) + 1);
                } else {
                    binom_test_calls.put(k, 1);
                }
                tot_num_calls_by_int_counter++;
                // REDUNDANT CALLS ============*/

                String k = new StringBuilder().append(tot_num_trials).append("_").append(num_cases).toString();
//                String k = tot_num_trials + "-" + num_cases;
                double resampled_binom_test_pvalue = 0;
                if (resampled_pvalues_cache.containsKey(k)) {
                    resampled_binom_test_pvalue = resampled_pvalues_cache.get(k);
                } else {
//                    Instant start_time = Instant.now();
                    resampled_binom_test_pvalue = binom_test.binomialTest(tot_num_trials, num_cases, p_success_a2, TAIL_TYPE);
//                    tot_duration_of_all_binom_calls_to_math_library = Duration.between(start_time, Instant.now()).plus(tot_duration_of_all_binom_calls_to_math_library);
                    resampled_pvalues_cache.put(k, resampled_binom_test_pvalue);
                }
                binom_test_stats_a2.add(resampled_binom_test_pvalue);

                // if resampled test stat equal or more extreme than obs test stat
                if (resampled_binom_test_pvalue <= curr_test_stat.obs_binom_pvalue_a2) {
                    curr_test_stat.r_a2++;
                }
            }
        }

        return new ArrayList[]{binom_test_stats_a1, binom_test_stats_a2};
    }


    /**
     * Point-wise null dist by segsite_ID.
     *
     * @param segsite_index
     * @param homoplasically_informative_sites
     * @param resampled_test_statistics
     * @param pointwise_nulldist
     */
    private void diagnostic_pointwise_resamplings(int segsite_index, ArrayList<Homoplasy_Events> homoplasically_informative_sites, phyC_Test_Statistic[] resampled_test_statistics,
                                                  ArrayList<int[]> pointwise_nulldist) {
        System.out.printf("%n%nDIAGNOSTIC:  point-wise resamplings:%n");
        System.out.println(Homoplasy_Events.COL_NAMES + "\t" + phyC_Test_Statistic.COL_NAMES);
        System.out.print(homoplasically_informative_sites.get(segsite_index).toString() + "\t" + resampled_test_statistics[segsite_index].toString());
        System.out.println("point-wise null dist at segsite_ID " + homoplasically_informative_sites.get(segsite_index).segsite_ID);
        pointwise_nulldist.forEach(site_counts -> System.out.println(Arrays.toString(site_counts) + "   " + (site_counts[0] + site_counts[1])));
    }


    /**
     * Resampled case:control ratio for a specific segsite.
     *
     * @param homoplasically_informative_sites
     * @param segsite_index
     * @param pointwise_nulldist
     */
    private void diagnostic_resampled_case_control_ratio_by_site(ArrayList<Homoplasy_Events> homoplasically_informative_sites, int segsite_index, ArrayList<int[]> pointwise_nulldist) {
        System.out.printf("%n%nDIAGNOSTIC:  case:control ratio by site%n");
        System.out.println("point-wise null dist at segsite_ID " + homoplasically_informative_sites.get(segsite_index).segsite_ID);
        System.out.println(homoplasically_informative_sites.get(segsite_index));

        int tot_cases_a1 = 0;
        int tot_controls_a1 = 0;

        int tot_cases_a2 = 0;
        int tot_controls_a2 = 0;
        for (int[] resampled_counts : pointwise_nulldist) {
            tot_cases_a1 += resampled_counts[0];
            tot_controls_a1 += resampled_counts[1];

            tot_cases_a2 += resampled_counts[2];
            tot_controls_a2 += resampled_counts[3];
        }

        System.out.println("tot_cases_a1 = " + tot_cases_a1);
        System.out.println("tot_controls_a1 = " + tot_controls_a1);
        System.out.println("tot_cases_a2 = " + tot_cases_a2);
        System.out.println("tot_controls_a2 = " + tot_controls_a2);
        float resampled_case_control_ratio_a1 = (float) tot_cases_a1 / (float) tot_controls_a1;
        float resampled_case_control_ratio_a2 = (float) tot_cases_a2 / (float) tot_controls_a2;
        System.out.println("resampled_case_control_ratio_a1 = " + resampled_case_control_ratio_a1);
        System.out.println("resampled_case_control_ratio_a2 = " + resampled_case_control_ratio_a2);
        System.out.println();
    }


    /**
     * This should really be a unit test, but for the sake of time right now I'll writing it here.
     * If EXTANT_NODES_ONLY code I just wrote works, then:
     * sum of all obs_counts must equal a1 + a2 extant only counts
     *
     * @param all_events
     */
    private void check_extant_nodes_only_code(ArrayList<Homoplasy_Events> all_events) {
        boolean sanity_check = true;
        Homoplasy_Events failed_site = null;

        for (Homoplasy_Events site : all_events) {
            int[] obs_counts = site.obs_counts;
            int tot_obs_counts = obs_counts[0] + obs_counts[1] + obs_counts[2] + obs_counts[3];

            if (tot_obs_counts != (site.a1_count_extant_only + site.a2_count_extant_only)) {
                sanity_check = false;
                failed_site = site;
                break;
            }
        }

        if (!sanity_check) {
            System.out.println("EXTANT_NODES_ONLY code failed sanity check");
            System.out.println(Homoplasy_Events.COL_NAMES);
            System.out.println("failed_site = " + failed_site);
            System.exit(-1);
        } else {
            System.out.println("EXTANT_NODES_ONLY code passed sanity check");
            System.exit(1);
        }
    }


    private int debug_get_index_by_segsite(int segsite_ID, ArrayList<Homoplasy_Events> homoplasically_informative_sites) {
        int index = -1;

        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            if (homoplasically_informative_sites.get(i).segsite_ID == segsite_ID) {
                index = i;
            }
        }

        return index;
    }


    /**
     * For all sites, for each allele:  compares obs vs maxT null dist.  Updates r_familywise (for all alleles) for each site.
     *
     * @param homoplasically_informative_sites
     * @param resampled_test_statistics
     * @param m
     * @param maxT_nulldist_a1
     * @param maxT_nulldist_a2
     * @param obs_homoplasy_counts
     */
    private void calc_familywise_estimates(ArrayList<Homoplasy_Events> homoplasically_informative_sites, phyC_Test_Statistic[] resampled_test_statistics, int m,
                                           float[] maxT_nulldist_a1, float[] maxT_nulldist_a2, ArrayList<int[]> obs_homoplasy_counts) {

        // this is only needed for diagnostics
        // sort all maxT null dists
        Arrays.parallelSort(maxT_nulldist_a1);
        Arrays.parallelSort(maxT_nulldist_a2);

        // DEBUG
        System.out.println("maxT_nulldist_a1.length = " + maxT_nulldist_a1.length);
        System.out.println("maxT_nulldist_a2.length = " + maxT_nulldist_a2.length);
        System.out.println("dn_maxTs for maxT nulldist a1: " + Arrays.toString(maxT_nulldist_a1));
        System.out.println("dn_maxTs for maxT nulldist a2: " + Arrays.toString(maxT_nulldist_a2));

/*      old site-centric code:
        // for each site's alleles, calc r familywise
        for (int curr_site = 0; curr_site < resampled_test_statistics.length; curr_site++) {
            int[] obs_site_counts = obs_homoplasy_counts.get(curr_site);

            float[] obs_dn = dn(obs_site_counts);
            float obs_dn_a1 = obs_dn[0];
            float obs_dn_a2 = obs_dn[1];

            // DEBUG:  obs values
//            System.out.println("obs_site_counts = " + Arrays.toString(obs_site_counts));
//            System.out.println("obs_dn[] = " + Arrays.toString(obs_dn));

            DoubleStream ds_a1 = IntStream.range(0, maxT_nulldist_a1.length)
                    .mapToDouble(i -> maxT_nulldist_a1[i]);
            DoubleStream ds_a2 = IntStream.range(0, maxT_nulldist_a2.length)
                    .mapToDouble(i -> maxT_nulldist_a2[i]);
            long r_maxT_a1 = ds_a1.filter(dn_maxT -> dn_maxT >= obs_dn_a1).count();
            long r_maxT_a2 = ds_a2.filter(dn_maxT -> dn_maxT >= obs_dn_a2).count();

            double familywise_pvalue_a1 = ((double)r_maxT_a1 + 1) / ((double)m + 1);
            double familywise_pvalue_a2 = ((double)r_maxT_a2 + 1) / ((double)m + 1);

            resampled_test_statistics[curr_site].obs_dn_a1 = obs_dn_a1;
            resampled_test_statistics[curr_site].obs_dn_a2 = obs_dn_a2;
            resampled_test_statistics[curr_site].r_maxT_a1 = Math.toIntExact(r_maxT_a1);  // TODO:  catch ArithmeticException to output a cleaner explanation of overflow?
            resampled_test_statistics[curr_site].r_maxT_a2 = Math.toIntExact(r_maxT_a2);
            resampled_test_statistics[curr_site].familywise_pvalue_a1 = familywise_pvalue_a1;
            resampled_test_statistics[curr_site].familywise_pvalue_a2 = familywise_pvalue_a2;
        }
*/

        // for each site, calc resampling-based FWER
        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(i);
            int[] obs_site_counts = obs_homoplasy_counts.get(i);

            float[] obs_dn = dn(obs_site_counts);  // TODO ALLELE-CENTRIC:  should dn() use both alleles always?

            if (site.a1_in_use) {
                float obs_dn_a1 = obs_dn[0];

                DoubleStream ds_a1 = IntStream.range(0, maxT_nulldist_a1.length)
                        .mapToDouble(j -> maxT_nulldist_a1[j]);
                long r_maxT_a1 = ds_a1.filter(dn_maxT -> dn_maxT >= obs_dn_a1).count();

                double familywise_pvalue_a1 = ((double) r_maxT_a1 + 1) / ((double) m + 1);

                resampled_test_statistics[i].obs_dn_a1 = obs_dn_a1;
                resampled_test_statistics[i].r_maxT_a1 = Math.toIntExact(r_maxT_a1);
                resampled_test_statistics[i].familywise_pvalue_a1 = familywise_pvalue_a1;
            }

            if (site.a2_in_use) {
                float obs_dn_a2 = obs_dn[1];

                DoubleStream ds_a2 = IntStream.range(0, maxT_nulldist_a2.length)
                        .mapToDouble(j -> maxT_nulldist_a2[j]);
                long r_maxT_a2 = ds_a2.filter(dn_maxT -> dn_maxT >= obs_dn_a2).count();

                double familywise_pvalue_a2 = ((double) r_maxT_a2 + 1) / ((double) m + 1);

                resampled_test_statistics[i].obs_dn_a2 = obs_dn_a2;
                resampled_test_statistics[i].r_maxT_a2 = Math.toIntExact(r_maxT_a2);
                resampled_test_statistics[i].familywise_pvalue_a2 = familywise_pvalue_a2;
            }
        }
    }


    private void calc_familywise_estimates_binom(ArrayList<Homoplasy_Events> homoplasically_informative_sites, Binomial_Test_Stat[] resampled_test_statistics, int m,
                                                 double[] maxT_nulldist_a1, double[] maxT_nulldist_a2) {

        // this is only needed for diagnostics
        // sort all maxT null dists
        Arrays.parallelSort(maxT_nulldist_a1);
        Arrays.parallelSort(maxT_nulldist_a2);

        // DEBUG
        System.out.println("maxT_nulldist_a1.length = " + maxT_nulldist_a1.length);
        System.out.println("maxT_nulldist_a2.length = " + maxT_nulldist_a2.length);
        System.out.println("binom pvalue maxTs for maxT nulldist a1: " + Arrays.toString(maxT_nulldist_a1));
        System.out.println("binom pvalue maxTs for maxT nulldist a2: " + Arrays.toString(maxT_nulldist_a2));

/*      old site-centric code:
        // for each site's alleles, calc r familywise
        for (int curr_site = 0; curr_site < resampled_test_statistics.length; curr_site++) {
            int[] obs_site_counts = obs_homoplasy_counts.get(curr_site);

            float[] obs_dn = dn(obs_site_counts);
            float obs_dn_a1 = obs_dn[0];
            float obs_dn_a2 = obs_dn[1];

            // DEBUG:  obs values
//            System.out.println("obs_site_counts = " + Arrays.toString(obs_site_counts));
//            System.out.println("obs_dn[] = " + Arrays.toString(obs_dn));

            DoubleStream ds_a1 = IntStream.range(0, maxT_nulldist_a1.length)
                    .mapToDouble(i -> maxT_nulldist_a1[i]);
            DoubleStream ds_a2 = IntStream.range(0, maxT_nulldist_a2.length)
                    .mapToDouble(i -> maxT_nulldist_a2[i]);
            long r_maxT_a1 = ds_a1.filter(dn_maxT -> dn_maxT >= obs_dn_a1).count();
            long r_maxT_a2 = ds_a2.filter(dn_maxT -> dn_maxT >= obs_dn_a2).count();

            double familywise_pvalue_a1 = ((double)r_maxT_a1 + 1) / ((double)m + 1);
            double familywise_pvalue_a2 = ((double)r_maxT_a2 + 1) / ((double)m + 1);

            resampled_test_statistics[curr_site].obs_dn_a1 = obs_dn_a1;
            resampled_test_statistics[curr_site].obs_dn_a2 = obs_dn_a2;
            resampled_test_statistics[curr_site].r_maxT_a1 = Math.toIntExact(r_maxT_a1);  // TODO:  catch ArithmeticException to output a cleaner explanation of overflow?
            resampled_test_statistics[curr_site].r_maxT_a2 = Math.toIntExact(r_maxT_a2);
            resampled_test_statistics[curr_site].familywise_pvalue_a1 = familywise_pvalue_a1;
            resampled_test_statistics[curr_site].familywise_pvalue_a2 = familywise_pvalue_a2;
        }
*/

        // for each site, calc resampling-based FWER
        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            Binomial_Test_Stat test_stat = resampled_test_statistics[s];


            if (site.a1_in_use) {
                DoubleStream ds_a1 = Arrays.stream(maxT_nulldist_a1);

                // TODO BINOM:  test idea that point estimates capture an important distributional shift that the obs does not
                long r_maxT_a1 = ds_a1.filter(maxT -> maxT <= test_stat.obs_binom_pvalue_a1).count();
//                long r_maxT_a1 = ds_a1.filter(maxT -> maxT <= test_stat.pointwise_pvalue_a1).count();

                double familywise_pvalue_a1 = ((double) r_maxT_a1 + 1) / ((double) m + 1);
                resampled_test_statistics[s].r_maxT_a1 = Math.toIntExact(r_maxT_a1);
                resampled_test_statistics[s].familywise_pvalue_a1 = familywise_pvalue_a1;
            }

            if (site.a2_in_use) {
                DoubleStream ds_a2 = Arrays.stream(maxT_nulldist_a2);

                // TODO BINOM:  test idea
                long r_maxT_a2 = ds_a2.filter(maxT -> maxT <= test_stat.obs_binom_pvalue_a2).count();
//                long r_maxT_a2 = ds_a2.filter(maxT -> maxT <= test_stat.pointwise_pvalue_a2).count();

                double familywise_pvalue_a2 = ((double) r_maxT_a2 + 1) / ((double) m + 1);
                resampled_test_statistics[s].r_maxT_a2 = Math.toIntExact(r_maxT_a2);
                resampled_test_statistics[s].familywise_pvalue_a2 = familywise_pvalue_a2;
            }
        }
    }


    private void calc_familywise_estimates_binom_combined_nulldists(ArrayList<Homoplasy_Events> homoplasically_informative_sites,
                                                                    Binomial_Test_Stat[] resampled_test_statistics,
                                                                    int m,
                                                                    double[] maxT_nulldist_a1_a2_combined) {

        System.out.printf(Ansi.AUTO.string("@|fg(123) Calculating family-wise estimates . . . |@"));

        // this is only needed for diagnostics
        // sort all maxT null dists
//        Arrays.parallelSort(maxT_nulldist_a1);
//        Arrays.parallelSort(maxT_nulldist_a2);
        Arrays.parallelSort(maxT_nulldist_a1_a2_combined);

        // DEBUG
//        System.out.println("maxT_nulldist_a1.length = " + maxT_nulldist_a1.length);
//        System.out.println("maxT_nulldist_a2.length = " + maxT_nulldist_a2.length);
//        System.out.println("binom pvalue maxTs for maxT nulldist a1: " + Arrays.toString(maxT_nulldist_a1));
//        System.out.println("binom pvalue maxTs for maxT nulldist a2: " + Arrays.toString(maxT_nulldist_a2));
        if (DEBUG_MODE) {
            try {
                outputOptions.debug.write("\ncalc_familywise_estimates_binom_combined_nulldists():\n");
                outputOptions.debug.write("maxT_nulldist_a1_a2_combined.length = " + maxT_nulldist_a1_a2_combined.length + "\n");
                outputOptions.debug.flush();
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the debug file.  Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }
//        System.out.println("binom pvalue maxTs for maxT nulldist a1 and a2 combined: " + Arrays.toString(maxT_nulldist_a1_a2_combined));


/*      old site-centric code:
        // for each site's alleles, calc r familywise
        for (int curr_site = 0; curr_site < resampled_test_statistics.length; curr_site++) {
            int[] obs_site_counts = obs_homoplasy_counts.get(curr_site);

            float[] obs_dn = dn(obs_site_counts);
            float obs_dn_a1 = obs_dn[0];
            float obs_dn_a2 = obs_dn[1];

            // DEBUG:  obs values
//            System.out.println("obs_site_counts = " + Arrays.toString(obs_site_counts));
//            System.out.println("obs_dn[] = " + Arrays.toString(obs_dn));

            DoubleStream ds_a1 = IntStream.range(0, maxT_nulldist_a1.length)
                    .mapToDouble(i -> maxT_nulldist_a1[i]);
            DoubleStream ds_a2 = IntStream.range(0, maxT_nulldist_a2.length)
                    .mapToDouble(i -> maxT_nulldist_a2[i]);
            long r_maxT_a1 = ds_a1.filter(dn_maxT -> dn_maxT >= obs_dn_a1).count();
            long r_maxT_a2 = ds_a2.filter(dn_maxT -> dn_maxT >= obs_dn_a2).count();

            double familywise_pvalue_a1 = ((double)r_maxT_a1 + 1) / ((double)m + 1);
            double familywise_pvalue_a2 = ((double)r_maxT_a2 + 1) / ((double)m + 1);

            resampled_test_statistics[curr_site].obs_dn_a1 = obs_dn_a1;
            resampled_test_statistics[curr_site].obs_dn_a2 = obs_dn_a2;
            resampled_test_statistics[curr_site].r_maxT_a1 = Math.toIntExact(r_maxT_a1);  // TODO:  catch ArithmeticException to output a cleaner explanation of overflow?
            resampled_test_statistics[curr_site].r_maxT_a2 = Math.toIntExact(r_maxT_a2);
            resampled_test_statistics[curr_site].familywise_pvalue_a1 = familywise_pvalue_a1;
            resampled_test_statistics[curr_site].familywise_pvalue_a2 = familywise_pvalue_a2;
        }
*/

        // for each site, calc resampling-based FWER
        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            Binomial_Test_Stat test_stat = resampled_test_statistics[s];


            if (site.a1_in_use) {
                DoubleStream ds_maxT_nulldist = Arrays.stream(maxT_nulldist_a1_a2_combined);

                // TODO BINOM:  test idea that point estimates capture an important distributional shift that the obs does not
                long r_maxT_a1 = ds_maxT_nulldist.filter(maxT -> maxT <= test_stat.obs_binom_pvalue_a1).count();
//                long r_maxT_a1 = ds_a1.filter(maxT -> maxT <= test_stat.pointwise_pvalue_a1).count();

                double familywise_pvalue_a1 = ((double) r_maxT_a1 + 1) / ((double) m + 1);
                resampled_test_statistics[s].r_maxT_a1 = Math.toIntExact(r_maxT_a1);
                resampled_test_statistics[s].familywise_pvalue_a1 = familywise_pvalue_a1;
            }

            if (site.a2_in_use) {
                DoubleStream ds_maxT_nulldist = Arrays.stream(maxT_nulldist_a1_a2_combined);

                // TODO BINOM:  test idea
                long r_maxT_a2 = ds_maxT_nulldist.filter(maxT -> maxT <= test_stat.obs_binom_pvalue_a2).count();
//                long r_maxT_a2 = ds_a2.filter(maxT -> maxT <= test_stat.pointwise_pvalue_a2).count();

                double familywise_pvalue_a2 = ((double) r_maxT_a2 + 1) / ((double) m + 1);
                resampled_test_statistics[s].r_maxT_a2 = Math.toIntExact(r_maxT_a2);
                resampled_test_statistics[s].familywise_pvalue_a2 = familywise_pvalue_a2;
            }
        }

        System.out.printf("family-wise estimates complete.%n");
    }


    /**
     * Calculates normalized d for this site, all alleles.
     *
     * @param site_counts
     * @return [a1, a2]
     */
    private float[] dn(int[] site_counts) {

        // int[4] counts := a1p1, a1p2, a2p1, a2p2
        float x1 = site_counts[0];
        float y1 = site_counts[1];
        float x2 = site_counts[2];
        float y2 = site_counts[3];

        // ALLELE-CENTRIC:  what is the default behavior of NaN values below?

        float dn_a1 = (x1 - y1) / (x1 + y1);
        float dn_a2 = (x2 - y2) / (x2 + y2);

        return new float[]{dn_a1, dn_a2};
    }


    /**
     * Returns the maxT for this replicate, one for each allele.
     * Ordering of test stats is based upon d normalized := (x - y) / (x + y)
     *
     * @param homoplasically_informative_sites
     * @param resampled_homoplasy_counts
     * @return
     */
    private float[] maxT(ArrayList<Homoplasy_Events> homoplasically_informative_sites, ArrayList<int[]> resampled_homoplasy_counts) {
        float max_a1 = -2.0f;  // as long as this initialized value is < min(range) of dn (range is -1 to 1)
        float max_a2 = -2.0f;

/*
        // DEBUG:  how many informative alleles are there?
        int informative_a1_alleles = 0;
        int informative_a2_alleles = 0;
        int i = 0;
*/

/*      old site-centric code
        for (int[] site_counts : resampled_homoplasy_counts) {
            float[] curr_maxT = dn(site_counts);  // calc d normalized for both alleles
            float curr_dn_a1 = curr_maxT[0];
            float curr_dn_a2 = curr_maxT[1];

            // DEBUG:  only use informative alleles (remember we're working with only extant nodes at this point).  Won't need this once I move to allele-centric code.
            if ((site_counts[0] + site_counts[1]) >= 6) {
                informative_a1_alleles++;

                if (curr_dn_a1 > max_a1) {
                    max_a1 = curr_dn_a1;
                }
            }

            // DEBUG:  only use informative alleles (remember we're working with only extant nodes at this point).  Won't need this once I move to allele-centric code.
            if ((site_counts[2] + site_counts[3]) >= 6) {
                informative_a2_alleles++;

                if (curr_dn_a2 > max_a2) {
                    max_a2 = curr_dn_a2;
                }
            }

            // DEBUG
//            if (max_a1 == 1) {
//                System.out.print("site = " + i);
//                System.out.print("     max_a1 = " + max_a1);
//                System.out.print("     resampled site_counts = " + Arrays.toString(site_counts));
//                System.out.print("     curr_dn_a1 = " + curr_dn_a1);
//                System.out.println("     curr_dn_a2 = " + curr_dn_a2);
//            }

//             DEBUG
//            i++;
        }
*/

        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(i);
            int[] resampled_site_counts = resampled_homoplasy_counts.get(i);

            float[] curr_maxT = dn(resampled_site_counts);  // calc d normalized for both alleles
            float curr_dn_a1 = curr_maxT[0];
            float curr_dn_a2 = curr_maxT[1];

            if (site.a1_in_use) {
                // DEBUG
//                informative_a1_alleles++;

                if (curr_dn_a1 > max_a1) {
                    max_a1 = curr_dn_a1;
                }
            }

            if (site.a2_in_use) {
                // DEBUG
//                informative_a2_alleles++;

                if (curr_dn_a2 > max_a2) {
                    max_a2 = curr_dn_a2;
                }
            }
        }

/*
        // DEBUG kill so I only see 1 replicate
        System.out.println("informative_a1_alleles = " + informative_a1_alleles);
        System.out.println("informative_a2_alleles = " + informative_a2_alleles);
        int tot_informative_alleles = informative_a1_alleles + informative_a2_alleles;
        System.out.println("tot # informative alleles = " + tot_informative_alleles);
        System.out.println("max_a1 = " + max_a1);
        System.out.println("max_a2 = " + max_a2);
        System.exit(-1);
*/

        return new float[]{max_a1, max_a2};
    }


    /**
     * Calculate resampling-derived point-estimates across all sites and all replicates.
     * <p>
     * point-wise estimate = (r + 1) / (m + 1)
     *
     * @param homoplasically_informative_sites
     * @param resampled_test_statistics
     * @param m
     */
    private void calc_pointwise_estimate(ArrayList<Homoplasy_Events> homoplasically_informative_sites, phyC_Test_Statistic[] resampled_test_statistics, int m) {

/*      old site-centric code
        for (phyC_Test_Statistic site : resampled_test_statistics) {
            site.pointwise_pvalue_a1 = ((double)(site.r_a1 + 1)) / ((double)(m + 1));
            site.pointwise_pvalue_a2 = ((double)(site.r_a2 + 1)) / ((double)(m + 1));
        }
*/

        for (int i = 0; i < homoplasically_informative_sites.size(); i++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(i);
            phyC_Test_Statistic resampled_test_stat = resampled_test_statistics[i];

            if (site.a1_in_use) {
                resampled_test_stat.pointwise_pvalue_a1 = ((double) (resampled_test_stat.r_a1 + 1)) / ((double) (m + 1));
            }

            if (site.a2_in_use) {
                resampled_test_stat.pointwise_pvalue_a2 = ((double) (resampled_test_stat.r_a2 + 1)) / ((double) (m + 1));
            }
        }
    }


    private void calc_pointwise_estimate_binom(ArrayList<Homoplasy_Events> homoplasically_informative_sites, Binomial_Test_Stat[] resampled_test_statistics, int m) {

        System.out.printf(Ansi.AUTO.string("@|fg(213) %nCalculating point-wise estimates . . . |@"));

/*      old site-centric code
        for (phyC_Test_Statistic site : resampled_test_statistics) {
            site.pointwise_pvalue_a1 = ((double)(site.r_a1 + 1)) / ((double)(m + 1));
            site.pointwise_pvalue_a2 = ((double)(site.r_a2 + 1)) / ((double)(m + 1));
        }
*/

        for (int s = 0; s < homoplasically_informative_sites.size(); s++) {
            Homoplasy_Events site = homoplasically_informative_sites.get(s);
            Binomial_Test_Stat resampled_test_stat = resampled_test_statistics[s];

            // TODO:  best practices for double treatment of integer division?
            if (site.a1_in_use) {
                resampled_test_stat.pointwise_pvalue_a1 = ((double) (resampled_test_stat.r_a1 + 1)) / ((double) (m + 1));
            }

            if (site.a2_in_use) {
                resampled_test_stat.pointwise_pvalue_a2 = ((double) (resampled_test_stat.r_a2 + 1)) / ((double) (m + 1));
            }
        }

        System.out.printf("point-wise estimates complete.%n");
    }


    /**
     * Currently, this version of the phyC test stat is meant to be close to original phyC.
     * <p>
     * For all sites in this current replicate/resample:  compares obs vs resampled homoplasy counts, update r (for both alleles).  Specifically, once for each allele:
     * phyC test stat :=  r++ if (resampled case count >= obs case count && resampled control count <= obs control count), else r is unchanged
     * <p>
     * homoplasy counts := int[4] := a1p1, a1p2, a2p1, a2p2 := {a1/case, a1/control}, {a2/case, a2/control}
     * TODO phyC:  should I enforce pheno labeling?  Is p1 always case and p2 always control?
     *
     * @param homoplasically_informative_sites
     * @param obs_homoplasy_counts
     * @param resampled_homoplasy_counts
     * @param resampled_test_statistics
     */
    private void phyC_test_stat(ArrayList<Homoplasy_Events> homoplasically_informative_sites, ArrayList<int[]> obs_homoplasy_counts,
                                ArrayList<int[]> resampled_homoplasy_counts, phyC_Test_Statistic[] resampled_test_statistics) {
        // DEBUG they should all equal each other
//        System.out.println("obs_homoplasy_counts.size() = " + obs_homoplasy_counts.size());
//        System.out.println("resampled_homoplasy_counts.size() = " + resampled_homoplasy_counts.size());
//        System.out.println("resampled_test_statistics.length = " + resampled_test_statistics.length);

        for (int s = 0; s < homoplasically_informative_sites.size(); s++) { // for each site
            Homoplasy_Events curr_site = homoplasically_informative_sites.get(s);
            phyC_Test_Statistic curr_test_stat = resampled_test_statistics[s];

            int[] curr_site_obs_counts = obs_homoplasy_counts.get(s);
            int[] curr_site_res_counts = resampled_homoplasy_counts.get(s);

            // TODO phyC:  should the comparison be a function in phyC_Test_Statistic?  eg:  curr_test_stat.is_more_extreme() or .compare()
            // if resampled homoplasy counts is as or more extreme than obs homoplasy counts:
            // a1
            if (curr_site.a1_in_use && ((curr_site_res_counts[0] >= curr_site_obs_counts[0]) && (curr_site_res_counts[1] <= curr_site_obs_counts[1]))) {
                curr_test_stat.update_r_a1();
            }

            // a2
            if (curr_site.a2_in_use && ((curr_site_res_counts[2] >= curr_site_obs_counts[2]) && (curr_site_res_counts[3] <= curr_site_obs_counts[3]))) {
                curr_test_stat.update_r_a2();
            }
        }
    }


    /**
     * Permute the phenos relative to the strains.
     *
     * @param phenos (strain id, pheno)
     * @return
     */
    private HashMap<String, String> permute_phenos(HashMap<String, String> phenos) {

        // DEBUG
//        System.out.println("phenos = " + phenos);

        // use entrySet() to guarantee the same pairing of keys and values (vs using values() and keySet())
        Set<Map.Entry<String, String>> pheno_mappings = phenos.entrySet();

        ArrayList<String> strain_ids = new ArrayList<>();
        ArrayList<String> pheno_labels = new ArrayList<>();

        for (Map.Entry<String, String> mapping : pheno_mappings) {

            // keep original key order
            strain_ids.add(mapping.getKey());

            // get original pheno order
            pheno_labels.add(mapping.getValue());
        }

        // permute phenos
        Collections.shuffle(pheno_labels);

        // reassociate phenos (now newly permuted) with strain ids (which are in original key ordering)
        HashMap<String, String> permuted_phenos = new HashMap<>();
        for (int i = 0; i < strain_ids.size(); i++) {
            permuted_phenos.put(strain_ids.get(i), pheno_labels.get(i));
        }

        // DEBUG
//        System.out.println("permuted_phenos = " + permuted_phenos);

        return permuted_phenos;
    }


    /**
     * This method simply returns all sites that are "homoplasically segregating" := both alleles show homoplasy events.
     * This method was originally created so that 2x2 association tests can only operate on homoplasically segregating sites.
     */
    private ArrayList<Homoplasy_Events> get_only_homoplasically_segregating_sites(ArrayList<Homoplasy_Events> all_events) {
        ArrayList<Homoplasy_Events> homoplasically_segregating_sites = new ArrayList<>();

        for (Homoplasy_Events segsite : all_events) {
            if (segsite.get_a1_count() >= 2 && segsite.get_a2_count() >= 2) {
                homoplasically_segregating_sites.add(segsite);
            }
        }

        // DEBUG
        System.out.println("get_only_homoplasically_segregating_sites()-------------------");
        System.out.println("all_events.size() = " + all_events.size());
        System.out.println("homoplasically_segregating_sites.size() = " + homoplasically_segregating_sites.size());

        return homoplasically_segregating_sites;
    }



/*
    private ArrayList<Fishers_Exact_Statistic> assoc_fishers(ArrayList<Homoplasy_Events> all_events, HashMap<String,String> phenos) {

        // Fisher's exact:
        // only associate homoplasically non-monomorphic alleles (i.e. 0, 1, or 2 alleles may be homoplasically segregating)
        // here we also define "derived" alleles relative to the MRCA (not some outgroup)
        for (Homoplasy_Events segsite : all_events) {
            if (segsite.get_a1_count() >= 2 && segsite.get_a2_count() >= 2) {

            }
        }

        return null;
    }
*/


    /**
     * // This method is now DEPRECATED
     *
     * Only multi-hypo corrections likely to be used are q-values and maxT, but need to test ad-hoc "sub-family" idea which should go here
     *
     *
     * Assess statistical significance of each segregating site.  Any statistical test can be used for this assessment.  Currently, using:
     * - fisher's exact test (no multi-hypo correction, I could dervie q-values from distribution of fisher's exact p-values)
     * - maxT resampling (FWER correction)
     */
/*
    private void assess_significance(ArrayList<Homoplasy_Events> all_events, HashMap<String, String> phenos) {
        ArrayList<Fishers_Exact_Statistic> test_statistics = fishers_exact(all_events, phenos);

        // chisq
        // maxT
        // other methods to assess statistical significance at each segregating site

        // pvalues
        double[] pvalues = new double[test_statistics.size()];
        for (int i = 0; i < pvalues.length; i++) {
            pvalues[i] = test_statistics.get(i).getPvalue();
        }

         // DEBUG
//        for (double p : pvalues) {
//            System.out.println(p);
//        }

        // multi-hypo correction (for right now, always calculate qvalues, later may include option to choose a different method like a FWER)
        QSet qset = qvalues_R(pvalues);

        // results:  each row := 1 seg site
        output_qvalues(all_events, test_statistics, qset);

    }
*/


    /**
     * Counts geno-pheno state pairs, and interfaces w R process to calculate fisher's exact (two-sided) for all sites.
     */
    private ArrayList<Fishers_Exact_Statistic> fishers_exact(ArrayList<Homoplasy_Events> homoplasically_segregating_sites, HashMap<String, String> phenos) {

        // TODO:  refactor count_phenos() into 2 methods:
        /*
        for (homoplasically segregating sites) {
            count_homoplasies()
            build_2x2()
        }
         */

        ArrayList<int[][]> contingency_tables = count_phenos(homoplasically_segregating_sites, phenos);
        ArrayList<Fishers_Exact_Statistic> test_statistics = fishers_exact_R(contingency_tables);

        return test_statistics;
    }


    /**
     * Builds a 2x2 contingency table for each segregating site.  These tables can be used for such tests as fishers exact, chisq, barnards exact, etc.
     * // TODO:  rename method:  count_homoplasies()?
     * // TODO:  consider breaking up this method in two:  1) count_homoplasies() returns counts[4]   2) constructing 2x2
     */
    private ArrayList<int[][]> count_phenos(ArrayList<Homoplasy_Events> homoplasically_segregating_sites, HashMap<String, String> phenos) {
        ArrayList<int[][]> contingency_tables = new ArrayList<>();

        for (Homoplasy_Events segsite : homoplasically_segregating_sites) {

            // for each homoplasic node, update one of four counters of the contingency table for this site:
            int a1_p1 = 0;  // p1 := case
            int a1_p2 = 0;  // p2 := control
            int a2_p1 = 0;
            int a2_p2 = 0;

            // TODO:  refactor the two for loops below with something like:
//            update_counters(HashMap, c1, c2)

            // allele 1
            for (Map.Entry<String, HashSet<String>> a1_MRCAs : segsite.allele1_MRCAs.entrySet()) {
                String MRCA_node_name = a1_MRCAs.getKey();

                // if no pheno is available for this node (i.e. this is an internal node), simply continue without updating any counters of the contingency table
                if (phenos.containsKey(MRCA_node_name)) {
                    String pheno = phenos.get(MRCA_node_name);

                    if (pheno.equals("1")) {  // case
                        a1_p1++;
                    } else if (pheno.equals("0")) {
                        a1_p2++;
                    }
                }
            }

            // allele 2
            for (Map.Entry<String, HashSet<String>> a2_MRCAs : segsite.allele2_MRCAs.entrySet()) {
                String MRCA_node_name = a2_MRCAs.getKey();

                if (phenos.containsKey(MRCA_node_name)) {
                    String pheno = phenos.get(MRCA_node_name);

                    if (pheno.equals("1")) {  // case
                        a2_p1++;
                    } else if (pheno.equals("0")) {
                        a2_p2++;
                    }
                }
            }

            // TODO:  set all 4 obs homoplasy counts inside this Homoplasy_Events obj, here . . . obs_counts[0] = a1_p1, etc


            int[][] contingency_table = new int[2][2];
            contingency_table[0][0] = a1_p1;             // a1p1 := a1 case
            contingency_table[0][1] = a1_p2;             // a1p2 := a1 control
            contingency_table[1][0] = a2_p1;             // a2p1 := a2 case
            contingency_table[1][1] = a2_p2;             // a2p2 := a2 control
            contingency_tables.add(contingency_table);
        }

        return contingency_tables;
    }


    /**
     * Process-to-process communication (java <-> Rscript)
     */
    private ArrayList<Fishers_Exact_Statistic> fishers_exact_R(ArrayList<int[][]> contingency_tables) {
        ArrayList<Fishers_Exact_Statistic> test_statistics = new ArrayList<>();  // TODO:  pull out to fishers_exact()?

//        ProcessBuilder pb = new ProcessBuilder("Rscript", "--vanilla", "--verbose", "piping_test.R");
        ProcessBuilder pb = new ProcessBuilder("Rscript", "--vanilla", "--verbose", R_dir + "fishers_exact.R");

        // DEBUG
//        pb.redirectErrorStream(true);  // combine stderr + stdout

        try {
            Process R_process = pb.start();

            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(R_process.getOutputStream()));

            // DEBUG
//            BufferedWriter bw = new BufferedWriter(new FileWriter("contingency_tables_as_rows.txt"));

/*
            // 3 fake contingency tables for testing:
            bw.write("10,10,10,10");
            bw.newLine();
            bw.write("60,0,10,60");
            bw.newLine();
            bw.write("50,8,5,49");
            bw.newLine();
*/


// COMMENTED OUT to allow a few fake contingency tables to be tested.  Code below is to process all seg sites:

            // pass to R process all contingency tables:  1 row representing 1 contingency table representing 1 segregating site
//            int cap = 1000000;  // DEBUG
//            int sites = 0;  // DEBUG
            for (int[][] curr_table : contingency_tables) {
                // TODO:  replace with StringBuilder?
                // row := a1p1, a1p2, a2p1, a2p2
                String table_as_row = String.valueOf(curr_table[0][0]) + "," +  // a1p1 := a1 case
                        String.valueOf(curr_table[0][1]) + "," +  // a1p2 := a1 control
                        String.valueOf(curr_table[1][0]) + "," +  // a2p1 := a2 case
                        String.valueOf(curr_table[1][1]);         // a2p2 := a2 control

                bw.write(table_as_row);
                bw.newLine();

                // DEBUG
//                sites++;
//                if (sites > cap) {
//                    break;
//                }
            }

            bw.close();

            // DEBUG
//            System.out.println("cap = " + cap);
            System.out.println("contingency_tables.size() = " + contingency_tables.size());

            // is waitFor() necessary?  Or go straight to blocking BufferedReader that waits on R process?  Ok, it's not necessary.
//            int err_code = R_process.waitFor();
//            boolean err_code = R_process.waitFor(1, TimeUnit.SECONDS);
//            System.out.println("err_code = " + err_code + "\n\n");


            // extracting results from R process:
            System.out.println("extracting results from R process");
            BufferedReader br = new BufferedReader(new InputStreamReader(R_process.getInputStream()));
            // DEBUG
            System.out.println("successfully opened input stream to R process");
            System.out.println("br.read() = " + br.ready());
            String R_process_output;
            String currline;
            StringBuilder sb = new StringBuilder(100000);  // set initial capacity to 100K chars to avoid fewer array expansions (expands at 2*current capacity)
            while ((currline = br.readLine()) != null) {  // read in ALL R output
//                System.out.println("R_process_output = " + R_process_output + "end_output");
//                R_process_output += currline + " ";
                sb.append(currline).append(" ");

                // DEBUG
//                System.out.println("currline = " + currline);
            }
            R_process_output = sb.toString();

            // parse R output
            String[] R_output = R_process_output.split(" ");

            // DEBUG
//            System.out.println("R_process_output = " + R_process_output);
            System.out.println("R_output.length = " + R_output.length);
//            for (String value : R_output) {
//                System.out.println("value = " + value);
//            }

            // DEBUG
//            int tot_num_sites = cap;  // for testing:  use 3 to test 3 contingency tables)
            int tot_num_sites = contingency_tables.size();  // i.e. tot # seg sites

            // inside R_output[]:
            // first s elements are pvalues (where s = tot # sites)
            // next s elements are ORs
            // next 2*s elements are the pairs of CIs for the ORs
            for (int pvalue_index = 0, or_index = tot_num_sites, ci_index = tot_num_sites * 2;
                 pvalue_index < tot_num_sites;
                 pvalue_index++, or_index++, ci_index += 2) {

                // add raw pheno counts as a row (same row passed into R).  contingency_tables and all_events are both in the same segsite order.
                int[][] curr_table = contingency_tables.get(pvalue_index);
                String counts_as_row = String.valueOf(curr_table[0][0]) + "," +
                        String.valueOf(curr_table[0][1]) + "," +
                        String.valueOf(curr_table[1][0]) + "," +
                        String.valueOf(curr_table[1][1]);

                // parse out pvalues from R output
                double p = Double.parseDouble(R_output[pvalue_index]);
                String or = R_output[or_index];
//                curr_stat.or_ci = R_output[j + tot_num_sites * 2] + "-" + R_output[j + tot_num_sites * 2 + 1];
                String or_ci = R_output[ci_index] + "-" + R_output[ci_index + 1];

                Fishers_Exact_Statistic curr_stat = new Fishers_Exact_Statistic(counts_as_row, p, or, or_ci);
                test_statistics.add(curr_stat);
            }

//            System.out.println("\n\nEnd Test.\n\n");
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        return test_statistics;
    }


    private QSet qvalues(ArrayList<Fishers_Exact_Statistic> test_statistics) {

        // get pvalues
        double[] pvalues = new double[test_statistics.size()];
        for (int i = 0; i < pvalues.length; i++) {
            pvalues[i] = test_statistics.get(i).getPvalue();
        }

        QSet qset = qvalues_R(pvalues);

//        output_qvalues(all_events, test_statistics, qset);

        return qset;
    }


    /**
     * This method interfaces with R process to calculate qvalue statistics for all sites.
     */
    private QSet qvalues_R(double[] pvalues) {
        // DEBUG
//        System.out.println("qvalues_R()-------------");
//        System.out.println("pvalues.length = " + pvalues.length);
//        for (double p : pvalues) {
//            System.out.print(p + "  ");
//        }


        QSet qset = null;

        ProcessBuilder pb = new ProcessBuilder("Rscript", "--vanilla", "--verbose", R_dir + "qvalues.R");

        try {
            Process R_process = pb.start();

            // write pvalues into pipe -> R process
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(R_process.getOutputStream()));
            for (double p : pvalues) {
                bw.write(String.valueOf(p));
                bw.newLine();
            }
            bw.close();

            // read results from pipe <- R process.  Parse out relevant variables from return value of Storey's qvalue() function.
            BufferedReader br = new BufferedReader(new InputStreamReader(R_process.getInputStream()));
            String R_process_output;
            String currline;
            StringBuilder sb = new StringBuilder(100000);  // set initial capacity to 100K chars to avoid fewer array expansions (expands at 2*current capacity)
            while ((currline = br.readLine()) != null) {  // read in ALL R output
//                System.out.println("R_process_output = " + R_process_output + "end_output");
//                R_process_output += currline + " ";
                sb.append(currline).append(" ");

                // DEBUG
//                System.out.println("currline = " + currline);
            }
            R_process_output = sb.toString();

            // parse R output
            String[] R_output = R_process_output.split(" ");

            // DEBUG
//            System.out.println("R_process_output = " + R_process_output + "end");
//            System.out.println("R_output.length = " + R_output.length);
//            for (String value : R_output) {
//                System.out.println("value = " + value);
//            }

            // DEBUG
//            int tot_num_sites = cap;  // for testing:  use 3 to test 3 contingency tables)

            // inside R_output[]:
            // 1st element := pi0
            // next s elements are qvalues (where s = tot # sites)
            // next s elements are lfdrs
            String pi_0 = R_output[0];
            ArrayList<String> qvalues = new ArrayList<>();
            ArrayList<String> lfdrs = new ArrayList<>();

            int tot_num_sites = pvalues.length;
            for (int qvalue_index = 1, lfdr_index = tot_num_sites + 1; qvalue_index < tot_num_sites + 1; qvalue_index++, lfdr_index++) {

                // DEBUG
//                System.out.println("tot_num_sites = " + tot_num_sites);
//                System.out.println("qvalue_index = " + qvalue_index);
//                System.out.println("lfdr_index = " + lfdr_index);
//                System.exit(-1);


                String qvalue = R_output[qvalue_index];
                String lfdr = R_output[lfdr_index];

                qvalues.add(qvalue);
                lfdrs.add(lfdr);
            }

            // DEBUG
//            System.out.println("qvalues.size() = " + qvalues.size());
//            System.out.println("lfdrs.size() = " + lfdrs.size());

            // encapsulate all results into one QSet object
            qset = new QSet(pi_0, qvalues, lfdrs);
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        return qset;
    }


    /**
     * each row := 1 seg site (outputs both test statistic (e.g. fisher's, phyC, etc) + qvalue)
     *
     * @param all_events
     * @param test_statistics
     * @param qset
     */
    private void output_qvalues(ArrayList<Homoplasy_Events> all_events, ArrayList<Fishers_Exact_Statistic> test_statistics, QSet qset) {
        /*
        // look at .forEach() java 8 construct:  http://www.javaworld.com/article/2461744/java-language/java-language-iterating-over-collections-in-java-8.html
        for (Fishers_Exact_Statistic curr_stat : test_statistics) {
            System.out.println(curr_stat);
        }

        // so something like this:
        test_statistics.forEach(curr_stat -> System.out.println(curr_stat));  // hmm but this is still only 1 data structure at a time?  although cleaner than above.
        */

        // print out 1 result row per seg site:  segsite ID, pos, a1, a2, a1 homoplasy count, a2 homoplasy count, pvalue, or, ci, qvalue, local fdr
        // all_events and test_statistics should both be in the same order (i.e. ordered by segsite_ID)
        // sanity check (both sizes must be equal)
        System.out.println("(# homoplasically segregating sites) all_events.size() = " + all_events.size());
        System.out.println("test_statistics.size() = " + test_statistics.size());

        // header row
        System.out.println(Homoplasy_Events.COL_NAMES + "\t" + Fishers_Exact_Statistic.COL_NAMES + "\tqvalue\tlfdr");
        ArrayList<String> qvalues = qset.getQvalues();
        ArrayList<String> lfdrs = qset.getLocal_fdrs();
        for (int i = 0; i < all_events.size(); i++) {
            System.out.println(all_events.get(i) + "\t" + test_statistics.get(i) + "\t" + qvalues.get(i) + "\t" + lfdrs.get(i));
        }
    }


    /**
     * Characterize various aspects of homoplasy events, ancestral reconstruction, etc.
     */
    private void characterize(ArrayList<Homoplasy_Events> all_events, NewickTree tree) {
//        output_num_events_by_segsite(all_events);
//        identify_segsites_w_back_mutations(all_events, tree);
        identify_internal_homoplasy_events(all_events, tree);

    }


    private void output_num_events_by_segsite(ArrayList<Homoplasy_Events> all_events) {
        for (int i = 0; i < all_events.size(); i++) {
            Homoplasy_Events homoplasy_events = all_events.get(i);
            homoplasy_events.output_num_of_events();
        }
    }


    /**
     * Identifies all seg sites that contain at least one back mutation.  I define a back mutation specifially as:
     * a reversion to a previous allelic state along the IBD path from root -> leaf.  As such, there are exactly p paths where p = # leaves.  If there exists 1 path
     * that contains a back mutation, then this seg site is counted as possessing a back mutation.
     * A back mutation is identified when 2 MRCAs are found along an IBD path, and any root MRCAs are not counted because they do not represent a parent with an allele
     * change simply because the root has no parent.
     */
    private void identify_segsites_w_back_mutations(ArrayList<Homoplasy_Events> all_events, NewickTree tree) {
        // DEBUG
//        System.out.println("all_events.size() = " + all_events.size());

        NewickTreeNode root = tree.getRoot();
        List leaves = tree.getLeafNodes();

//        for (Homoplasy_Events segsite_events : all_events) {
        for (int site = 0; site < all_events.size(); site++) {  // each seg site has exactly 1 Homoplasy_Events object, i.e. all_events.size() == # seg sites
            int segsite_num = site + 1;  // make this 1-based
            Homoplasy_Events curr_site_events = all_events.get(site);

            for (int i = 0; i < leaves.size(); i++) {
                int MRCA_count = 0;  // # MRCAs in this IBD path from leaf -> root.  >= 3 := at least 1 back mutation.
                NewickTreeNode curr_leaf = (NewickTreeNode) leaves.get(i);

                for (NewickTreeNode curr_node = curr_leaf; curr_node != null; curr_node = curr_node.getParent()) {  // leaf -> root
                    if (curr_site_events.isMRCA(curr_node) && !curr_node.equals(root)) {  // if there is an allelic change in the parent of the curr_node
                        MRCA_count++;
                    }
                }

                if (MRCA_count >= 2) {  // at least 1 back mutation found
                    System.out.println(segsite_num + "\t" + MRCA_count + "\t" + curr_leaf.getId());  // 1 IBD path containing at least 1 back mutation (keep in mind
                    // there can exist multiple of these paths per seg site)
                }
            }
        }
    }


    /**
     * This method counts the # of homoplasic snps that occur on internal nodes.  These internal node stats help quantify how much information is gained by using
     * phenotypes assigned to internal nodes vs using only extant phenotypes.  Remember, genotype ancestral reconstruction is still needed to identify homoplasy events
     * on leaves (not just to identify homoplasic snps occurring on internal nodes).
     * <p>
     * For example, many homoplasy events occuring on internal nodes suggests that ancestral reconstruction of phenotypes may be worthwhile.
     * <p>
     * All MRCAs that are not leaves are considered internal homoplasy events.
     */
    private void identify_internal_homoplasy_events(ArrayList<Homoplasy_Events> all_events, NewickTree tree) {
        // create HashSet of leaves
        final HashSet<String> leaf_names = new HashSet<>();
        List<NewickTreeNode> leaves = tree.getLeafNodes();
        for (NewickTreeNode leaf : leaves) {
            leaf_names.add(leaf.getId());
        }

        // for (each site):  check if each MRCA is a leaf (both MRCA HashMaps)
        for (int site = 0; site < all_events.size(); site++) {
            Homoplasy_Events site_events = all_events.get(site);
            site_events.output_internal_events(leaf_names);
        }
    }


    /**
     * All debug method calls are from here.
     */
    private void debug(ArrayList<Homoplasy_Events> all_events, HashMap<String, char[]> seg_sites, int[] segsite_indexes) {
//        homoplaser.test_newick_parser(args);
//        homoplaser.test_tree_methods(args);

//        output_segsite_events_as_txt(all_events, 3);
//        output_segsite_events_as_branch_colors(all_events, 3383);
//        output_segsite_MRCAs_and_clades_as_branch_colors(all_events, 80715);

//        output_segsite_alleles(all_events, seg_sites, segsite_indexes, 68362);
    }


    /**
     * Output ith segsite (1-based) in plain text.
     */
    private void output_segsite_events_as_txt(ArrayList<Homoplasy_Events> all_events, int i) {
        Homoplasy_Events segsite_events = all_events.get(i - 1);
        segsite_events.output_events_as_txt();
    }


    /**
     * Output ith segsite (1-based) in itol branch color format.  This helps me quickly compare the MRCAs my code identifies
     * against the MRCAs I find visually by drawing the ancestrally reconstructed alleles using Build_Branch_Color_Annotation_Files.
     */
    private void output_segsite_events_as_branch_colors(ArrayList<Homoplasy_Events> all_events, int i) {
//        Homoplasy_Events segsite_events = all_events.get(i - 1);  // now explicitly using the index into ArrayList of Homoplasy_Events so no longer need a minus 1
        Homoplasy_Events segsite_events = all_events.get(i);

        // DEBUG
        System.out.println("segsite_events.toString()\t: " + segsite_events.toString());

        segsite_events.output_events_as_branch_colors();
    }


    /**
     * Output ith segsite (1-based) in itol branch color format.  This method not only outputs the MRCA nodes that represent the homoplasy events, but also the clades
     * that belong to each MRCA node.  The final drawing should equal the drawing based upon the alleles of the ancestral reconstruction for this seg site (output from
     * Build_Branch_Color_Annotation_File).
     */
    private void output_segsite_MRCAs_and_clades_as_branch_colors(ArrayList<Homoplasy_Events> all_events, int i) {
//        Homoplasy_Events segsite_events = all_events.get(i - 1);  // now explicitly using the index into ArrayList of Homoplasy_Events so no longer need a minus 1
        Homoplasy_Events segsite_events = all_events.get(i);

        // DEBUG
        System.out.println("segsite_events.toString()\t: " + segsite_events.toString());

        segsite_events.output_events_and_their_clades_as_branch_colors();
    }


    /**
     * Simply shows all the alleles at this ith site (1-based).
     *
     * I'm using this method to compare the alleles I pull out here by using CFML files against the alleles
     * I pull out from the map/ped world.  This check is necessary so that I make sure I am using the proper
     * physical positions from the input map file (remember, for some asinine reason CFML uses an additional
     * layer in the form of the position_cross_reference file to index into the ML_sequence.fasta file).
     */
    /*private void output_segsite_alleles(ArrayList<Homoplasy_Events> all_events, HashMap<String, char[]> seg_sites, int[] segsite_indexes, int i) {
        int segsite_index = segsite_indexes[i - 1] - 1;  // segsite_index is 0-based (while values in segsite_indexes[] from position_cross_references.txt are 1-based)
        HashMap<String, Character> mapped_snps = map_seg_site_to_tree(seg_sites, segsite_index);
        mapped_snps.forEach((k, v) -> System.out.println(k + "\t" + v));
    }*/


    /**
     * An inner class that makes it easier to group all homoplasy events together and work with them.
     * Each instance represents 1 seg site.
     * <p>
     * Currently, if allele1_MRCAs has less than 2 entries, they are removed and the HashMap has size = 0 simply because there needs to be 2 MRCA events to satisfy the
     * definition of a homoplasy event.
     * <p>
     * // TODO:  should probably rename Homoplasy_Events to Segsite_Homoplasy_Events
     */
    private class Homoplasy_Events {
        private final int segsite_ID;  // Just an ID for which segsite this is.  1-based.
        // NOTE: only non-missing genotypes can be homoplasic mutations, and as such no missing geno data is added to either a1 or a2 homoplasic mutation sets.
        private HashMap<String, HashSet<String>> allele1_MRCAs = null;  // key := node name of MRCA, value := members of this MRCA: HashSet<node names>
        private HashMap<String, HashSet<String>> allele2_MRCAs = null;
        // Straight-forward extension to tri/quad-allelic sites:  change allele MRCAs into an ArrayList to allow looping through a variable # of alleles at each seg site

        private char allele1 = 'Z';
        private char allele2 = 'Z';
        private int physical_pos;  // currently relative to some refseq, but in the future I should explore reference-free identification of variants
        private int a1_count;  // # homoplasy events for allele1
        private int a2_count;
        //        private int a1_count_internal;  // # of homoplasy events for allele1 that occurs on an internal node
//        private int a2_count_internal;
        private int a1_count_extant_only = 0;
        private int a2_count_extant_only = 0;
        private int[] obs_counts;  // int[4] counts := a1p1, a1p2, a2p1, a2p2  (same as count_phenos())

        // "maf filter".  true := a1_count (or a1_count_extant_only if EXTANT_NODES_ONLY == true) > min_hcount (default min_hcount = 0, ie: ignore 0 homoplasic mutation counts)
        private boolean a1_in_use = false;  // true := a1 at this segsite is in use
        private boolean a2_in_use = false;

        // corresponds to toString() cols
        private static final String COL_NAMES = "segsite_ID\tphysical_pos\tallele1\tallele2\ta1_count\ta2_count\ta1_count_extant_only\ta2_count_extant_only\tobs_homoplasy_counts";

        // true := this segsite is covered by at least 1 recomb tract (from any sample)
        private boolean isRecombinant = false;

        Homoplasy_Events(HashMap<String, HashSet<String>> allele1_MRCAs, HashMap<String, HashSet<String>> allele2_MRCAs, int segsite_ID, int p,
                         HashMap<String, Character> major_minor_alleles, HashSet<String> leaf_names) {
            this.allele1 = major_minor_alleles.get("a1"); // a1 := major
            this.allele2 = major_minor_alleles.get("a2"); // a2 := minor
            this.allele1_MRCAs = allele1_MRCAs;
            this.allele2_MRCAs = allele2_MRCAs;
            this.segsite_ID = segsite_ID;
            this.physical_pos = p;
            this.a1_count = allele1_MRCAs.size();
            this.a2_count = allele2_MRCAs.size();

//            calc_num_extant_homoplasic_mutations(leaves);
            calc_num_extant_homoplasic_mutations(leaf_names);
        }


        @Override
        public String toString() {
            return segsite_ID + "\t" + physical_pos + "\t" + allele1 + "\t" + allele2 + "\t" + a1_count + "\t" + a2_count + "\t" +
                    a1_count_extant_only + "\t" + a2_count_extant_only + "\t" + Arrays.toString(obs_counts);
        }


        /**
         * Set the counts for extant homoplasic mutations only (for all alleles).
         *
         * @param leaf_names
         */
        private void calc_num_extant_homoplasic_mutations(HashSet<String> leaf_names) {
            // Option #1:  Most efficient way: use HashSet of leaf names, iterate through homoplasic mutations only (instead of all samples) for fewer hashes
            // for each a1 homoplasic mutation by node name
            for (String homoplasic_mutation : allele1_MRCAs.keySet()) {
                if (leaf_names.contains(homoplasic_mutation)) {
                    a1_count_extant_only++;
                }
            }

            // for each a2 homoplasic mutation by node name
            for (String homoplasic_mutation : allele2_MRCAs.keySet()) {
                if (leaf_names.contains(homoplasic_mutation)) {
                    a2_count_extant_only++;
                }
            }

/*          Option #2:  A slightly more efficient way compared to the exhaustive search below.

            // for each a1 homoplasic mutation by node name
            for (String homoplasic_mutation : allele1_MRCAs.keySet()) {
                if (leaves.stream().map(leaf -> leaf.getId()).anyMatch(leaf_name -> homoplasic_mutation.equals(leaf_name))) {  // if (contains)
                    a1_count_extant_only++;
                }
            }

            // for each a2 homoplasic mutation by node name, is it also a leaf node?
            for (String homoplasic_mutation : allele2_MRCAs.keySet()) {
                if (leaves.stream().map(leaf -> leaf.getId()).anyMatch(leaf_name -> homoplasic_mutation.equals(leaf_name))) {
                    a2_count_extant_only++;
                }
            }
*/

/*          Option #3:  Another way to do this, but not as efficient as using streams above because streams uses short-circuiting, ie not all leaves are looked at.
            // for each leaf (i.e. iterates through all samples)
            for (NewickTreeNode leaf : leaves) {
                if (allele1_MRCAs.containsKey(leaf.getId())) {
                    a1_count_extant_only++;
                }
            }
*/
        }


        /**
         * output just the # of events at this segsite: e.g. "407:  7   3"
         */
        void output_num_of_events() {
            int allele1_num_homoplasy_events = 0;
            int allele2_num_homoplasy_events = 0;

            allele1_num_homoplasy_events = allele1_MRCAs.size();
            allele2_num_homoplasy_events = allele2_MRCAs.size();

//        System.out.print("# homoplasy events (allele1) = " + allele1_num_homoplasy_events);
//        System.out.println("\t# homoplasy events (allele2) = " + allele2_num_homoplasy_events);
            System.out.println(segsite_ID + ":" + "\t" + allele1_num_homoplasy_events + "\t" + allele2_num_homoplasy_events);

        }


        /**
         * Output all homoplasy events at this segsite in plain txt (just node names for now).
         */
        void output_events_as_txt() {
            System.out.println("segsite " + segsite_ID);
            System.out.println("# homoplasy events for allele1 = " + allele1_MRCAs.size());
            System.out.println("allele1 homoplasy events (i.e. MRCA node names):");
            for (Map.Entry<String, HashSet<String>> mapping : allele1_MRCAs.entrySet()) {  // for (each MRCA for allele1)
                String MRCA_nodename = mapping.getKey();
                System.out.println(MRCA_nodename);
            }

            System.out.printf("%n# homoplasy events for allele2 = " + allele2_MRCAs.size() + "%n");
            System.out.println("allele2 homoplasy events (i.e. MRCA node names):");
            for (Map.Entry<String, HashSet<String>> mapping : allele2_MRCAs.entrySet()) {
                String MRCA_nodename = mapping.getKey();
                System.out.println(MRCA_nodename);
            }
        }


        /**
         * Outputs all homoplasy events at this segsite as itol branch colors txt format.  For now, just output the MRCA nodes (will draw children later).
         */
        void output_events_as_branch_colors() {
            final String RED = "#ff0000";
            final String GREEN = "#00ff00";

            // write headers
            System.out.println("TREE_COLORS");
            System.out.println("SEPARATOR TAB");
            System.out.println("DATA");

            // allele1 MRCAs
            for (Map.Entry<String, HashSet<String>> mapping : allele1_MRCAs.entrySet()) {  // for (each MRCA for allele1)
                String MRCA_nodename = mapping.getKey();
                System.out.print(MRCA_nodename);
                System.out.println("\tbranch" + "\t" + RED + "\t" + "normal\t4");
            }

            // allele2 MRCAs
            for (Map.Entry<String, HashSet<String>> mapping : allele2_MRCAs.entrySet()) {
                String MRCA_nodename = mapping.getKey();
                System.out.print(MRCA_nodename);
                System.out.println("\tbranch" + "\t" + GREEN + "\t" + "normal\t4");
            }
        }


        /**
         * The same as output_events_as_branch_colors() except that this method also draws the clade members of each MRCA.
         */
        void output_events_and_their_clades_as_branch_colors() {
            final String RED = "#ff0000";
            final String GREEN = "#00ff00";

            // write headers
            System.out.println("TREE_COLORS");
            System.out.println("SEPARATOR TAB");
            System.out.println("DATA");

            // allele1 MRCAs
            for (Map.Entry<String, HashSet<String>> mapping : allele1_MRCAs.entrySet()) {  // for (each MRCA for allele1)
                String MRCA_nodename = mapping.getKey();
                System.out.print(MRCA_nodename);
                System.out.println("\tbranch" + "\t" + RED + "\t" + "normal\t4");

                // draw the clade belonging to this MRCA
                HashSet<String> clade = mapping.getValue();
                for (String clade_member : clade) {
                    System.out.print(clade_member);
                    System.out.println("\tbranch" + "\t" + RED + "\t" + "normal\t4");
                }
            }

            // allele2 MRCAs
            for (Map.Entry<String, HashSet<String>> mapping : allele2_MRCAs.entrySet()) {
                String MRCA_nodename = mapping.getKey();
                System.out.print(MRCA_nodename);
                System.out.println("\tbranch" + "\t" + GREEN + "\t" + "normal\t4");

                // draw the clade belonging to this MRCA
                HashSet<String> clade = mapping.getValue();
                for (String clade_member : clade) {
                    System.out.print(clade_member);
                    System.out.println("\tbranch" + "\t" + GREEN + "\t" + "normal\t4");
                }
            }
        }


        /**
         * return true if the curr_node is a MRCA for either allele.
         */
        public boolean isMRCA(NewickTreeNode curr_node) {
            return (allele1_MRCAs.containsKey(curr_node.getId()) || allele2_MRCAs.containsKey(curr_node.getId()));
        }


        /**
         * Each row represents 1 internal homoplasy event for this seg site.
         */
        public void output_internal_events(HashSet<String> leaf_names) {

            // allele1 MRCAs
            for (Map.Entry<String, HashSet<String>> mapping : allele1_MRCAs.entrySet()) {  // for (each MRCA for allele1)
                String MRCA_nodename = mapping.getKey();

                if (!leaf_names.contains(MRCA_nodename)) {  // if (an internal node)
                    System.out.println(segsite_ID + "\t" + MRCA_nodename);
                }
            }

            // allele2 MRCAs
            for (Map.Entry<String, HashSet<String>> mapping : allele2_MRCAs.entrySet()) {
                String MRCA_nodename = mapping.getKey();

                if (!leaf_names.contains(MRCA_nodename)) {  // if (an internal node)
                    System.out.println(segsite_ID + "\t" + MRCA_nodename);
                }
            }
        }


        /**
         * For this seg site, count the # of homoplasy events occuring on internal nodes
         */
        public int get_num_internal_events(NewickTree t) {
            int num_internal_events = 0;

            return num_internal_events;
        }

        public int get_a1_count() {
            return a1_count;
        }

        public int get_a2_count() {
            return a2_count;
        }
    }


    private class Fishers_Exact_Statistic {
        private String counts;  // "a1p1, a1p2, a2p1, a2p2"
        private double pvalue;  // using double in case I want to sort by pvalue
        private String or;
        private String or_ci;
        private final static String COL_NAMES = "obs_counts\tpvalue\tor\tor_ci";

        private Fishers_Exact_Statistic(String counts, double pvalue, String or, String or_ci) {
            this.counts = counts;
            this.pvalue = pvalue;
            this.or = or;
            this.or_ci = or_ci;
        }

        public String toString() {
            return counts + "\t" + pvalue + "\t" + or + "\t" + or_ci;
        }

        private double getPvalue() {
            return pvalue;
        }
    }


    // TODO:  best practices for inner classes
    private class phyC_Test_Statistic {

        // observed
        private float obs_dn_a1;
        private float obs_dn_a2;

        // resampled
        private int r_a1 = 0;  // # replicates as or more extreme than observed counts for allele a1 (for point-wise null dist)
        private int r_a2 = 0;
        private int r_maxT_a1 = 0;  // r for maxT null dists
        private int r_maxT_a2 = 0;
        private double pointwise_pvalue_a1;
        private double pointwise_pvalue_a2;
        private double familywise_pvalue_a1;
        private double familywise_pvalue_a2;

        private final static String COL_NAMES = "r_a1\tr_a2\tpointwise_pvalue_a1\tpointwise_pvalue_a2\tobs_dn_a1\tobs_dn_a2\tr_maxT_a1\tr_maxT_a2\tfamilywise_pvalue_a1\tfamilywise_pvalue_a2";

        public phyC_Test_Statistic(boolean a1_in_use, boolean a2_in_use) {
            if (!a1_in_use) {
                r_a1 = -1;
                r_maxT_a1 = -1;
                pointwise_pvalue_a1 = Double.NaN;
                obs_dn_a1 = Float.NaN;
                familywise_pvalue_a1 = Double.NaN;
            }

            if (!a2_in_use) {
                r_a2 = -1;
                r_maxT_a2 = -1;
                pointwise_pvalue_a2 = Double.NaN;
                obs_dn_a2 = Float.NaN;
                familywise_pvalue_a2 = Double.NaN;
            }
        }

        private void update_r_a1() {
            r_a1++;
        }

        private void update_r_a2() {
            r_a2++;
        }


        @Override
        public String toString() {
            return r_a1 + "\t" + r_a2 + "\t" + pointwise_pvalue_a1 + "\t" + pointwise_pvalue_a2 + "\t" +  // point-wise
                    obs_dn_a1 + "\t" + obs_dn_a2 + "\t" + r_maxT_a1 + "\t" + r_maxT_a2 + "\t" + familywise_pvalue_a1 + "\t" + familywise_pvalue_a2;  // family-wise
        }
    }


    /**
     * One test stat obj for each site.
     * <p>
     * Note that this obj does not actually store each replicate's resampled binom pvalues as that would be memory prohibitive, only a counter r is used.
     */
    private class Binomial_Test_Stat {

        // observed
        private double obs_binom_pvalue_a1;
        private double obs_binom_pvalue_a2;

        // resampled
        private int r_a1 = 0;  // # replicates as or more extreme than observed counts for allele a1 (for point-wise null dist)
        private int r_a2 = 0;
        private int r_maxT_a1 = 0;  // r for maxT null dists
        private int r_maxT_a2 = 0;
        private double pointwise_pvalue_a1;
        private double pointwise_pvalue_a2;
        private double familywise_pvalue_a1;
        private double familywise_pvalue_a2;

        private final static String COL_NAMES = "r_a1\tr_a2\tpointwise_pvalue_a1\tpointwise_pvalue_a2\tobs_binom_pvalue_a1\tobs_binom_pvalue_a2\tr_maxT_a1\tr_maxT_a2\tfamilywise_pvalue_a1\tfamilywise_pvalue_a2";

        public Binomial_Test_Stat(boolean a1_in_use, boolean a2_in_use) {
            if (!a1_in_use) {
                r_a1 = -1;
                r_maxT_a1 = -1;
                pointwise_pvalue_a1 = Double.NaN;
                obs_binom_pvalue_a1 = Double.NaN;
                familywise_pvalue_a1 = Double.NaN;
            }

            if (!a2_in_use) {
                r_a2 = -1;
                r_maxT_a2 = -1;
                pointwise_pvalue_a2 = Double.NaN;
                obs_binom_pvalue_a2 = Double.NaN;
                familywise_pvalue_a2 = Double.NaN;
            }
        }

        private void update_r_a1() {
            r_a1++;
        }

        private void update_r_a2() {
            r_a2++;
        }


        @Override
        public String toString() {
            return r_a1 + "\t" + r_a2 + "\t" + pointwise_pvalue_a1 + "\t" + pointwise_pvalue_a2 + "\t" +  // point-wise
                    obs_binom_pvalue_a1 + "\t" + obs_binom_pvalue_a2 + "\t" + r_maxT_a1 + "\t" + r_maxT_a2 + "\t" + familywise_pvalue_a1 + "\t" + familywise_pvalue_a2;  // family-wise
        }
    }


    /**
     * This inner class represents the resampling space used.
     */
    private class Resampling_Space {

        private ArrayList<Integer> a1_permutation_array = new ArrayList<>();
        private ArrayList<Integer> a2_permutation_array = new ArrayList<>();
        private Iterator<Integer> a1_iterator;
        private Iterator<Integer> a2_iterator;

        private double p_success_a1;  // obs tot # a1 cases / tot # a1 homoplasic mutations
        private double p_success_a2;

        // DEBUG:  how many elements of the permutation arrays are actually used?
//        int a1_array_elements_used_count = 0;
//        int a2_array_elements_used_count = 0;


        public Resampling_Space(ArrayList<int[]> obs_homoplasy_counts) {
            construct_permutation_arrays(obs_homoplasy_counts);
        }


        /**
         * permutation array := {0,1}, size = # homoplasic mutations
         *
         * @param obs_homoplasy_counts int[4] := a1p1, a1p2, a2p1, a2p2   p1 := case, p2 := control
         */
        private void construct_permutation_arrays(ArrayList<int[]> obs_homoplasy_counts) {

            int tot_obs_a1_cases = 0;
            int tot_obs_a2_cases = 0;

            for (int[] site_counts : obs_homoplasy_counts) {

                tot_obs_a1_cases += site_counts[0];
                tot_obs_a2_cases += site_counts[2];

                // a1 cases
                for (int i = 0; i < site_counts[0]; i++) {
                    a1_permutation_array.add(1);
                }

                // a1 controls
                for (int i = 0; i < site_counts[1]; i++) {
                    a1_permutation_array.add(0);
                }

                // a2 cases
                for (int i = 0; i < site_counts[2]; i++) {
                    a2_permutation_array.add(1);
                }

                // a2 controls
                for (int i = 0; i < site_counts[3]; i++) {
                    a2_permutation_array.add(0);
                }
            }

            // minimize memory footprint
            a1_permutation_array.trimToSize();
            a2_permutation_array.trimToSize();

            // set initial iterators
            a1_iterator = a1_permutation_array.iterator();
            a2_iterator = a2_permutation_array.iterator();

            // calc expected background p_success_a1 and a2 values
            p_success_a1 = (double) tot_obs_a1_cases / (double) a1_permutation_array.size();
            p_success_a2 = (double) tot_obs_a2_cases / (double) a2_permutation_array.size();

            // DEBUG:  use obs p_success or all_mutations p_success?
//            p_success_a1 = 0.379;
//            p_success_a2 = 0.379;

            // DEBUG
            System.out.println("\nDEBUG:  permutation arrays");
            System.out.println("tot_obs_a1_cases = " + tot_obs_a1_cases);
            System.out.println("tot_obs_a2_cases = " + tot_obs_a2_cases);
            System.out.println("a1_permutation_array.size() = " + a1_permutation_array.size());
            System.out.println("a2_permutation_array.size() = " + a2_permutation_array.size());
            System.out.println("obs_p_success_a1 (used in all a1 binomial tests) = " + p_success_a1);
            System.out.println("obs_p_success_a2 (used in all a2 binomial tests) = " + p_success_a2);
            System.out.println();
        }


        private void permute() {
            Collections.shuffle(a1_permutation_array);
            Collections.shuffle(a2_permutation_array);
        }


        public int next_pheno(Allele_Type allele_type) {
            int pheno = -1;

            switch (allele_type) {
                case A1:
                    if (a1_iterator.hasNext()) {
                        pheno = a1_iterator.next();

                        // DEBUG:
//                        a1_array_elements_used_count++;
                    } else {  // else no more elements left to consume in permutation array (i.e.: start of new replicate)
                        // DEBUG:  check how many elements used in the replicate
//                        System.out.println("a1_array_elements_used_count = " + a1_array_elements_used_count);
//                        a1_array_elements_used_count = 0;  // reset because new Iterator at start of new replicate
//                        System.out.println("\na1_iterator.hasNext() = " + a1_iterator.hasNext());
//                        System.out.println("shuffle array!");
//                        Collections.shuffle(a1_permutation_array);
//                        System.out.println("a1_iterator.hasNext() = " + a1_iterator.hasNext());
//                        System.out.println("shuffle array!");
//                        Collections.shuffle(a1_permutation_array);
//                        System.out.println("a1_iterator.hasNext() = " + a1_iterator.hasNext());

//                        System.out.println("Reset a1_iterator");
                        a1_iterator = a1_permutation_array.iterator();  // reset iterator for next replicate
//                        System.out.println("remaining elements in a1_iterator after iterator reset:");
//                        a1_iterator.forEachRemaining(e -> System.out.println("remaining_element: " + e));
//                        System.exit(-1);

                        pheno = a1_iterator.next();
                        // DEBUG
//                        a1_array_elements_used_count++;

                    }
                    break;

                case A2:
                    if (a2_iterator.hasNext()) {
                        pheno = a2_iterator.next();

                        // DEBUG
//                        a2_array_elements_used_count++;
                    } else {
                        // DEBUG
//                        System.out.println("a2_array_elements_used_count = " + a2_array_elements_used_count);
//                        a2_array_elements_used_count = 0;

                        a2_iterator = a2_permutation_array.iterator();
                        pheno = a2_iterator.next();
                        // DEBUG
//                        a2_array_elements_used_count++;
                    }
                    break;
            }

            return pheno;
        }
    }


    /**
     * This inner class contains relevant statistics computed by Storey's QValue package for a given set of pvalues.
     * Currently, each QSet obj represents 1 allele (i.e. a1 or a2) across all segsites.  A QSet[] represents all alleles.
     */
    private class QSet {
        private String pi_0;
        private ArrayList<String> qvalues;
        private ArrayList<String> local_fdrs;

        private QSet(String pi_0, ArrayList<String> qvalues, ArrayList<String> local_fdrs) {
            this.pi_0 = pi_0;
            this.qvalues = qvalues;
            this.local_fdrs = local_fdrs;
        }

        // TODO:  consider removing accessor methods as this is an inner class
        private ArrayList<String> getQvalues() {
            return qvalues;
        }

        private ArrayList<String> getLocal_fdrs() {
            return local_fdrs;
        }
    }


    private class Resampling_Space_Diagnostics {

        private long tot_resampled_cases_a1 = 0;
        private long tot_resampled_controls_a1 = 0;
        private long tot_resampled_cases_a2 = 0;
        private long tot_resampled_controls_a2 = 0;


        public void tally_current_replicate(ArrayList<int[]> replicate) {
            for (int[] site : replicate) {
                tot_resampled_cases_a1 += site[0];
                tot_resampled_controls_a1 += site[1];
                tot_resampled_cases_a2 += site[2];
                tot_resampled_controls_a2 += site[3];
            }
        }


        private void output_resampled_pheno_counts() {
            try {
                outputOptions.debug.write("\noutput_resampled_pheno_counts():\n");
                outputOptions.debug.write("Resampled case and control counts over all replicates:\n");
                outputOptions.debug.write("tot_resampled_cases_a1 = " + tot_resampled_cases_a1 + "\n");
                outputOptions.debug.write("tot_resampled_controls_a1 = " + tot_resampled_controls_a1 + "\n");
                outputOptions.debug.write("tot_resampled_cases_a2 = " + tot_resampled_cases_a2 + "\n");
                outputOptions.debug.write("tot_resampled_controls_a2 = " + tot_resampled_controls_a2 + "\n");
                outputOptions.debug.write("tot_resampled_a1_counts = " + (tot_resampled_cases_a1 + tot_resampled_controls_a1) + "\n");
                outputOptions.debug.write("tot_resampled_a2_counts = " + (tot_resampled_cases_a2 + tot_resampled_controls_a2) + "\n");

                float resampled_case_control_ratio_a1 = (float) tot_resampled_cases_a1 / (float) tot_resampled_controls_a1;
                float resampled_case_control_ratio_a2 = (float) tot_resampled_cases_a2 / (float) tot_resampled_controls_a2;
                outputOptions.debug.write("resampled_case_control_ratio_a1 = " + resampled_case_control_ratio_a1 + "\n");
                outputOptions.debug.write("resampled_case_control_ratio_a2 = " + resampled_case_control_ratio_a2 + "\n");

                float resampled_p_success_both_alleles = (float) (tot_resampled_cases_a1 + tot_resampled_cases_a2) / (float) (tot_resampled_cases_a1 + tot_resampled_cases_a2 + tot_resampled_controls_a1 + tot_resampled_controls_a2);
                float resampled_p_success_a1 = (float) tot_resampled_cases_a1 / (float) (tot_resampled_cases_a1 + tot_resampled_controls_a1);
                float resampled_p_success_a2 = (float) tot_resampled_cases_a2 / (float) (tot_resampled_cases_a2 + tot_resampled_controls_a2);
                outputOptions.debug.write("Resampled p(success|both alleles) := # cases / (# cases + # controls) = " + resampled_p_success_both_alleles + "\n");
                outputOptions.debug.write("Resampled p(success|a1 only) = " + resampled_p_success_a1 + "\n");
                outputOptions.debug.write("Resampled p(success|a2 only) = " + resampled_p_success_a2 + "\n");
                outputOptions.debug.flush();
            } catch (IOException e) {
                System.out.printf(Ansi.AUTO.string("@|red,bold %nAn error has occurred while writing out to the debug file.  Please see the stack trace in the log file for more information on this error.%n|@"));
                e.printStackTrace();
                System.exit(-1);
            }
        }


        //  replicate is an ArrayList of sites counts, each element represents 1 site
        public void output_obs_p_success(ArrayList<int[]> replicate) {
            int a1_case_counts = 0;
            int a1_control_counts = 0;
            int a2_case_counts = 0;
            int a2_control_counts = 0;

            for (int[] site : replicate) {  // tally up all sites
                a1_case_counts += site[0];
                a1_control_counts += site[1];
                a2_case_counts += site[2];
                a2_control_counts += site[3];
            }

            float obs_p_success_both_alleles = (float) (a1_case_counts + a2_case_counts) / (float) (a1_case_counts + a2_case_counts + a1_control_counts + a2_control_counts);
            System.out.println(obs_p_success_both_alleles);
        }
    }



    private class Row {
        private Homoplasy_Events homoplasy_events_for_one_site;
        private Binomial_Test_Stat test_stat_for_one_site;

        public Row(Homoplasy_Events homoplasy_events_for_one_site, Binomial_Test_Stat test_stat_for_one_site) {
            this.homoplasy_events_for_one_site = homoplasy_events_for_one_site;
            this.test_stat_for_one_site = test_stat_for_one_site;
        }

//        public Double get_maxT() {
//            return test_stat_for_one_site.familywise_pvalue_a2;
//        }
    }


    // --Commented out by Inspection START (9/17/16, 8:20 PM):
//    /**
//     * Stats for 1 seg site.
//     */
//    private void gather_stats(HashMap<String, HashSet<String>> allele1_mrcAs, HashMap<String, HashSet<String>> allele2_mrcAs) {
//        int allele1_num_homoplasy_events = 0;
//        int allele2_num_homoplasy_events = 0;
//
//        if (allele1_mrcAs.size() > 1) {
//            allele1_num_homoplasy_events = allele1_mrcAs.size();
//        }
//
//        if (allele2_mrcAs.size() > 1) {
//            allele2_num_homoplasy_events = allele2_mrcAs.size();
//        }
//
////        System.out.print("# homoplasy events (allele1) = " + allele1_num_homoplasy_events);
////        System.out.println("\t# homoplasy events (allele2) = " + allele2_num_homoplasy_events);
//        System.out.println(allele1_num_homoplasy_events + "\t" + allele2_num_homoplasy_events);
//    }
// --Commented out by Inspection STOP (9/17/16, 8:20 PM)


// --Commented out by Inspection START (9/13/16, 6:22 PM):
//    /**
//     * Testing the coevolution API only.
//     */
//    private void test_newick_parser(String[] args) {
//
//        try (BufferedReader br = new BufferedReader(new FileReader(args[0]))) {
//            // read in a newick tree
//            NewickTreeTokenizer tokenizer = new NewickTreeTokenizer(br);
//            NewickTreeReader tree_reader = new NewickTreeReader(tokenizer);
//            NewickTree tree = tree_reader.readTree();
//
//            // spit the tree back out in newick txt format
//            String newick_as_string = tree.toString();
//            System.out.println(newick_as_string);
//        } catch (IOException | DataFormatException e) {
//            e.printStackTrace();
//        }
//    }
// --Commented out by Inspection STOP (9/13/16, 6:22 PM)


// --Commented out by Inspection START (9/13/16, 6:23 PM):
//    /**
//     * Testing the coevolution API only.
//     */
//    private void test_tree_methods(String[] args) {
//        NewickTree tree = null;
//
//        try (BufferedReader br = new BufferedReader(new FileReader(args[0]))) {
//            // read in a newick tree
//            NewickTreeTokenizer tokenizer = new NewickTreeTokenizer(br);
//            NewickTreeReader tree_reader = new NewickTreeReader(tokenizer);
//            tree = tree_reader.readTree();
//        } catch (IOException | DataFormatException e) {
//            e.printStackTrace();
//            System.exit(-1);
//        }
//
//        // test some tree methods:
//        NewickTreeNode root = tree.getRoot();
//        System.out.println("\nroot: " + root.getId());
//
//        // verify # of children of root node:
//        System.out.println("\nverifying # children of root node:");
//        System.out.println("root.getChildrenCount() = " + root.getChildrenCount());
//        ArrayList children = (ArrayList) root.getChildren();
//        System.out.println("children.size() = " + children.size());
//
//        // verify children node names:
//        System.out.println("\nverifying node names of children of root directly with ArrayList:");
//        for (int i = 0; i < children.size(); i++) {
//            NewickTreeNode curr_child_node =  (NewickTreeNode) children.get(i);
//            System.out.println("curr_child_node.getId() = " + curr_child_node.getId());
//        }
//
//        // verify children node names using API's iterator method:
//        System.out.println("\nverifying node names of children of root node using API's iterator:");
//        Iterator it = root.getChildrenIterator();
//        while (it.hasNext()) {
//            NewickTreeNode curr_child_node = (NewickTreeNode) it.next();
//            System.out.println("curr_child_node.getId() = " + curr_child_node.getId());
//        }
//
//        // verify leaf method
//        System.out.println("\nverifying leaf nodes:");
//        ArrayList leaves = (ArrayList) tree.getLeafNodes();
//        System.out.println("leaves.size() = " + leaves.size());
//        System.out.println("leaves = " + leaves);
//
//
//    }
// --Commented out by Inspection STOP (9/13/16, 6:23 PM)
}

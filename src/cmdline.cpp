#include "cmdline.hpp"

#include <args.hxx>
#include "version.hpp"

class Version {};


CommandLineOptions parse_command_line_arguments(int argc, char **argv) {

    args::ArgumentParser parser("strobelign " + version_string());
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "strobealign";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});
    args::ActionFlag version(parser, "version", "Print version and exit", {"version"}, []() { throw Version(); });

    // Threading
    args::ValueFlag<int> threads(parser, "INT", "Number of threads [3]", {'t', "threads"});
    args::ValueFlag<int> chunk_size(parser, "INT", "Number of reads processed by a worker thread at once [10000]", {"chunk-size"}, args::Options::Hidden);

    args::Group io(parser, "Input/output:");
    args::ValueFlag<std::string> o(parser, "PATH", "redirect output to file [stdout]", {'o'});
    args::Flag v(parser, "v", "Verbose output", {'v'});
    args::Flag x(parser, "x", "Only map reads, no base level alignment (produces PAF file)", {'x'});
    args::ValueFlag<std::string> rgid(parser, "ID", "Read group ID", {"rg-id"});
    args::ValueFlagList<std::string> rg(parser, "TAG:VALUE", "Add read group metadata to SAM header (can be specified multiple times). Example: SM:samplename", {"rg"});

    args::ValueFlag<int> N(parser, "INT", "Retain at most INT secondary alignments (is upper bounded by -M and depends on -S) [0]", {'N'});
    args::ValueFlag<std::string> L(parser, "PATH", "Print statistics of indexing to PATH", {'L'});
    args::Flag i(parser, "index", "Do not map reads; only generate the strobemer index and write it to disk. If read files are provided, they are used to estimate read length", {"create-index", 'i'});
    args::Flag use_index(parser, "use_index", "Use a pre-generated index previously written with --create-index.", { "use-index" });

    args::Group seeding(parser, "Seeding:");
    //args::ValueFlag<int> n(parser, "INT", "Number of strobes [2]", {'n'});
    args::ValueFlag<int> r(parser, "INT", "Mean read length. This parameter is estimated from the first 500 records in each read file. No need to set this explicitly unless you have a reason.", {'r'});
    args::ValueFlag<int> m(parser, "INT", "Maximum seed length. Defaults to r - 50. For reasonable values on -l and -u, the seed length distribution is usually determined by parameters l and u. Then, this parameter is only active in regions where syncmers are very sparse.", {'m'});

    args::ValueFlag<int> k(parser, "INT", "Strobe length, has to be below 32. [20]", {'k'});
    args::ValueFlag<int> l(parser, "INT", "Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]", {'l'});

    args::ValueFlag<int> u(parser, "INT", "Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]", {'u'});
    args::ValueFlag<int> c(parser, "INT", "Bitcount length between 2 and 63. [8]", {'c'});
    args::ValueFlag<int> s(parser, "INT", "Submer size used for creating syncmers [k-4]. Only even numbers on k-s allowed. A value of s=k-4 roughly represents w=10 as minimizer window [k-4]. It is recommended not to change this parameter unless you have a good understanding of syncmers as it will drastically change the memory usage and results with non default values.", {'s'});

    args::Group alignment(parser, "Alignment:");
    args::ValueFlag<int> A(parser, "INT", "Matching score [2]", {'A'});
    args::ValueFlag<int> B(parser, "INT", "Mismatch penalty [8]", {'B'});
    args::ValueFlag<int> O(parser, "INT", "Gap open penalty [12]", {'O'});
    args::ValueFlag<int> E(parser, "INT", "Gap extension penalty [1]", {'E'});

    args::Group search(parser, "Search parameters:");
    args::ValueFlag<float> f(parser, "FLOAT", "Top fraction of repetitive strobemers to filter out from sampling [0.0002]", {'f'});
    args::ValueFlag<float> S(parser, "FLOAT", "Try candidate sites with mapping score at least S of maximum mapping score [0.5]", {'S'});
    args::ValueFlag<int> M(parser, "INT", "Maximum number of mapping sites to try [20]", {'M'});
    args::ValueFlag<int> R(parser, "INT", "Rescue level. Perform additional search for reads with many repetitive seeds filtered out. This search includes seeds of R*repetitive_seed_size_filter (default: R=2). Higher R than default makes strobealign significantly slower but more accurate. R <= 1 deactivates rescue and is the fastest.", {'R'});

    args::Positional<std::string> ref_filename(parser, "reference", "Reference in FASTA format", args::Options::Required);
    args::Positional<std::string> reads1_filename(parser, "reads1", "Reads 1 in FASTA or FASTQ format, optionally gzip compressed");
    args::Positional<std::string> reads2_filename(parser, "reads2", "Reads 2 in FASTA or FASTQ format, optionally gzip compressed");

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e) {
        std::cout << e.what();
        exit(EXIT_SUCCESS);
    }
    catch (const args::Help&) {
        std::cout << parser;
        exit(EXIT_SUCCESS);
    }
    catch (const Version& e) {
        std::cout << version_string() << std::endl;
        exit(EXIT_SUCCESS);
    }
    catch (const args::Error& e) {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    CommandLineOptions opt;

    // Threading
    if (threads) { opt.n_threads = args::get(threads); }
    if (chunk_size) { opt.chunk_size = args::get(chunk_size); }

    // Input/output
    if (o) { opt.output_file_name = args::get(o); opt.write_to_stdout = false; }
    if (v) { opt.verbose = true; }
    if (x) { opt.is_sam_out = false; }
    if (rgid) { opt.read_group_id = args::get(rgid); }
    if (rg) { opt.read_group_fields = args::get(rg); }
    if (N) { opt.max_secondary = args::get(N); }
    if (L) { opt.logfile_name = args::get(L); }
    if (i) { opt.only_gen_index = true; }
    if (use_index) { opt.use_index = true; }

    // Seeding
    if (r) { opt.r = args::get(r); opt.r_set = true; }
    if (m) { opt.max_seed_len = args::get(m); opt.max_seed_len_set = true; }
    if (k) { opt.k = args::get(k); opt.k_set = true; }
    if (l) { opt.l = args::get(l); }
    if (u) { opt.u = args::get(u); }
    if (s) { opt.s = args::get(s); opt.s_set = true; }
    if (c) { opt.c = args::get(c); opt.c_set = true; }

    // Alignment
    // if (n) { n = args::get(n); }
    if (A) { opt.A = args::get(A); }
    if (B) { opt.B = args::get(B); }
    if (O) { opt.O = args::get(O); }
    if (E) { opt.E = args::get(E); }

    // Search parameters
    if (f) { opt.f = args::get(f); }
    if (S) { opt.dropoff_threshold = args::get(S); }
    if (M) { opt.maxTries = args::get(M); }
    if (R) { opt.R = args::get(R); }

    // Reference and read files
    opt.ref_filename = args::get(ref_filename);
    opt.reads_filename1 = args::get(reads1_filename);

    if (reads2_filename) {
        opt.reads_filename2 = args::get(reads2_filename);
        opt.is_SE = false;
    } else {
        opt.reads_filename2 = std::string();
        opt.is_SE = true;
    }

    if (opt.use_index && opt.only_gen_index) {
        std::cerr << "Error: Options -i and --use-index cannot be used at the same time" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (opt.reads_filename1.empty() && !opt.only_gen_index) {
        std::cerr << "At least one file with reads must be specified." << std::endl;
        exit(EXIT_FAILURE);
    }

    return opt;
}

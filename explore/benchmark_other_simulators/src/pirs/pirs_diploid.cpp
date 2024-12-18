#define __STDC_FORMAT_MACROS
#include <assert.h>
#include <ctime>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "Random.h"
#include "util.h"
#include "InputStream.h"
#include "OutputStream.h"

using std::string;
using std::vector;

extern bool read_scaffold(InputStream &in, string &id, vector<char> &seq);

static const char *ref_filename = NULL;

// Probability of simulating a SNP at a given position.
static double snp_rate = 0.001;

// Probability of simulating a small indel at a given position.
// small_insertion_rate and small_deletion_rate are currently each set at 1/2
// the small_indel_rate.
static double small_indel_rate = 0.0001;
static double small_insertion_rate;
static double small_deletion_rate;

// Probability of simulating a big structural variation at a given position.  
// sv_insertion_rate, sv_deletion_rate, and sv_inversion_rate are currently each
// set at 1/3 the sv_rate.
static double sv_rate = 0.000001;
static double sv_insertion_rate;
static double sv_deletion_rate;
static double sv_inversion_rate;

// Ratio of transitions to transversions in simulated SNPs.
static double snp_transition_to_transversion_ratio = 2; 

// Prefix of output files.
static string output_prefix = "pirs_diploid";

// Output file name (set to default value if not provided on the command line)
static const char *output_filename = NULL;

// Whether to write the log files or not.
static bool write_log_files = true;

// Cumulative probability distribution for SNP biases; the first entry indicates
// a transition, while the remaining two entries indicate the two possible
// transversions for a given base.
static double snp_bias[3];

// Cumulative probability distribution for the lengths of small indels, and the
// corresponding indel lengths.
// indel-length 1(64.82%), 2(17.17%), 3(7.2%), 4(7.29%), 5(2.18%), 6(1.34%)
static const size_t small_indel_len_tab[6] = {1, 2, 3, 4, 5, 6};
static const double small_indel_len_cum_prob_tab[6] = {
	0.6482, 0.8199, 0.8919, 0.9648, 0.9866, 2.0
};

// Cumulative probability distribution for the lengths of large structural
// variation sites, and the corresponding lengths.
// SV-length 100(70%),200(20%),500(7%),1000(2%),2000(1%)
static const size_t sv_len_tab[5] = {100, 200, 500, 1000, 2000};
static const double sv_len_cum_prob_tab[5] = {0.70, 0.90, 0.97, 0.99, 2.0};


// This array is initialized as a cumulative probabality array for
// heterozygosity types.
//
// The types are:
// 	{ None, SNP, small insertion, small deletion, large insertion, large
// 	deletion, large inversion}
//
// The idea is to set up a cumulative probability distribution for the
// probability of each event occuring, but have the probabilities be 64-bit
// unsigned integers such that a probability of 100% is (2**64 - 1).  Then,
// given a random number, we can search the table for the first index that
// contains a number less than or equal to the random number.  This index is the
// heterozygosity type we choose.
static uint64_t het_type_tab[7];

// Random number generator.
// pirs_diploid() is not multi-threaded so there is just one random number
// generator.
static Random rgen;

// Reference genome file
static InputStream in_file;

// Output genome file, with all the SNPs, indels, etc. made to it.
static OutputStream out_file;

// Log files for SNPs, indels, and inversions.
static OutputStream snp_log_file;
static OutputStream indel_log_file;
static OutputStream inversion_log_file;

static void pirs_diploid_usage()
{
	const char *usage_str =
"Usage: pirs diploid [OPTIONS...] REFERENCE\n"
"\n"
"Simulate a diploid genome by creating a copy of a haploid genome with\n"
"heterozygosity introduced.  REFERENCE specifies a FASTA file containing\n"
"the reference genome.  It may be compressed (gzip).  It may contain multiple\n"
"sequences (scaffolds or chromosomes), each marked with a separate FASTA tag\n"
"line.  The introduced heterozygosity takes the form of SNPs, indels, and\n"
"large-scale structural variation (insertions, deletions and inversions).\n"
"If REFERENCE is '-', the reference sequence is read from stdin, but it must be\n"
"uncompressed.\n"
"\n"
"The probabilities of SNPs, indels, and large-scale structural variation can be\n"
"specified with the '-s', '-d', and '-v' options, respectively.  You can also\n"
"set the ratio of transitions to transversions (for SNPs) with the '-R' option.\n"
"\n"
"Indels are split evenly between insertions and deletions. The length\n"
"distribution of the indels is as follows and is derived from panda\n"
"re-sequencing data:\n"
"\t1bp\t64.82%\n"
"\t2bp\t17.17%\n"
"\t3bp\t7.20%\n"
"\t4bp\t7.29%\n"
"\t5bp\t2.18%\n"
"\t6bp\t1.34%\n"
"\n"
"Large-scale structural variation is split evenly among large-scale insertions,\n"
"deletions, and inversions.  By default, the length distribution of these\n"
"large-scale features is as follows:\n"
"\t100bp\t70%\n"
"\t200bp\t20%\n"
"\t500bp\t7%\n"
"\t1000bp\t2%\n"
"\t2000bp\t1%\n"
"\n"
"`pirs diploid' does not use multiple threads, even if pIRS was configured with\n"
"--enable-multiple threads.\n"
"\n"
"OPTIONS:\n"
"  -s, --snp-rate=RATE    A floating-point number in the interval [0, 1] that\n"
"                           specifies the heterozygous SNP rate.  Default: 0.001\n"
"\n"
"  -d, --indel-rate=RATE  A floating-point number in the interval [0, 1] that\n"
"                           specifies the heterozygous indel rate.\n"
"                           Default: 0.0001\n"
"\n"
"  -v, --sv-rate=RATE     A floating-point number in the interval [0, 1] that\n"
"                         specifies the large-scale structural variation\n"
"                         (insertion, deletion, inversion) rate in the diploid\n"
"                           genome. Default: 0.000001\n"
"\n"
"  -R, --transition-to-transversion-ratio=RATIO\n"
"                         In a SNP, a transition is when a purine or pyrimidine\n"
"                           is changed to one of the same (A <=> G, C <=> T)\n"
"                           while a transversion is when a purine is changed\n"
"                           into a pyrimidine or vice-versa.  This option\n"
"                           specifies a floating-point number RATIO that gives\n"
"                           the ratio of the transition probability to the\n"
"                           transversion probability for simulated SNPs.\n"
"                           Default: 2.0\n"
"\n"
"  -o, --output-prefix=PREFIX\n"
"                         Use PREFIX as the prefix of the output file and logs.\n"
"                           Default: \"pirs_diploid\"\n"
"\n"
"  -O, --output-file=FILE\n"
"                        Use FILE as the name of the output file. Use '-'\n"
"                           for standard output; this also moves the\n"
"                           informational messages from stdout to stderr.\n"
"\n"
"  -c, --output-file-type=TYPE\n"
"                         The string \"text\" or \"gzip\" to specify the type of\n"
"                           the output FASTA file containing the diploid copy\n"
"                           of the genome, as well as the log files.\n"
"                           Default: \"text\"\n"
"\n"
"  -n, --no-logs          Do not write the log files.\n"
"\n"
"  -S, --random-seed=SEED Use SEED as the random seed. Default:\n"
"                           time(NULL) * getpid()\n"
"\n"
"  -q, --quiet            Do not print informational messages.\n"
"\n"
"  -h, --help             Show this help and exit.\n"
"\n"
"  -V, --version          Show version information and exit.\n"
"\n"
"EXAMPLE:\n"
"  ./pirs diploid ref_sequence.fa -s 0.001 -d 0.0001 -v 0.000001\\\n"
"                 -o ref_sequence >pirs.out 2>pirs.err"
	;
	puts(usage_str);
}

static void pirs_diploid_usage_short()
{
	const char *usage_str_short = 
"Usage: pirs diploid [OPTIONS...] REFERENCE\n"
"Try `pirs diploid -h' for more information.\n"
	;
	puts(usage_str_short);
}

static const char *optstr = "s:d:v:R:o:O:c:nS:qhV";
static const struct option longopts[] = {
	{"snp-rate",			required_argument, NULL, 's'},
	{"indel-rate",			required_argument, NULL, 'd'},
	{"sv-rate",			required_argument, NULL, 'v'},
	{"transition-to-transversion-ratio",required_argument, NULL, 'R'},
	{"output-prefix",		required_argument, NULL, 'o'},
	{"output-file",			required_argument, NULL, 'O'},
	{"output-file-type",		required_argument, NULL, 'c'},
	{"no-logs",			no_argument,       NULL, 'n'},
	{"random-seed",		 	required_argument, NULL, 'S'},
	{"quiet",		 	no_argument, 	   NULL, 'q'},
	{"help",			no_argument,	   NULL, 'h'},
	{"version",			no_argument,	   NULL, 'V'},
	{NULL, 0, NULL, 0}
};

/* 
 * Parses the command line given to pirs diploid program.  There is no structure
 * to hold the parameters; we just set static variables.
 */
static void pirs_diploid_parse_command_line(int argc, char *argv[])
{
	int c;
	bool have_seed = false;
	uint64_t seed = 0;
	char *tmp;
	while ((c = getopt_long(argc, argv, optstr, longopts, NULL)) != -1) {
		switch(c) {
		case 's': 
			snp_rate = strtod(optarg, &tmp);
			if (tmp == optarg || !in_unit_interval(snp_rate))
				fatal_error("The SNP rate must be a double in "
					    "the interval [0, 1]!");
			break;
		case 'd': 
			small_indel_rate = strtod(optarg, &tmp);
			if (tmp == optarg || !in_unit_interval(snp_rate))
				fatal_error("The indel rate must be a double "
					    "in the interval [0, 1]");
			break;
		case 'v': 
			sv_rate = strtod(optarg, &tmp);
			if (tmp == optarg || !in_unit_interval(sv_rate))
				fatal_error("The large structural variation rate "
					    "must be a double in the interval [0, 1]");
			break;
		case 'R': 
			snp_transition_to_transversion_ratio = 
					strtod(optarg, &tmp);
			if (tmp == optarg || snp_transition_to_transversion_ratio < 0.0)
				fatal_error("The rate of transition to transversion "
					    "cannot be negative!");
			break;
		case 'c': 
			OutputStream::set_default_output_type(optarg);
			break;
		case 'o': 
			output_prefix = optarg;
			break;
		case 'O':
			output_filename = optarg;
			if (strcmp(output_filename, "-") == 0)
				info_messages_fp = stderr;
			break;
		case 'n':
			write_log_files = false;
			break;
		case 'S': {
				seed = strtoull(optarg, &tmp, 10);
				if (tmp == optarg)
					fatal_error("The random seed must be "
						    "an integer.");
				have_seed = true;
			}
			break;
		case 'q':
			info_messages_fp = NULL;
			break;
		case 'h': 
		case 'V':
			pirs_diploid_usage();
			exit(0);
		default:
			pirs_diploid_usage_short();
			exit(2);
		}
	}
	argv += optind;
	argc -= optind;
	if (argc > 1)
		fatal_error("Too many arguments! Run `pirs -h' "
			    "for usage information.");

	if (argc == 0)
		fatal_error("Reference sequence not supplied.  Run `pirs diploid -h' "
			    "for usage information.");
	
	ref_filename = argv[0];

	if (snp_rate == 0.0 && small_indel_rate == 0.0 && sv_rate == 0.0) {
		fatal_error("The SNP rate, Indel rate, and big SV rate were all 0.\n"
			    "Please input at least one parameter greater than 0.");
	}

	if (snp_rate > 0.0) {
		double snp_sum = snp_transition_to_transversion_ratio + 1.0;
		double transition_rate = snp_transition_to_transversion_ratio / snp_sum;
		double transversion_rate = 1.0 - transition_rate;
		snp_bias[0] = transition_rate;
		snp_bias[1] = transition_rate + transversion_rate * 0.5;
		snp_bias[2] = 2.0;
	}

	if (!have_seed)
		seed = (uint64_t)time(NULL) * (uint64_t)getpid();
	rgen.seed(seed);
}


static char get_snp(char base, const double snp_bias[])
{
	static const char snp_tab[4][3] = {
		{'G', 'T', 'C'}, //A->G transition, A->T/C transversion
		{'T', 'A', 'G'}, //C->T transition, C->A/T transversion
		{'A', 'T', 'C'}, //G->A transition, G->T/C transversion
		{'C', 'G', 'A'}, //T->C transition, T->G/A transversion
	};

	double r = rgen.next_double();
	size_t i = 0;
	while (r > snp_bias[i])
		i++;
	return snp_tab[dna_char_to_bin(base)][i];
}

// Counter structure for the number of occurrences of a SNP, indel, or
// inversion, and the total number of base pairs involved in the occurrences.
class Counter {
public:
	size_t count;
	size_t total_len;
	Counter() : count(0), total_len(0) 
	{ }

	void put(size_t len)
	{
		count++;
		total_len += len;
	}
};

// 
// Try doing an indel.  If successful, the indel is logged to the
// indel_log_file, and insertions are added to the diploid sequence, while
// deletions cause the index in the reference sequence to be advanced.
//
// @indel_cum_prob_tab, @indel_len_tab
// 	Cumulative probability distribution for the indel lengths, and the
// 	corresponding lengths.
// @ref_seq, @ref_seq_idx, @ref_seq_len
// 	The reference sequence, the current index (counting from 0), and the
// 	total length of the reference sequence.  @ref_seq_idx is advanced if a
// 	deletion is done.
// @diploid_seq
// 	The second half of the diploid sequence, under construction.
// @is_deletion
// 	#true if we want to do an insertion; #false if we want to do a deletion.
// @id
// 	ID of the current scaffold or chromosome, from the input FASTA file.
// @counter
// 	Counter object to count the number of times an insertion or deletion has
// 	occurred.
// 
// @return
// 	#true if an indel was done; #false if it was not, due to the desired
// 	length being too large for the remaining reference sequence, or due to
// 	the needed part of the reference sequence containing non-DNA characters.
//
static bool try_indel(const double indel_cum_prob_tab[],
		      const size_t indel_len_tab[],
		      const vector<char> &ref_seq, size_t &ref_seq_idx,
		      vector<char> &diploid_seq,
		      bool is_deletion, const char *id, Counter &counter)
{
	double r = rgen.next_double();
	size_t i = 0;

	while (r > indel_cum_prob_tab[i])
		i++;

	size_t indel_len = indel_len_tab[i];
	if (is_deletion && ref_seq_idx + indel_len > ref_seq.size())
		return false;

	char insert[indel_len];
	const char *indel;
	uint64_t log_ref_seq_idx;
	uint64_t log_diploid_seq_idx;
	if (is_deletion) {
		// Log reference index as the 1-based index of the first deleted
		// base.
		log_ref_seq_idx = ref_seq_idx + 1;

		// Log diploid index as the 1-based index of the diploid base
		// corresponding to the reference base directly before the
		// deletion (which is the previous base added to the diploid
		// sequence)
		log_diploid_seq_idx = diploid_seq.size();
		indel = &ref_seq[ref_seq_idx];
		ref_seq_idx += indel_len;
	} else {
		// Log diploid index as the 1-based index of the first inserted
		// base.
		log_diploid_seq_idx = diploid_seq.size() + 1;
		
		// Log reference index as the 1-based index of the reference
		// base corresponding to the diploid base directly before the
		// insertion (which is the previous base in the reference
		// sequence)
		log_ref_seq_idx = ref_seq_idx;

		rgen.random_dna_seq(insert, indel_len);
		indel = insert;
		diploid_seq.insert(diploid_seq.end(), &indel[0], &indel[indel_len]);
	}

	if (indel_log_file.is_open()) {
		indel_log_file.printf("%s\t%zu\t%zu\t%c\t%zu\t",
				      id, log_ref_seq_idx, log_diploid_seq_idx,
				      (is_deletion) ? '-' : '+',
				      indel_len);
		indel_log_file.write(indel, indel_len);
		indel_log_file.putc('\n');
	}
	counter.put(indel_len);
	return true;
}

static bool try_inversion(const double inversion_len_cum_prob_tab[],
			  const size_t inversion_len_tab[], 
			  const vector<char> &ref_seq, size_t &ref_seq_idx,
			  vector<char> &diploid_seq,
			  const char *id, Counter &counter)
{
	double r = rgen.next_double();
	size_t inversion_len;
	size_t i = 0;

	while (r > inversion_len_cum_prob_tab[i])
		i++;
	inversion_len = inversion_len_tab[i];
	if (ref_seq_idx + inversion_len > ref_seq.size())
		return false;

	const char *inversion = &ref_seq[ref_seq_idx];

	if (seq_contains_non_dna_chars(inversion, inversion_len))
		return false;

	if (inversion_log_file.is_open()) {
		inversion_log_file.printf("%s\t%zu\t%zu\t%zu\t",
					  id, ref_seq_idx + 1, 
					  diploid_seq.size() + 1, inversion_len);
	}

	ref_seq_idx += inversion_len;
	counter.put(inversion_len);

	const char *p = inversion + inversion_len - 1;
	while (inversion_len--) {
		diploid_seq.push_back(dna_char_complement(*p));
		if (inversion_log_file.is_open())
			inversion_log_file.putc(*p);
		p++;
	}
	if (inversion_log_file.is_open())
		inversion_log_file.putc('\n');
	return true;
}

//
// Simulates the second half of a diploid sequence, given a reference sequence.
//
// @diploid_seq
// 	An empty char vector into which the diploid sequence will be returned.
// @ref_seq
// 	The reference scaffold or chromosome.
// @id 
// 	ID of the current scaffold or chromosome, from the input FASTA file.
static void make_diploid_seq(vector<char> &diploid_seq,
			     const vector<char> &ref_seq, const char *id)
{
	Counter snp_counter;
	Counter small_insert_counter;
	Counter small_delet_counter;
	Counter sv_insert_counter;
	Counter sv_delet_counter;
	Counter sv_invert_counter;

	size_t ref_seq_idx = 0;

	while (ref_seq_idx < ref_seq.size()) {
		char ref_base = ref_seq[ref_seq_idx];
		uint64_t r = rgen.next_uint64();

		if (r < het_type_tab[0]) {
			// No heterozygosity here
			diploid_seq.push_back(ref_base);
			ref_seq_idx++;
			continue;
		} else if (r < het_type_tab[1]) {
			// Do a SNP, unless it's an unknown base pair here.
			char new_base;
			if (is_valid_dna_char(ref_base)) {
				new_base = get_snp(ref_base, snp_bias);
				snp_counter.put(1);
				if (snp_log_file.is_open()) {
					snp_log_file.printf("%s\t%zu\t%zu\t%c\t%c\n",
							    id, ref_seq_idx + 1 ,
							    diploid_seq.size() + 1,
							    ref_base, new_base);
				}
			} else {
				new_base = ref_base;
			}
			diploid_seq.push_back(new_base);
			ref_seq_idx++;
			continue;
		} 

		// There is heterozygosity here, but not a SNP.  To make sure we
		// always are making progress, and not adding a whole bunch of
		// short indels right after each other, add a homozygous base
		// first.
		diploid_seq.push_back(ref_base);
		ref_seq_idx++;

		if (r < het_type_tab[2]) {
			// Do a small insertion.
			try_indel(small_indel_len_cum_prob_tab,
				  small_indel_len_tab,
				  ref_seq, ref_seq_idx,
				  diploid_seq, false, id, 
				  small_insert_counter);

		} else if (r < het_type_tab[3]) {
			// Try doing a small deletion.
			try_indel(small_indel_len_cum_prob_tab,
				  small_indel_len_tab,
				  ref_seq, ref_seq_idx,
				  diploid_seq, true, id, 
				  small_delet_counter);
		} else if (r < het_type_tab[4]) {
			// Do a large structural variation (SV) insertion.
			try_indel(sv_len_cum_prob_tab, sv_len_tab,
				  ref_seq, ref_seq_idx,
				  diploid_seq, false, id, 
				  sv_insert_counter);
		} else if (r < het_type_tab[5]) {
			// Try doing a large structural variation (SV) deletion.
			try_indel(sv_len_cum_prob_tab, sv_len_tab,
				  ref_seq, ref_seq_idx,
				  diploid_seq, true, id,
				  sv_delet_counter);
		} else {
			// Try doing a large structural variation (SV) inversion.
			try_inversion(sv_len_cum_prob_tab, sv_len_tab,
				      ref_seq, ref_seq_idx,
				      diploid_seq, 
				      id, sv_invert_counter);
		}
	}
	
	/* Print some statistics about the diploid genome we created. */
	info("\n");
	info("Heterozygosity introduced to sequence \"%s\":\n", id);
	info("   Original sequence length:  %zu\n", ref_seq.size());
	info("   Diploid sequence length:   %zu\n", diploid_seq.size());
	info("   Small-scale variation:\n");
	info("      Number of SNPs:                    %zu\n", snp_counter.count);
	info("      Number of small insertions:        %zu\n", small_insert_counter.count);
	info("      Number of small deletions:         %zu\n", small_delet_counter.count);
	info("      Total length of small insertions:  %zu\n", small_insert_counter.total_len);
	info("      Total length of small deletions:   %zu\n", small_delet_counter.total_len);
	info("   Structural variation:\n");
	info("      Number of large-scale insertions:        %zu\n", sv_insert_counter.count);
	info("      Number of large-scale deletions:         %zu\n", sv_delet_counter.count);
	info("      Number of large-scale inversions:        %zu\n", sv_invert_counter.count);
	info("      Total length of large-scale insertions:  %zu\n", sv_insert_counter.total_len);
	info("      Total length of large-scale deletions:   %zu\n", sv_delet_counter.total_len);
	info("      Total length of large-scale inversions:  %zu\n", sv_invert_counter.total_len);
	info("\n");
}

/* Write a sequence to a FASTA file. */
static void write_seq(vector<char> &seq, string &id)
{
	static const int LINE_LEN = 70;
	size_t i;

	out_file.putc('>');
	out_file.puts(id);
	out_file.putc('\n');
	for (i = 0; i + LINE_LEN < seq.size(); i += LINE_LEN) {
		out_file.write(&seq[i], LINE_LEN);
		out_file.putc('\n');
	}
	out_file.write(&seq[i], seq.size() - i);
	out_file.putc('\n');
}

/*
 * Begin simulating the diploid genome.
 */
static void do_pirs_diploid()
{
	info("\n");
	info("Beginning diploid simulation with the following parameters:\n");
	info("   SNP rate = %g\n", snp_rate);
	info("   Small insertion rate = %g\n", small_insertion_rate);
	info("   Small deletion rate = %g\n", small_deletion_rate);
	info("   Large insertion rate = %g\n", sv_insertion_rate);
	info("   Large deletion rate = %g\n", sv_deletion_rate);
	info("   Large inversion rate = %g\n", sv_inversion_rate);
	info("   SNP transition to transversion ratio = %g\n",
			snp_transition_to_transversion_ratio);

	info("The command line was:\n");
	info("    %s\n", get_command_line());
	info("The date is: %s\n", get_timestamp());
	info("The actual cumulative heterozygosity probability table that will be used is:\n");
	info("   [%#" PRIx64", %#" PRIx64", %#" PRIx64"\n",
	     het_type_tab[0], het_type_tab[1], het_type_tab[2]);
	info("    %#" PRIx64", %#" PRIx64", %#" PRIx64"\n",
	     het_type_tab[3], het_type_tab[4], het_type_tab[5]);
	info("    %#" PRIx64"]\n", het_type_tab[6]);
	info("\n");

	vector<char> diploid_seq;
	vector<char> ref_seq;
	string id;

	while (read_scaffold(in_file, id, ref_seq)) {
		info("Processing scaffold \"%s\" (length = %zu)\n", 
					id.c_str(), ref_seq.size());
		diploid_seq.clear();
		diploid_seq.reserve(ref_seq.size());
		make_diploid_seq(diploid_seq, ref_seq, id.c_str());
		info("Outputting diploid sequence \"%s\"\n", id.c_str());
		write_seq(diploid_seq, id);
	}
	info("Done generating diploid scaffolds\n");
}

/*
 * Prepares the cumulative probabality array for heterozygosity types.
 */
static void prepare_het_type_tab()
{
	small_insertion_rate = small_indel_rate / 2.0;
	small_deletion_rate  = small_indel_rate / 2.0;
	sv_insertion_rate    = sv_rate / 3.0;
	sv_deletion_rate     = sv_rate / 3.0;
	sv_inversion_rate    = sv_rate / 3.0;

	het_type_tab[1] = ~0ULL * snp_rate;
	het_type_tab[2] = ~0ULL * small_insertion_rate;
	het_type_tab[3] = ~0ULL * small_deletion_rate;
	het_type_tab[4] = ~0ULL * sv_insertion_rate;
	het_type_tab[5] = ~0ULL * sv_deletion_rate;
	het_type_tab[6] = ~0ULL * sv_inversion_rate;

	uint64_t total_het_rate_u64 = het_type_tab[1];
	for (size_t i = 2; i < ARRAY_LEN(het_type_tab); i++) {
		if (total_het_rate_u64 + het_type_tab[i] < total_het_rate_u64) {
			fatal_error("Cannot have the sum of probabilities of "
				    "heterozygosity types add up to more than 1.0!");
		}
		total_het_rate_u64 += het_type_tab[i];
	}
	het_type_tab[0] = ~0ULL - total_het_rate_u64;
	uint64_t cum_sum = 0;
	for (size_t i = 0; i < ARRAY_LEN(het_type_tab); i++) {
		cum_sum += het_type_tab[i];
		het_type_tab[i] = cum_sum;
	}
}

/*
 * This is the main procedure that is executed when the `pirs diploid' command
 * is given.
 */
void pirs_diploid(int argc, char *argv[])
{
	argc--;
	argv++;

	pirs_diploid_parse_command_line(argc, argv);

	info("Preparing heterozygosity probability table\n");
	prepare_het_type_tab();

	in_file.open(ref_filename);

	string out_filename;
	if (output_filename) {
		out_filename = output_filename;
	} else {
		out_filename = output_prefix;
		if (snp_rate > 0.0)
			out_filename += ".snp";

		if (small_indel_rate > 0.0 || sv_rate > 0.0)
			out_filename += ".indel";

		if (sv_rate > 0.0)
			out_filename += ".inversion";

		out_filename += ".fa";
	}
	out_file.open(out_filename);

	if (write_log_files) {
		if (snp_rate > 0.0) {
			snp_log_file.open(output_prefix + ".snp.lst");
			snp_log_file.write_header();
			snp_log_file.puts("# seq_id\torig_position\t"
					  "new_position\torig_base\t"
					  "new_base\n");
		}
		
		if ((small_indel_rate > 0.0 || sv_rate > 0.0)) {
			indel_log_file.open(output_prefix + ".indel.lst");
			indel_log_file.write_header();
			indel_log_file.puts("# seq_id\torig_position\t"
					    "new_position\ttype\t"
					    "indel_length\tindel_sequence\n");
		}

		if (sv_rate > 0.0) {
			inversion_log_file.open(output_prefix + ".inversion.lst");
			inversion_log_file.write_header();
			inversion_log_file.puts("# seq_id\torig_position\t"
						"new_position\t"
						"inversion_length\n");
		}
	}

	do_pirs_diploid();
}

// genodsp.c-- (GENOme Digital Signal Processing) read a list of genomic
//             intervals, with values, and perform base-granularity digital
//             signal processing operations on the waveforms.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include "utilities.h"

// program revision vitals (not the best way to do this!))

char* programName = "genodsp";

#define programVersionMajor    "0"
#define programVersionMinor    "0"
#define programVersionSubMinor "8"
#define programRevisionDate    "20220615"

// dsp operators

#define  globals_owner			// (make this the owner of the global variables)
#include "genodsp_interface.h"
#include "sum.h"
#include "clump.h"
#include "percentile.h"
#include "add.h"
#include "multiply.h"
#include "mask.h"
#include "logical.h"
#include "minmax.h"
#include "morphology.h"
#include "map.h"
#include "opio.h"
#include "variables.h"

//----------
//
// global data and types--
//
//----------

// command line options

dspop*		pipeline         = NULL;
dspop*		tailOp           = NULL;

int			valColumn        = 4-1;
int			noOutputValues   = false;
int			valPrecision     = 0;
int			collapseRuns     = true;
int			showUncovered    = uncovered_hide;  // (one of uncovered_show, etc.)
int			clipToLength     = false;
int			originOne        = false;
int			inhibitOutput    = false;

int			dbgInput         = false;
int			dbgPipe          = false;
int			dbgGlobals       = false;

// linked lists for scratch vectors

typedef struct svspec			// vector of values
	{
	struct svspec* next;		// next spec in a linked list
	int			inUse;			// true => someone is using this vector
	valtype*	vector;			// vector of values
	} svspec;

typedef struct svispec			// vector of integers
	{
	struct svispec* next;		// next spec in a linked list
	int			inUse;			// true => someone is using this vector
	s32*		vector;			// vector of values
	} svispec;

// linked lists for named global variables

typedef struct namedglobal
	{
	struct namedglobal* next;	// next variable in a linked list
	char*	name;				// variable's name (allocated within this block)
	valtype	v;					// variable's value
	} namedglobal;

// miscellany

#define min_of(a,b) ((a <= b)? a: b)

//----------
//
// prototypes--
//
//----------

int main (int argc, char** argv);

// private functions

static void  parse_options              (int _argc, char** _argv);
static int   process_operator_options   (int _argc, char** _argv);
static void  read_chromosome_lengths    (char* filename);
static int   add_chromosome_spec        (char* name,
                                         u32 chromStart, u32 chromLength);
static void  sort_chromosomes_by_length (void);
static void  init_scratch_vectors       (u32 scratchLength);
static void  free_scratch_vectors       (void);
static void  init_named_globals         (void);
static void  free_named_globals         (void);

// dsp operations table

dspinfo dspTable[] =
	{dspinforecord("sum"           , op_window_sum)     ,
	 dspinfoalias ("window_sum")                        ,
	 dspinforecord("slidingsum"    , op_sliding_sum)    ,
	 dspinfoalias ("sliding_sum")                       ,
	 dspinforecord("smooth"        , op_smooth)         ,
	 dspinforecord("cumulativesum" , op_cumulative_sum) ,
	 dspinfoalias ("cumulative")                        ,
	 dspinfoalias ("integrate")                         ,
	 dspinforecord("clump"         , op_clump)          ,
	 dspinforecord("anticlump"     , op_skimp)          ,
	 dspinfoalias ("anti_clump")                        ,
	 dspinfoalias ("skimp")                             ,
	 dspinforecord("percentile"    , op_percentile)     ,
	 dspinforecord("add"           , op_add)            ,
	 dspinforecord("subtract"      , op_subtract)       ,
	 dspinforecord("addconst"      , op_add_constant)   ,
	 dspinfoalias ("add_const")                         ,
	 dspinforecord("invert"        , op_invert)         ,
	 dspinforecord("multiply"      , op_multiply)       ,
	 dspinforecord("divide"        , op_divide)         ,
	 dspinforecord("abs"           , op_absolute_value) ,
	 dspinforecord("mask"          , op_mask)           ,
	 dspinforecord("masknot"       , op_mask_not)       ,
	 dspinfoalias ("mask_not")                          ,
	 dspinforecord("clip"          , op_clip)           ,
	 dspinforecord("erase"         , op_erase)          ,
	 dspinforecord("binarize"      , op_binarize)       ,
	 dspinforecord("or"            , op_or)             ,
	 dspinforecord("and"           , op_and)            ,
	 dspinforecord("maxover"       , op_max_in_interval),
	 dspinfoalias ("max_over")                          ,
	 dspinforecord("minover"       , op_min_in_interval),
	 dspinfoalias ("min_over")                          ,
	 dspinforecord("localmin"      , op_local_minima)   ,
	 dspinfoalias ("local_min")                         ,
	 dspinforecord("localmax"      , op_local_maxima)   ,
	 dspinfoalias ("local_max")                         ,
	 dspinforecord("bestmin",        op_best_local_min) ,
	 dspinfoalias ("best_min")                          ,
	 dspinfoalias ("bestlocalmin")                      ,
	 dspinfoalias ("best_local_min")                    ,
	 dspinforecord("bestmax",        op_best_local_max) ,
	 dspinfoalias ("best_max")                          ,
	 dspinfoalias ("bestlocalmax")                      ,
	 dspinfoalias ("best_local_max")                    ,
	 dspinforecord("minwith",        op_min_with)       ,
	 dspinfoalias ("min_with")                          ,
	 dspinforecord("maxwith",        op_max_with)       ,
	 dspinfoalias ("max_with")                          ,
	 dspinforecord("close"         , op_close)          ,
	 dspinforecord("open"          , op_open)           ,
	 dspinforecord("dilate"        , op_dilate)         ,
	 dspinforecord("erode"         , op_erode)          ,
	 dspinforecord("map"           , op_map)            ,
	 dspinforecord("input"         , op_input)          ,
	 dspinforecord("output"        , op_output)         ,
	 dspinforecord("variables"     , op_show_variables) };

#define dspTableLen (sizeof(dspTable)/sizeof(dspinfo))

//----------
//
// option parsing--
//
//----------

#define specialPipeChar '='

static opfunc_usage chastiseUsage     = NULL;
static char*        chastiseUsageName = NULL;

static void usage            (char* message);
static void usage_operations (void);


void chastise (const char* format, ...)
	{
	va_list	args;

	va_start (args, format);
	if (format != NULL)
		vfprintf (stderr, format, args);
	va_end (args);

	if (chastiseUsage != NULL)
		{
		(*chastiseUsage) (chastiseUsageName, stderr, "  ");
		exit (EXIT_FAILURE);
		}

	usage (NULL);
	}


static void usage (char* message)
	{
	if (message != NULL) fprintf (stderr, "%s\n", message);

	fprintf (stderr, "usage: [cat <file>] | %s --chromosomes=<filename> [options] [operations]\n", programName);
	fprintf (stderr, "\n");
	//                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (stderr, "  --chromosomes=<filename>  (required) read chromosome names and lengths from\n");
	fprintf (stderr, "                            a file\n");
	fprintf (stderr, "  --value=<col>             input intervals contain a value in the specified\n");
	fprintf (stderr, "                            column;  by default we assume this is in column 4\n");
	fprintf (stderr, "  --novalue                 input intervals have no value (value given is 1)\n");
	fprintf (stderr, "  --nooutputvalue           don't write value with output intervals\n");
	fprintf (stderr, "  --precision=<number>      number of digits to round output values to\n");
	fprintf (stderr, "                            (by default, output is rounded to integers)\n");
	fprintf (stderr, "  --nocollapse              in output, don't collapse runs of identical values\n");
	fprintf (stderr, "                            to intervals\n");
	fprintf (stderr, "  --uncovered:hide          don't output intervals that have no coverage\n");
	fprintf (stderr, "                            (this is the default)\n");
	fprintf (stderr, "  --uncovered:show          in output, include intervals that have no coverage\n");
	fprintf (stderr, "  --uncovered:NA            in output, mark uncovered intervals as NA\n");
	fprintf (stderr, "  --cliptochromosome        clip interals to chromosome length\n");
	fprintf (stderr, "                            (default is to report such intervals as errors)\n");
	fprintf (stderr, "  --origin=one              input/output intervals are origin-one, closed\n");
	fprintf (stderr, "  --origin=zero             input/output intervals are origin-zero, half-open\n");
	fprintf (stderr, "                            (this is the default)\n");
	fprintf (stderr, "  --nooutput                don't output the resulting intervals/values\n");
	fprintf (stderr, "                            (by default these are written to stdout)\n");
	fprintf (stderr, "  --window=<length>         (W=) size of window\n");
	fprintf (stderr, "                            (for operators that have a window size)\n");
	fprintf (stderr, "  --help[=<operator>]       get detail about a particular operator\n");
	fprintf (stderr, "  ?                         list available operators with brief descriptions\n");
	fprintf (stderr, "  ?<operator>               same as --help=<operator>\n");
	fprintf (stderr, "  --report=comments         copy comments from the input to stderr. Comments\n");
	fprintf (stderr, "                            are lines beginning with a \"#\". This can be\n");
	fprintf (stderr, "                            helpful in tracking progress during a long run.\n");
	fprintf (stderr, "  --progress=input:<n>      report processing of every nth input line\n");
	fprintf (stderr, "  --progress=operations     report each operation as it begins\n");
	fprintf (stderr, "  --version                 report the program version and quit\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Note that if input intervals overlap, their values are summed.\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Input is usually piped in on stdin. However, if the first operator is \"input\"\n");
	fprintf (stderr, "stdin is ignored.\n");
	fprintf (stderr, "\n");

	fprintf (stderr, "For a list of available operations, do \"genodsp ?\".\n");
	fprintf (stderr, "For more detailed descriptions of the operations, do \"genodsp --help\".\n");

	exit (EXIT_FAILURE);
	}


static void usage_operations (void)
	{
	dspinfo*	opInfo;
	u32			dspIx;

	fprintf (stderr, "Operations (general form is %c <operator> [arguments]):\n",
	                 specialPipeChar);

	for (dspIx=0 ; dspIx<dspTableLen ; dspIx++)
		{
		opInfo = &dspTable[dspIx];
		if (opInfo->funcShort == NULL) continue;
		(*opInfo->funcShort) (opInfo->name, 12, stderr, "  ");
		}

	exit (EXIT_FAILURE);
	}


static void parse_options (int _argc, char** _argv)
	{
	int			argc;
	char**		argv;
	char*		arg, *argVal, *argScan, *argScan2;
	char*		chromsFilename = NULL;
	u32			chromStart = 0;
	u32			chromLength;
	dspinfo*	opInfo, *realOpInfo;
	u32			dspIx;
	int			argsConsumed;
	int			tempInt;

	// skip program name

	//programName = _argv[0];
	argv = _argv+1;  argc = _argc - 1;

	//////////
	// scan arguments
	//////////

	if (argc == 0)
		chastise (NULL);

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// operation

		if (arg[0] == specialPipeChar)
			{
			if (dbgPipe)
				fprintf (stderr, "pipe char in arg[%d]\n", (int) (argv-_argv));

			if ((argc == 1) && (arg[0] == specialPipeChar) && (arg[1] == 0))
				chastise ("%c at end of command line, with no operation\n",
				          specialPipeChar);
			argsConsumed = process_operator_options (argc, argv);
			argv += argsConsumed - 1;
			argc -= argsConsumed - 1;
			goto next_arg;
			}

		// --chromosomes=<filename>

		if ((strcmp_prefix (arg, "--chromosomes=") == 0)
		 || (strcmp_prefix (arg, "--chroms=")      == 0))
			{
			chromsFilename = argVal;
			goto next_arg;
			}

		// --value=<col>

		if ((strcmp (arg, "--novalue") == 0)
		 || (strcmp (arg, "--novalues") == 0)
		 || (strcmp (arg, "--value=none") == 0))
			{
			valColumn = -1;
			set_named_global ("valColumn", (valtype) valColumn);
			goto next_arg;
			}

		if (strcmp_prefix (arg, "--value=") == 0)
			{
			valColumn = string_to_int (argVal) - 1;
			if (valColumn == -1)
				chastise ("value column can't be 0 (\"%s\")\n", arg);
			if (valColumn < 0)
				chastise ("value column can't be negative (\"%s\")\n", arg);
			if (valColumn < 3)
				chastise ("value column can't be 1, 2 or 3 (\"%s\")\n", arg);
			set_named_global ("valColumn", (valtype) valColumn);
			goto next_arg;
			}

		if ((strcmp (arg, "--nooutputvalue")  == 0)
		 || (strcmp (arg, "--nooutputvalues") == 0))
			{
			noOutputValues = true;
			set_named_global ("noOutputValues", (valtype) noOutputValues);
			goto next_arg;
			}

		// --precision=<col>

		if (strcmp_prefix (arg, "--precision=") == 0)
			{
			valPrecision = string_to_int (argVal);
			if (valPrecision < 0)
				chastise ("precision can't be negative (\"%s\")\n", arg);
			set_named_global ("valPrecision", (valtype) valPrecision);
			goto next_arg;
			}

		// --nocollapse

		if (strcmp (arg, "--nocollapse") == 0)
			{
			collapseRuns = false;
			set_named_global ("collapseRuns", (valtype) collapseRuns);
			goto next_arg;
			}

		// --uncovered:hide

		if ((strcmp (arg, "--uncovered:hide") == 0)
		 || (strcmp (arg, "--hide:uncovered") == 0))
			{
			showUncovered = uncovered_hide;
			set_named_global ("showUncovered", (valtype) showUncovered);
			goto next_arg;
			}

		// --uncovered:show

		if ((strcmp (arg, "--uncovered:show") == 0)
		 || (strcmp (arg, "--show:uncovered") == 0))
			{
			showUncovered = uncovered_show;
			set_named_global ("showUncovered", (valtype) showUncovered);
			goto next_arg;
			}

		// --uncovered:NA

		if ((strcmp (arg, "--uncovered:NA") == 0)
		 || (strcmp (arg, "--uncovered:mark") == 0)
		 || (strcmp (arg, "--mark:uncovered") == 0)
		 || (strcmp (arg, "--markgaps")       == 0))
			{
			showUncovered = uncovered_NA;
			set_named_global ("showUncovered", (valtype) showUncovered);
			goto next_arg;
			}

		// --cliptochromosome

		if ((strcmp (arg, "--cliptochromosome") == 0)
		 || (strcmp (arg, "--cliptochrom")      == 0)
		 || (strcmp (arg, "--cliptolength")     == 0)
		 || (strcmp (arg, "--clip")             == 0))
			{
			clipToLength = true;
			goto next_arg;
			}

		// --origin=one, --origin=zero

		if ((strcmp (arg, "--origin=one") == 0)
		 || (strcmp (arg, "--origin=1")   == 0))
			{
			originOne = true;
			set_named_global ("originOne", (valtype) originOne);
			goto next_arg;
			}

		if ((strcmp (arg, "--origin=zero") == 0)
		 || (strcmp (arg, "--origin=0")    == 0))
			{
			originOne = false;
			set_named_global ("originOne", (valtype) originOne);
			goto next_arg;
			}

		// --nooutput

		if (strcmp (arg, "--nooutput") == 0)
			{ inhibitOutput = true;  goto next_arg; }

		// --window=<length> or W=<length>

		if ((strcmp_prefix (arg, "--window=") == 0)
		 || (strcmp_prefix (arg, "W=")        == 0)
		 || (strcmp_prefix (arg, "--W=")      == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt == 0)
				chastise ("window size can't be zero (\"%s\")\n", arg);
			if (tempInt < 0)
				chastise ("window size can't be negative (\"%s\")\n", arg);
			set_named_global ("windowSize", (valtype) tempInt);
			goto next_arg;
			}

		// --help and --help=<operator> (and variations starting with '?')

		if (strcmp (arg, "?") == 0)
			{
			usage_operations ();
			exit (EXIT_SUCCESS);
			}

		if (strcmp_prefix (arg, "?") == 0)
			{ argVal = arg+1;  goto help_for_one; }

		if ((strcmp_prefix (arg, "--help=") == 0)
		 || (strcmp_prefix (arg, "?=")      == 0))
			{
			if (strcmp (argVal, "*") == 0)
				goto help_for_all_operators;
		help_for_one:
			opInfo = realOpInfo = NULL;
			for (dspIx=0 ; dspIx<dspTableLen ; dspIx++)
				{
				if (dspTable[dspIx].funcShort != NULL) realOpInfo = &dspTable[dspIx];
				if (strcmp (argVal,dspTable[dspIx].name) != 0) continue;
				opInfo = realOpInfo;
				break;
				}
			if (opInfo == NULL) goto help_for_unknown;
			fprintf (stderr, "=== %s ===\n", opInfo->name);
			(*opInfo->funcUsage) (opInfo->name, stderr, "  ");
			exit (EXIT_SUCCESS);
			}

		if (strcmp (arg, "--help") == 0)
			{
		help_for_all_operators:
			for (dspIx=0 ; dspIx<dspTableLen ; dspIx++)
				{
				opInfo = &dspTable[dspIx];
				if (opInfo->funcShort == NULL) continue;
				fprintf (stderr, "=== %s ===\n", opInfo->name);
				(*opInfo->funcUsage) (opInfo->name, stderr, "  ");
				}
			exit (EXIT_SUCCESS);
			}

		// --report=comments

		if ((strcmp (arg, "--report=comments") == 0)
		 || (strcmp (arg, "--report:comments") == 0))
			{ reportComments = true;  goto next_arg; }

		// --report=progress=input:<n>

		if ((strcmp_prefix (arg, "--progress=input:") == 0)
		 || (strcmp_prefix (arg, "--progress:input=") == 0)
		 || (strcmp_prefix (arg, "--progress:input:") == 0))
			{
			if (strcmp_prefix (argVal, "input:") == 0)
				argVal = strchr(arg,':') + 1;
			reportInputProgress = string_to_unitized_int (argVal, /*thousands*/ true);
			goto next_arg;
			}

		// --progress

		if ((strcmp (arg, "--progress=operations") == 0)
		 || (strcmp (arg, "--progress:operations") == 0)
		 || (strcmp (arg, "--debug=operations")    == 0))  // backward compatibility
			{ trackOperations = true;  goto next_arg; }

		// --version

		if (strcmp (arg, "--version") == 0)
			{
			fprintf (stderr, "%s (version %s.%s.%s released %s)\n",
			                 programName,
			                 programVersionMajor, programVersionMinor, programVersionSubMinor, programRevisionDate);
			exit (EXIT_SUCCESS);
			}

		// --debug arguments

		if (strcmp (arg, "--debug=input") == 0)
			{ dbgInput = true;  goto next_arg; }

		if (strcmp (arg, "--debug=pipe") == 0)
			{ dbgPipe = true;  goto next_arg; }

		if (strcmp (arg, "--debug=globals") == 0)
			{ dbgGlobals = true;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("Can't understand \"%s\"\n", arg);

		// (undocumented) <chromosome>:<length> or <chromosome>:<start>:<end>
		// .. in this case <start> and <end> are origin-zero half-open,
		// .. regardless of any user setting

		argScan = strchr(arg,':');
		if (argScan == NULL) goto no_length;
		argScan2 = strchr(argScan+1,':');
		if (argScan2 == NULL)
			{
			chromStart   = 0;
			*(argScan++) = 0;
			chromLength  = string_to_u32 (argScan);
			}
		else
			{
			*(argScan++)  = 0;
			*(argScan2++) = 0;
			chromStart    = string_to_u32 (argScan);
			chromLength   = string_to_u32 (argScan2) - chromStart;
			}

		if (!add_chromosome_spec (arg, chromStart, chromLength))
			chastise ("can't specify %s more than once\n", arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	//////////
	// read chromosome lengths
	//////////

	if (chromsFilename != NULL)
		read_chromosome_lengths (chromsFilename);

	//////////
	// sanity checks
	//////////

	// make sure we got at least one chromosome-of-interest

	if (chromsOfInterest == NULL)
		chastise ("gotta give me some chromosome names\n");

	return;

	//////////
	// failure exits
	//////////

help_for_unknown:
	fprintf (stderr, "\"%s\" is not a known operation\n",
	                 argVal);
	exit(EXIT_FAILURE);

no_length:
	fprintf (stderr, "\"%s\" contains no chromosome length\n"
	                 "(expected \"chromosome:length\" or \"chromosome:start:end\")\n",
			 arg);
	exit (EXIT_FAILURE);
	}


// process_operator_options--

static int process_operator_options (int _argc, char** _argv)
	{
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg;
	int			argsConsumed, dspArgC, argIx;
	char*		dspName;
	dspinfo*	opInfo, *realOpInfo;
	u32			dspIx;
	dspop*		op;

	argv = _argv;  argc = _argc;

	// process the operator name

	arg          = argv[0];
	argsConsumed = 0;
	dspArgC      = 0;

	if ((arg[0] == specialPipeChar) && (arg[1] == 0))
		{
		argv++;
		argc--;			// (this will never go to zero)
		argsConsumed++;
		dspName = argv[0];
		}
	else
		dspName = skip_whitespace (arg+1);

	if (dbgPipe)
		fprintf (stderr, "  dspName=\"%s\"\n", dspName);

	opInfo = realOpInfo = NULL;
	for (dspIx=0 ; dspIx<dspTableLen ; dspIx++)
		{
		if (dspTable[dspIx].funcShort != NULL) realOpInfo = &dspTable[dspIx];
		if (strcmp (dspName,dspTable[dspIx].name) != 0) continue;
		opInfo = realOpInfo;
		break;
		}
	if (opInfo == NULL)
		chastise ("\"%s\" is not a known operation\n", dspName);

	argv++; argc--; argsConsumed++; // skip operator name

	// scan ahead to locate the terminating arg (next arg that starts with the
	// internal pipe character, or end of arg list)

	for (argIx=0 ; argIx<argc ; argIx++)
		{
		arg = argv[argIx];
		if (arg[0] == specialPipeChar) break;
		argsConsumed++;
		}

	dspArgC = argIx;

	if (dbgPipe)
		{
		fprintf (stderr, "  args:");
		for (argIx=0 ; argIx<argc ; argIx++)
			fprintf (stderr, " %s", argv[argIx]);
		fprintf (stderr, "\n");
		}

	// tell operator to parse its arguments

	chastiseUsage     = opInfo->funcUsage;
	chastiseUsageName = opInfo->name;
	op = (*opInfo->funcParse) (opInfo->name, dspArgC, argv);
	chastiseUsage     = NULL;
	chastiseUsageName = NULL;

	op->name      = copy_string (opInfo->name);
	op->funcApply = opInfo->funcApply;
	op->funcFree  = opInfo->funcFree;
	// op->atRandom must be set by the parse function

	// add operation item to the pipeline

	if (tailOp == NULL) pipeline     = op;
				   else tailOp->next = op;
	op->next = NULL;
	tailOp   = op;

	if (dbgPipe)
		fprintf (stderr, "  argsConsumed=%d\n", argsConsumed);

	return argsConsumed;
	}


// read_chromosome_lengths--

static void read_chromosome_lengths (char* filename)
	{
	FILE*	f;
	char	lineBuffer[1001];
	u32		lineNumber = 0;
	int		missingEol = false;
	int		lineLen;
	char*	scan, *mark, *field;
	char*	chrom;
	u32		chromLength;

	f = fopen (filename, "rt");
	if (f == NULL) goto cant_open_file;

	lineNumber = 0;
	while (true)
		{
		if (fgets (lineBuffer, sizeof(lineBuffer), f) == NULL)
			break;

		lineNumber++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way)

		if (missingEol) goto missing_eol;

		lineLen = strlen(lineBuffer);
		if (lineLen != 0)
			missingEol = (lineBuffer[lineLen-1] != '\n');

		// parse the line

		scan = skip_whitespace(lineBuffer);
		if (*scan == 0)   continue;  // empty line
		if (*scan == '#') continue;  // comment line

		chrom = scan = lineBuffer;
		if (*scan == 0) goto no_chrom;
		mark = skip_darkspace(scan);
		scan = skip_whitespace(mark);
		if (*mark != 0) *mark = 0;

		if (*scan == 0) goto no_length;
		field = scan;
		mark = skip_darkspace(scan);
		scan = skip_whitespace(mark);
		if (*mark != 0) *mark = 0;
		chromLength = string_to_u32 (field);

		if (!add_chromosome_spec (chrom, 0, chromLength))
			goto same_chromosome;
		}

	// success

	fclose (f);
	return;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "can't open \"%s\" for reading\n", filename);
	exit (EXIT_FAILURE);

missing_eol:
	fprintf (stderr, "problem at line %u, line is longer than internal buffer\n",
			 lineNumber-1);
	exit (EXIT_FAILURE);

no_chrom:
	fprintf (stderr, "problem at line %u, line contains no chromosome or begins with whitespace\n",
			 lineNumber);
	exit (EXIT_FAILURE);

no_length:
	fprintf (stderr, "problem at line %u, line contains no chromosome length\n",
			 lineNumber);
	exit (EXIT_FAILURE);

same_chromosome:
	fprintf (stderr, "problem at line %u, chromosome \"%s\" appears more than once\n",
			 lineNumber, chrom);
	exit (EXIT_FAILURE);
	}

//----------
//
// main program--
//
//----------

int main
   (int			argc,
	char**		argv)
	{
	char*		chrom;
	u32			vLen;
	valtype*	v;
	spec*		chromSpec;
	dspop*		firstOp, *stopOp, *op, *nextOp;
	u32			maxLength;
	u32			ix, chromIx;
	opfunc_free	funcFree;

	init_named_globals ();
	set_named_global ("valColumn",     (valtype) valColumn);
	set_named_global ("valPrecision",  (valtype) valPrecision);
	set_named_global ("collapseRuns",  (valtype) collapseRuns);
	set_named_global ("showUncovered", (valtype) showUncovered);
	set_named_global ("originOne",     (valtype) originOne);

	parse_options (argc, argv);

	//////////
	// allocate vectors
	//////////

	sort_chromosomes_by_length ();

	// determine the length of any scratch vectors;  every scratch vector will
	// be as long as the longest chromosome

	maxLength = 0;
	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		if (chromSpec->length == 0) goto no_length_specified;
		if (chromSpec->length > maxLength) maxLength = chromSpec->length;
		}

	init_scratch_vectors (maxLength);

	// allocate chromosome value vectors

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];

		if (trackOperations)
			tracking_report ("allocate(%s / %s bytes)\n",
			                  chromSpec->chrom, ucommatize(chromSpec->length));

		chromSpec->valVector = (valtype*) calloc (chromSpec->length, sizeof(valtype));
		if (chromSpec->valVector == NULL) goto cant_allocate_val;
		v = chromSpec->valVector;
		for (ix=0 ; ix<chromSpec->length ; ix++)
			v[ix] = 0.0;
		}

	if (trackOperations)
		tracking_report ("allocate(--done--)\n");

	//////////
	// process intervals
	//////////

	// read intervals and values;  note that if the first operator is an
	// input operator, we avoid this step (since the first operator will
	// overwrite whatever intervals we might have read)

	op = pipeline;
	if ((op == NULL) || (strcmp (op->name, "input") != 0))
		read_intervals (stdin, valColumn, originOne, ri_overlapSum, /*clear*/ false, 0.0);

	// perform operations;  when possible, we perform a series of operations on
	// one chromosome before moving onto the next;  it is expected that this
	// will have some improvement in cache performance, at least for chromosomes
	// smaller than the cache

	firstOp = pipeline;
	while (firstOp != NULL)
		{
		for (stopOp=firstOp ; stopOp!=NULL ; stopOp=stopOp->next)
			{ if (stopOp->atRandom) break; }

		if (stopOp != firstOp)
			{
			// run several operations serially on each chromosome
			for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
				{
				chromSpec = chromsSorted[chromIx];
				for (op=firstOp ; op!=stopOp ; op=op->next)
					{
					chrom = chromSpec->chrom;
					vLen  = chromSpec->length;
					v     = chromSpec->valVector;
					if (trackOperations)
						fprintf (stderr, "%s(%s)\n", op->name, chrom);
					(*op->funcApply) (op, chrom, vLen, v);
					}
				}
			}

		if (stopOp == NULL)
			firstOp = NULL;
		else
			{
			// run one operation on all chromosomes 'simultaneously'
			op = stopOp;

			if (trackOperations)
				tracking_report ("%s(*)\n", op->name);
			(*op->funcApply) (op, "*", maxLength, NULL);
			firstOp = stopOp->next;
			}
		}

	// report intervals

	if (!inhibitOutput)
		report_intervals (stdout, valPrecision, noOutputValues,
		                  collapseRuns, showUncovered, originOne);

	//////////
	// success
	//////////

	// deallocate

	free_scratch_vectors ();
	free_named_globals   ();

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		if (chromSpec->chrom     != NULL) free (chromSpec->chrom);
		if (chromSpec->valVector != NULL) free (chromSpec->valVector);
		free (chromSpec);
		}
	chromsOfInterest = NULL;

	free (chromsSorted);
	chromsSorted = NULL;

	for (op=pipeline ; op!=NULL ; op=nextOp)
		{
		nextOp = op->next;
		if (op->name != NULL) { free (op->name);  op->name = NULL; }
		funcFree = op->funcFree;
		funcFree (op);
		}

	// proclaim success

	return EXIT_SUCCESS;

	//////////
	// failure exits
	//////////

no_length_specified:
	fprintf (stderr, "no length was specified for %s\n",
	                 chromSpec->chrom);
	return EXIT_FAILURE;

cant_allocate_val:
	fprintf (stderr, "failed to allocate %d-base vector for %s, %d bytes per base\n",
	                 chromSpec->length, chromSpec->chrom, (int) sizeof(valtype));
	return EXIT_FAILURE;
	}

//----------
//
// add_chromosome_spec--
//	Add a chromosome to our list of chromosome specifications.
//
//----------
//
// Arguments:
//	char*	name:			The name of the chromsome (e.g. "chr1").
//	u32		chromStart:		The first base of interest on the chromosome
//							.. (origin zero);  usually this is zero.
//	u32		chromLength:	The number of bases to allocate for chromosome;
//							.. usually this is the chromosome's length.  Note
//							.. that if this is zero, we don't bother adding
//							.. a record for the chromosome.
//
// Returns:
//	true if successful;  false if the chromosome is already in the list;  all
//	other failures result in program termination.
//
//----------

static int add_chromosome_spec
   (char*	name,
	u32		chromStart,
	u32		chromLength)
	{
	spec*	scanSpec, *tailSpec, *newSpec;

	if (chromLength == 0) return true;

	// check whether this chromosome name is already in use, and find the tail
	// of the list if the name is not in use

	tailSpec = NULL;
	for (scanSpec=chromsOfInterest ; scanSpec!=NULL ; scanSpec=scanSpec->next)
		{
		tailSpec = scanSpec;
		if (strcmp (name, scanSpec->chrom) == 0) return false;
		}

	// create a new spec for this chromosome and add it to the list

	newSpec = (spec*) malloc (sizeof(spec));
	if (newSpec == NULL) goto cant_allocate_spec;
	if (tailSpec == NULL) chromsOfInterest = newSpec;
					 else tailSpec->next   = newSpec;
	newSpec->next      = NULL;
	newSpec->chrom     = copy_string (name);
	newSpec->start     = chromStart;
	newSpec->length    = chromLength;
	newSpec->valVector = NULL;

	return true;

	//////////
	// failure exits
	//////////

cant_allocate_spec:
	fprintf (stderr, "failed to allocate spec for \"%s\", %d bytes\n",
	                 name, (int) sizeof(spec));
	exit(EXIT_FAILURE);
	return false; // (never reaches here)
	}


//----------
//
// find_chromosome_spec--
//	Locate the spec for a specific chromosome name.
//
//----------
//
// Arguments:
//	char*	chrom:	name of the chromosome to look for.
//
// Returns:
//	a pointer to the spec for the chromosome;  NULL if the chromosome is not
//	in our list.
//
//----------

spec* find_chromosome_spec
   (char*	chrom)
	{
	spec*	scanSpec;

	for (scanSpec=chromsOfInterest ; scanSpec!=NULL ; scanSpec=scanSpec->next)
		{ if (strcmp (chrom, scanSpec->chrom) == 0) return scanSpec; }

	return NULL;
	}


//----------
//
// sort_chromosomes_by_length--
//	Make an array of our list of chromosome specifications, sorted by
//	decreasing length.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	(nothing)
//
//----------

int spec_length_descending (const void* _v1, const void* _v2);
int spec_length_descending (const void* _v1, const void* _v2)
	{
	const u32 v1 = (*(const spec**) _v1)->length;
	const u32 v2 = (*(const spec**) _v2)->length;
	return (v1 < v2) - (v1 > v2);
	}

// sort_chromosomes_by_length--

static void sort_chromosomes_by_length
   (void)
	{
	spec*	chromSpec;
	int		numChroms, chromIx;
	u32		numBytes;

	// count the chromosomes

	numChroms = 0;
	for (chromSpec=chromsOfInterest ; chromSpec!=NULL ; chromSpec=chromSpec->next)
		numChroms++;

	// allocate an array of pointers to the chromosome specs;  we allocate one
	// extra entry at the end, to use as a terminator.

	numBytes = (numChroms+1) * sizeof(spec*);
	chromsSorted = (spec**) malloc (numBytes);
	if (chromsSorted == NULL) goto cant_allocate;

	// copy the pointers into the array

	chromIx = 0;
	for (chromSpec=chromsOfInterest ; chromSpec!=NULL ; chromSpec=chromSpec->next)
		chromsSorted[chromIx++] = chromSpec;

	chromsSorted[chromIx] = NULL; // (terminator)

	// sort 'em

	qsort (chromsSorted, numChroms, sizeof(spec*), spec_length_descending);

	return;

	//////////
	// failure exits
	//////////

cant_allocate:
	fprintf (stderr, "failed to allocate sorted chromosome list, %d bytes\n",
	                 numBytes);
	exit(EXIT_FAILURE);
	}


//----------
//
// read_intervals--
//	Read all intervals from a file.
//
//----------
//
// Arguments:
//	FILE*	f:			File to read from.
//	int		valCol:		The column that contains interval value;  -1
//						.. indicates no such column.
//	int		originOne:	true  => interpret intervals as origin-one, closed
//						false => interpret intervals as origin-zero, half-open
//	int		overlapOp:	The operation to perform when a position is covered
//						.. by more than one interval;  one of ri_overlapSum",
//						.. etc.  When in doubt, sum is usually the appropriate
//						.. choice.
//	int		clear:		true  => zero all vectors before reading
//	valtype	missingVal:	value to clear to (only valid if clear is true).
//
// Returns:
//	nothing;  input failures result in program termination.
//
// Side effects:
//	The chromsOfInterest vectors are written, representing the values of the
//	incoming intervals.
//
//----------

void read_intervals
   (FILE*		f,
	int			valCol,
	int			originOne,
	int			overlapOp,
	int			clear,
	valtype		missingVal)
	{
	char		lineBuffer[1001];
	char		prevChrom[1001];
	valtype*	v = NULL;
	char*		chrom;
	spec*		chromSpec;
	u32			start, end, o, adjStart, adjEnd;
	valtype		val;
	u32			ix, chromIx;
	int			ok;

	// $$$ reset read_interval's lineNumber

	if (trackOperations)
		{
		for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
			{
			chromSpec = chromsSorted[chromIx];
			chromSpec->flag = false;
			}
		}

	// clear all chromosomes

	if (clear)
		{
		for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
			{
			chromSpec = chromsSorted[chromIx];

			chrom = chromSpec->chrom;
			v     = chromSpec->valVector;
			start = 0;
			end   = chromSpec->length;

			for (ix=start ; ix<end ; ix++)
				v[ix] = missingVal;
			}

		}

	// read intervals and values

	if (originOne) o = 1;
	          else o = 0;

	prevChrom[0] = 0;
	chromSpec    = NULL;

	v = NULL;
	while (true)
		{
		ok = read_interval (f, lineBuffer, sizeof(lineBuffer), valCol,
	                        &chrom, &start, &end, &val);
		if (!ok) break;

		//fprintf (stderr, "%s %u %u %f\n", chrom, start, end, val);

		if (strcmp (chrom, prevChrom) != 0)
			{
			v = NULL;
			chromSpec = find_chromosome_spec (chrom);
			if (chromSpec != NULL) v = chromSpec->valVector;
			strncpy (prevChrom, chrom, sizeof(prevChrom)-1);
			}

		if (chromSpec == NULL) continue;

		if ((trackOperations) && (!chromSpec->flag))
			{
			tracking_report ("input(%s)\n", chrom);
			chromSpec->flag = true;
			}

		start -= o;
		adjStart = start;
		adjEnd   = end;

		if (clipToLength)
			{
			// if the use has told us to clip intervals to the chromosome
			// length, do so

			if (start > chromSpec->start + chromSpec->length)
				adjStart = start = chromSpec->start + chromSpec->length;

			if (end > chromSpec->start + chromSpec->length)
				adjEnd = end = chromSpec->start + chromSpec->length;
			}

		if (chromSpec->start == 0)
			{
			// if only length has been specified, we *reject* intervals beyond
			// the end

			if (end > chromSpec->length) goto chrom_too_short;
			}
		else
			{
			// if start and end have been specified, we *ignore* intervals, or
			// portions of intervals, beyond the end

			if (end <= chromSpec->start) continue;

			adjEnd = end - chromSpec->start;
			if (start <= chromSpec->start) adjStart = 0;
			                          else adjStart = start - chromSpec->start;
			if (adjStart >= chromSpec->length) continue;
			if (adjEnd   >= chromSpec->length) adjEnd = chromSpec->length;
			}

		// "write" the value into the vector, across the interval

		if (overlapOp == ri_overlapMin)
			{
			for (ix=adjStart ; ix<adjEnd ; ix++)
				{
				if      ((clear) && (v[ix] == missingVal)) v[ix] = val;
				else if (val < v[ix])                      v[ix] = val;
				}
			}
		else if (overlapOp == ri_overlapMax)
			{
			for (ix=adjStart ; ix<adjEnd ; ix++)
				{
				if      ((clear) && (v[ix] == missingVal)) v[ix] = val;
				else if (val > v[ix])                      v[ix] = val;
				}
			}
		else // if (overlapOp == ri_overlapSum)
			{
			for (ix=adjStart ; ix<adjEnd ; ix++)
				{
				if ((clear) && (v[ix] == missingVal)) v[ix] =  val;
												 else v[ix] += val;
				}
			}
		}

	if (trackOperations)
		tracking_report ("input(--done--)\n");

	//////////
	// success
	//////////

	return;

	//////////
	// failure exits
	//////////

chrom_too_short:
	fprintf (stderr, "%s %d %d is beyond the end of the chromosome (L=%d)\n",
	                 chrom, start, end, chromSpec->length);
	exit (EXIT_FAILURE);
	}

//----------
//
// read_interval--
//	Read the next interval from a file.
//
//----------
//
// Arguments:
//	FILE*		f:			File to read from.
//	char*		buffer:		Buffer to read the line into.  Note that the caller
//							.. should not expect anything about the contents of
//							.. this buffer upon return.
//	int			bufferLen:	Number of bytes allocated for the buffer.
//	int			valCol:		The column that contains interval value;  -1
//							.. indicates no such column.
//	char**		chrom:		Place to return a pointer to the chromosome.  The
//							.. returned value will point into the line buffer,
//							.. and to a zero-terminated string.
//	u32*		start:		Place to return the start.
//	u32*		end:		Place to return the end.
//	valtype*	val:		Place to return the value.
//
// Returns:
//	true if we were successful;  false if there are no more lines in the file.
//	Failures result in program termination.
//
//----------

// $$$ need a way for the caller to reset lineNumber, perhaps just add a routine
//     .. reset_read_interval;  note that this is called from several op routines
//     .. in addition to read_intervals()

// $$$ also add an option to report "progress" to the user

int read_interval
   (FILE*		f,
	char*		buffer,
	int			bufferLen,
	int			valCol,
	char**		_chrom,
	u32*		_start,
	u32*		_end,
	valtype*	_val)
	{
	static u32	lineNumber = 0;
	static int	missingEol = false;
	int			reportProgressNow;
	int			lineLen;
	char*		scan, *mark, *field;
	int			col;
	char*		chrom;
	u32			start, end;
	valtype		val;

	// read the next line

try_again:

	if (fgets (buffer, bufferLen, f) == NULL)
		return false;

	lineNumber++;

	// check for lines getting split by fgets (the final line in the file might
	// not have a newline, but no internal lines can be that way)

	if (missingEol) goto missing_eol;

	lineLen = strlen(buffer);
	if (lineLen != 0)
		missingEol = (buffer[lineLen-1] != '\n');

	if (dbgInput)
		fprintf (stderr, "input = \"%s\"\n", buffer);

	// ignore track lines (so that the input may be a bedgraph file)

	if (strcmp_prefix (buffer, "track ") == 0)
		goto try_again;

	// parse the line

	reportProgressNow = ((reportInputProgress != 0)
	                  && ((lineNumber == 1) || (lineNumber % reportInputProgress == 0)));

	scan = skip_whitespace(buffer);
	if (*scan == 0)                    // empty line
		{
		if (reportProgressNow)
			fprintf (stderr, "progress: input line %u\n", lineNumber);
		goto try_again;
		}
	if (*scan == '#')                  // comment line
		{
		if (reportComments)
			fprintf (stderr, "input line %u: %s", lineNumber,scan);
		else if (reportProgressNow)
			fprintf (stderr, "progress: input line %u\n", lineNumber);
		goto try_again;
		}

	if (reportProgressNow)
		fprintf (stderr, "progress: input line %u\n", lineNumber);

	chrom = scan = buffer;
	if (*scan == ' ') goto no_chrom;
	mark = skip_darkspace(scan);
	scan = skip_whitespace(mark);
	if (*mark != 0) *mark = 0;

	if (*scan == 0) goto no_start;
	field = scan;
	mark = skip_darkspace(scan);
	scan = skip_whitespace(mark);
	if (*mark != 0) *mark = 0;
	start = string_to_u32 (field);

	if (*scan == 0) goto no_end;
	field = scan;
	mark = skip_darkspace(scan);
	scan = skip_whitespace(mark);
	if (*mark != 0) *mark = 0;
	end = string_to_u32 (field);

	if ((valCol == -1) || (_val == NULL))
		val = 1.0;
	else
		{
		for (col=3 ; col<=valCol ; col++)
			{
			if (*scan == 0) goto no_value;
			field = scan;
			mark = skip_darkspace(scan);
			scan = skip_whitespace(mark);
			}
		if (*mark != 0) *mark = 0;
		val = (valtype) string_to_valtype (field);
		}

	//////////
	// success
	//////////

	if (_chrom != NULL) *_chrom = chrom;
	if (_start != NULL) *_start = start;
	if (_end   != NULL) *_end   = end;
	if (_val   != NULL) *_val   = val;

	return true;

	//////////
	// failure exits
	//
	// nota bene: I'd like to report the line's contents with these messages,
	//            but the line has usually been modified before we get here
	//////////

missing_eol:
	fprintf (stderr, "problem at line %u, line is longer than internal buffer\n",
			 lineNumber-1);
	exit (EXIT_FAILURE);

no_chrom:
	fprintf (stderr, "problem at line %u, line contains no chromosome or begins with whitespace\n",
			 lineNumber);
	exit (EXIT_FAILURE);

no_start:
	fprintf (stderr, "problem at line %u, line contains no interval start\n"
	                 "(expected \"chromosome start end ...\", but there are fewer than 2 fields)\n",
			 lineNumber);
	exit (EXIT_FAILURE);

no_end:
	fprintf (stderr, "problem at line %u, line contains no interval end\n"
	                 "(expected \"chromosome start end ...\", but there are fewer than 3 fields)\n",
			 lineNumber);
	exit (EXIT_FAILURE);

no_value:
	fprintf (stderr, "problem at line %u, line contains no interval value\n"
	                 "(expected \"chromosome start end value\", but there are fewer than 4 fields)\n",
			 lineNumber);
	exit (EXIT_FAILURE);
	}

//----------
//
// report_intervals--
//	Report all intervals to a file.
//
//----------
//
// Arguments:
//	FILE*	f:				File to write to.
//	int		precision:		Number of digits of precision.
//	int		noOutputValues:	true  => just output intervals, without values
//	int		collapseRuns:	true  => collapse runs into a single interval
//							false => output every entry as its own interval
//	int		showUncovered:	how to report uncovered gaps are reported as zeros
//							  uncovered_hide => gaps are not reported
//							  uncovered_show => gaps are reported as zeros
//							  uncovered_NA   => gaps are reported as NAs
//	int		originOne:		true  => report intervals as origin-one, closed
//							false => report intervals as origin-zero, half-open
//
// Returns:
//	(nothing)
//
//----------

void report_intervals
   (FILE*		f,
	int			precision,
	int			noOutputValues,
	int			collapseRuns,
	int			showUncovered,
	int			originOne)
	{
	valtype*	v = NULL;
	spec*		chromSpec;
	u32			o;
	u32			start, outputStart, outputEnd, prevOutputEnd;
	valtype		val;
	u32			ix;
	int			active;

	if (originOne) o = 1;
	          else o = 0;

	for (chromSpec=chromsOfInterest ; chromSpec!=NULL ; chromSpec=chromSpec->next)
		{
		v = chromSpec->valVector;

		if (trackOperations)
			tracking_report ("output(%s)\n", chromSpec->chrom);

		active = (showUncovered != uncovered_hide);

		start = prevOutputEnd = 0;
		val = 0.0;
		for (ix=0 ; ix<chromSpec->length ; ix++)
			{
			// if this value is zero and we're hiding uncovered intervals,
			// output the previous interval (if there was one), and set the
			// current state to 'inactive'

			if ((v[ix] == 0) && (showUncovered != uncovered_show))
				{
				if ((active) && (ix != start))
					{
					outputStart = chromSpec->start+start;
					outputEnd   = chromSpec->start+ix;
					if ((showUncovered == uncovered_NA)
					 && (outputStart != prevOutputEnd))
						fprintf (f, "%s\t%d\t%d\tNA\n",
						            chromSpec->chrom, prevOutputEnd+o, outputStart);
					if (noOutputValues)
						fprintf (f, "%s\t%d\t%d\n",
						            chromSpec->chrom, outputStart+o, outputEnd);

					else
						fprintf (f, "%s\t%d\t%d\t" valtypeFmtPrec "\n",
						            chromSpec->chrom, outputStart+o, outputEnd,
						            precision, val);
					prevOutputEnd = outputEnd;
					}
				active = false;  start = 0;  val = 0.0;
				continue;
				}

			// if we weren't previously 'active', mark the start of this active
			// interval

			if (!active)
				{
				active = true;  start = ix;  val = v[ix];
				continue;
				}

			// if the new value is the same as the active run, and if we're
			// supposed to collapse each equal-value runs into single output
			// records, there's nothing else to do

			if ((v[ix] == val) && (collapseRuns))
				continue;

			// otherwise, we need to output the previous active interval and
			// mark the start of this new one 

			if (ix != start)
				{
				outputStart = chromSpec->start+start;
				outputEnd   = chromSpec->start+ix;
				if ((showUncovered == uncovered_NA)
				 && (outputStart != prevOutputEnd))
					fprintf (f, "%s\t%d\t%d\tNA\n",
								chromSpec->chrom, prevOutputEnd+o, outputStart);
				if (noOutputValues)
					fprintf (f, "%s\t%d\t%d\n",
								chromSpec->chrom, outputStart+o, outputEnd);
				else
					fprintf (f, "%s\t%d\t%d\t" valtypeFmtPrec "\n",
								chromSpec->chrom, outputStart+o, outputEnd,
								precision, val);
				prevOutputEnd = outputEnd;
				}
			active = true;  start = ix;  val = v[ix];
			}

		// show the final active interval (if there is one)

		if ((active) && (chromSpec->length != start))
			{
			outputStart = chromSpec->start+start;
			outputEnd   = chromSpec->start+chromSpec->length;
			if ((showUncovered == uncovered_NA)
			 && (outputStart != prevOutputEnd))
				fprintf (f, "%s\t%d\t%d\tNA\n",
							chromSpec->chrom, prevOutputEnd+o, outputStart);
			if (noOutputValues)
				fprintf (f, "%s\t%d\t%d\n",
							chromSpec->chrom, outputStart+o, outputEnd);
			else
				fprintf (f, "%s\t%d\t%d\t" valtypeFmtPrec "\n",
							chromSpec->chrom, outputStart+o, outputEnd,
							precision, val);
			prevOutputEnd = outputEnd;
			}
		else if ((showUncovered == uncovered_NA)
		      && (chromSpec->start+chromSpec->length != prevOutputEnd))
			{
			outputEnd = chromSpec->start+chromSpec->length;
			fprintf (f, "%s\t%d\t%d\tNA\n",
						chromSpec->chrom, prevOutputEnd+o, outputEnd);
			}
		}

	if (trackOperations)
		tracking_report ("output(--done--)\n");

	}

//----------
//
// read_all_chromosomes, write_all_chromosomes--
//	Read/write the current data set from/to a file.  Operators can use these as
//	a pair, to save and preserve an incomping dataset when the operation is
//	destructive.
//
// No claim is made about the file format.  Files written by one run of the
// program should not be expected to work when read into another run of the
// program.
//
// $$$ this could be improved by using a smarter file format
//
//----------
//
// Arguments:
//	char*	filename:	Name of the file to read from or write to.
//
// Returns:
//	(nothing)
//
//----------

// read_all_chromosomes--

void read_all_chromosomes
   (char*	filename)
	{
	FILE*	f;
	int		saveTrackOperations;

	f = fopen (filename, "rt");
	if (f == NULL) goto cant_open_file;

	if (trackOperations)
		fprintf (stderr, "read_all(%s)\n", filename);

	saveTrackOperations = trackOperations;
	trackOperations     = false;
	read_intervals (f,
	                /* val column */ 4-1,
	                /* origin one */ false,
	                /* overlap op */ ri_overlapSum,
	                /* clear      */ true, 0.0);
	trackOperations = saveTrackOperations;

	// success!

	fclose (f);
	return;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "can't open \"%s\" for reading\n",
	                 filename);
	exit (EXIT_FAILURE);
	}


// write_all_chromosomes--

void write_all_chromosomes
   (char*	filename)
	{
	FILE*	f;
	int		saveTrackOperations;

	f = fopen (filename, "wt");
	if (f == NULL) goto cant_open_file;

	if (trackOperations)
		fprintf (stderr, "write_all(%s)\n", filename);

	saveTrackOperations = trackOperations;
	trackOperations     = false;
	report_intervals (f,
	                  /* precision        */ 10,
	                  /* no output values */ false,
	                  /* collapse runs    */ true,
	                  /* show uncovered   */ false,
	                  /* origin one       */ false);
	trackOperations = saveTrackOperations;

	// success!

	fclose (f);
	return;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "can't open \"%s\" for writing\n",
	                 filename);
	exit (EXIT_FAILURE);
	}

// binary-file versions;  these have not been debugged, and they have a problem

//// read_all_chromosomes--
//
//void read_all_chromosomes
//   (char*	filename)
//	{
//	FILE*	f;
//	spec*	chromSpec;
//	u32		chromIx, bytesToRead, bytesRead;
//
//	f = fopen (filename, "rb");
//	if (f == NULL) goto cant_open_file;
//
//	if (trackOperations)
//		fprintf (stderr, "read_all(%s)\n", filename);
//
//	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
//		{
//		chromSpec = chromsSorted[chromIx];
//
//		bytesToRead = chromSpec->length * sizeof(valtype);
//		bytesRead   = fread (chromSpec->valVector, 1, bytesToRead, f);
//		if (bytesRead != bytesToRead) goto read_failure;
//		}
//
//	// success!
//
//	fclose (f);
//	return;
//
//	//////////
//	// failure exits
//	//////////
//
//cant_open_file:
//	fprintf (stderr, "can't open \"%s\" for reading\n",
//	                 filename);
//	exit (EXIT_FAILURE);
//
//read_failure:
//	fprintf (stderr, "problem reading from %s (attempted %d bytes, read %d)\n",
//	                 filename, bytesToRead, bytesRead);
//	}
//
//
//// write_all_chromosomes--
//
//void write_all_chromosomes
//   (char*	filename)
//	{
//	FILE*	f;
//	spec*	chromSpec;
//	u32		chromIx, bytesToWrite, bytesWritten;
//
//	f = fopen (filename, "wb");
//	if (f == NULL) goto cant_open_file;
//
//	if (trackOperations)
//		fprintf (stderr, "write_all(%s)\n", filename);
//
//	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
//		{
//		chromSpec = chromsSorted[chromIx];
//
//		bytesToWrite = chromSpec->length * sizeof(valtype);
//		bytesWritten = fwrite (chromSpec->valVector, 1, bytesToWrite, f);
//		if (bytesWritten != bytesToWrite) goto write_failure;
//		}
//
//	// success!
//
//	fclose (f);
//	return;
//
//	//////////
//	// failure exits
//	//////////
//
//cant_open_file:
//	fprintf (stderr, "can't open \"%s\" for writing\n",
//	                 filename);
//	exit (EXIT_FAILURE);
//
//write_failure:
//	fprintf (stderr, "problem writing to %s (attempted %d bytes, wrote %d)\n",
//	                 filename, bytesToWrite, bytesWritten);
//	}

//----------
//
// init_scratch_vectors, get_scratch_vector, release_scratch_vector, free_scratch_vectors--
//	Allocate and reuse scratch vectors.
//
//----------

static u32		scratchLength;
static svspec*	scratchVectorHead    = NULL;
static svispec*	scratchVectorIntHead = NULL;


static void init_scratch_vectors
   (u32 _scratchLength)
	{
	scratchLength        = _scratchLength;
	scratchVectorHead    = NULL;
	scratchVectorIntHead = NULL;
	}


valtype* get_scratch_vector
   (void)
	{
	svspec*	svSpec;

	for (svSpec=scratchVectorHead ; svSpec!=NULL ; svSpec=svSpec->next)
		{
		if (svSpec->inUse) continue;
		svSpec->inUse = true;
		//fprintf (stderr, "re-using scratch vector: %p\n", svSpec->vector);
		return svSpec->vector;
		}

	svSpec = malloc (sizeof(svspec));
	if (svSpec == NULL) goto cant_allocate_spec;
	//fprintf (stderr, "allocated %u bytes for scratch element: %p\n",
	//                 (u32) sizeof(svspec), svSpec);
	svSpec->next = scratchVectorHead;
	scratchVectorHead = svSpec;
	svSpec->inUse = true;

	svSpec->vector = (valtype*) calloc (scratchLength, sizeof(valtype));
	if (svSpec->vector == NULL) goto cant_allocate_scratch;
	//fprintf (stderr, "allocated %u bytes for scratch vector: %p\n",
	//                 (u32) (scratchLength*sizeof(valtype)), svSpec->vector);
	return svSpec->vector;

cant_allocate_spec:
	fprintf (stderr, "failed to allocate spec for scratch vector, %d bytes\n",
	                 (int) sizeof(svSpec));
	exit(EXIT_FAILURE);

cant_allocate_scratch:
	fprintf (stderr, "failed to allocate %d-base scratch vector, %d bytes per base\n",
	                 scratchLength, (int) sizeof(valtype));
	exit(EXIT_FAILURE);
	}


s32* get_scratch_ints
   (void)
	{
	svispec*	sviSpec;

	for (sviSpec=scratchVectorIntHead ; sviSpec!=NULL ; sviSpec=sviSpec->next)
		{
		if (sviSpec->inUse) continue;
		sviSpec->inUse = true;
		//fprintf (stderr, "re-using scratch vector: %p\n", sviSpec->vector);
		return sviSpec->vector;
		}

	sviSpec = malloc (sizeof(svispec));
	if (sviSpec == NULL) goto cant_allocate_spec;
	//fprintf (stderr, "allocated %u bytes for scratch element: %p\n",
	//                 (u32) sizeof(svispec), sviSpec);
	sviSpec->next = scratchVectorIntHead;
	scratchVectorIntHead = sviSpec;
	sviSpec->inUse = true;

	sviSpec->vector = (s32*) calloc (scratchLength, sizeof(s32));
	if (sviSpec->vector == NULL) goto cant_allocate_scratch;
	//fprintf (stderr, "allocated %u bytes for scratch vector: %p\n",
	//                 (u32) (scratchLength*sizeof(s32)), sviSpec->vector);
	return sviSpec->vector;

cant_allocate_spec:
	fprintf (stderr, "failed to allocate spec for int vector, %d bytes\n",
	                 (int) sizeof(sviSpec));
	exit(EXIT_FAILURE);

cant_allocate_scratch:
	fprintf (stderr, "failed to allocate %d-base int vector, %d bytes per base\n",
	                 scratchLength, (int) sizeof(u32));
	exit(EXIT_FAILURE);
	}


void release_scratch_vector
   (valtype*	v)
	{
	svspec*		svSpec;

	for (svSpec=scratchVectorHead ; svSpec!=NULL ; svSpec=svSpec->next)
		{
		if (svSpec->vector != v) continue;
		svSpec->inUse = false;
		//fprintf (stderr, "releasing scratch vector: %p\n", v);
		break;
		}
	}


void release_scratch_ints
   (s32*		v)
	{
	svispec*	sviSpec;

	for (sviSpec=scratchVectorIntHead ; sviSpec!=NULL ; sviSpec=sviSpec->next)
		{
		if (sviSpec->vector != v) continue;
		sviSpec->inUse = false;
		//fprintf (stderr, "releasing scratch vector: %p\n", v);
		break;
		}
	}


static void free_scratch_vectors
   (void)
	{
	svspec*		svSpec,  *nextSpec;
	svispec*	sviSpec, *nextISpec;

	for (svSpec=scratchVectorHead ; svSpec!=NULL ; svSpec=nextSpec)
		{
		nextSpec = svSpec->next;
		//fprintf (stderr, "about to free scratch element: %p\n", svSpec);
		//fprintf (stderr, "about to free scratch vector:  %p\n", svSpec->vector);
		if (svSpec->vector != NULL) free (svSpec->vector);
		free (svSpec);
		}
	scratchVectorHead = NULL;

	for (sviSpec=scratchVectorIntHead ; sviSpec!=NULL ; sviSpec=nextISpec)
		{
		nextISpec = sviSpec->next;
		//fprintf (stderr, "about to free scratch element: %p\n", sviSpec);
		//fprintf (stderr, "about to free scratch vector:  %p\n", sviSpec->vector);
		if (sviSpec->vector != NULL) free (sviSpec->vector);
		free (sviSpec);
		}
	scratchVectorIntHead = NULL;
	}

//----------
//
// init_named_globals, set_named_global, get_named_global, free_named_globals--
//	Maintain a collection of named global variables.
//
//----------

static namedglobal*	namedGlobalHead    = NULL;


static void init_named_globals
   (void)
	{
	namedGlobalHead = NULL;
	}


void set_named_global
   (char*			name,
	valtype			val)
	{
	namedglobal*	ng, *ngVar;
	u32				numBytes, nameOffset;

	if (dbgGlobals)
		fprintf (stderr, "set_named_global(%s," valtypeFmt ")\n", name, val);

	ngVar = NULL;
	for (ng=namedGlobalHead ; ng!=NULL ; ng=ng->next)
		{
		if ((ng->name == NULL) || (strcmp (name, ng->name) != 0)) continue;
		ngVar = ng;
		break;
		}

	if (ngVar == NULL)
		{
		numBytes   =  sizeof(namedglobal);
		nameOffset =  numBytes;
		numBytes   += strlen(name) + 1;

		ngVar = (namedglobal*) malloc (numBytes);
		if (ngVar == NULL) goto cant_allocate;
		ngVar->next     = namedGlobalHead;
		namedGlobalHead = ngVar;

		ngVar->name = ((char*) ngVar) + nameOffset;
		strcpy (ngVar->name, name);
		}

	ngVar->v = val;

	return;

	//////////
	// failure exits
	//////////

cant_allocate:
	fprintf (stderr, "failed to allocate named global \"%s\", %d bytes\n",
	                 name, numBytes);
	exit(EXIT_FAILURE);
	}


valtype get_named_global
   (char*			name,
	valtype			defaultVal)
	{
	namedglobal*	ng, *ngVar;

	if (dbgGlobals)
		fprintf (stderr, "get_named_global(%s) = ", name);

	ngVar = NULL;
	for (ng=namedGlobalHead ; ng!=NULL ; ng=ng->next)
		{
		if ((ng->name == NULL) || (strcmp (name, ng->name) != 0)) continue;
		ngVar = ng;
		break;
		}

	if (ngVar == NULL)
		{
		if (dbgGlobals)
			fprintf (stderr, valtypeFmt " (default)\n", defaultVal);
		return defaultVal;
		}
	else
		{
		if (dbgGlobals)
			fprintf (stderr, valtypeFmt "\n", ngVar->v);
		return ngVar->v;
		}

	}


int named_global_exists
   (char*			name,
	valtype*		v)
	{
	namedglobal*	ng, *ngVar;

	if (dbgGlobals)
		fprintf (stderr, "named_global_exists(%s) = ", name);

	ngVar = NULL;
	for (ng=namedGlobalHead ; ng!=NULL ; ng=ng->next)
		{
		if ((ng->name == NULL) || (strcmp (name, ng->name) != 0)) continue;
		ngVar = ng;
		break;
		}

	if (ngVar == NULL)
		{
		if (dbgGlobals)
			fprintf (stderr, " (not found)\n");
		return false;
		}
	else
		{
		if (dbgGlobals)
			fprintf (stderr, valtypeFmt "\n", ngVar->v);
		if (v != NULL) (*v) = ngVar->v;
		return true;
		}

	}


static void free_named_globals
   (void)
	{
	namedglobal*	ng,  *ngNext;

	for (ng=namedGlobalHead ; ng!=NULL ; ng=ngNext)
		{ ngNext = ng->next;  free (ng); }
	namedGlobalHead = NULL;
	}


void report_named_globals
   (FILE*	f,
	char*	indent)
	{
	namedglobal*	ng;
	int				nameLen, nameW;

	if (indent == NULL) indent = "";

	nameW = 1;
	for (ng=namedGlobalHead ; ng!=NULL ; ng=ng->next)
		{
		nameLen = strlen(ng->name);
		if (nameLen > nameW) nameW = nameLen;
		}
	if (nameW > 20) nameW = 20;

	for (ng=namedGlobalHead ; ng!=NULL ; ng=ng->next)
		fprintf (f, "%s%*s = " valtypeFmt "\n",
		            indent, nameW, ng->name, ng->v);
	}

//----------
//
// tracking_report--
//	Make a tracking report to stderr.
//
// Each report is written as a single line.  If the line does NOT end with a
// new-line, the next report will be written over it.
//
//----------

static char	trLineBuff[1001];
static int  trLineLen     = 0;
static int  trPrevLineLen = 0;
static int  trNewLine     = false;

void tracking_report (const char* format, ...)
	{
	va_list	args;

	va_start (args, format);
	trLineBuff[0] = 0;
	if (format != NULL)
		vsprintf (trLineBuff, format, args);
	va_end (args);

	trLineLen = strlen(trLineBuff);
	trNewLine = (trLineLen > 0) && (trLineBuff[trLineLen-1] == '\n');
	if (trNewLine) trLineBuff[--trLineLen] = 0;

	fprintf (stderr, "%s", trLineBuff);
	if (trPrevLineLen > trLineLen)
		fprintf (stderr, "%*s", trPrevLineLen - trLineLen, "");
	if (trNewLine) { fprintf (stderr, "\n");  trPrevLineLen = 0;         }
			  else { fprintf (stderr, "\r");  trPrevLineLen = trLineLen; }
	}

//----------
// [[-- comparison function for the standard c function qsort --]]
//
// valtype_ascending--
//	Compare two values, to sort them into increasing order.
//
// The implemantation is a modified version of one found here:
//	www.gnu.org/software/libc/manual/html_node/Comparison-Functions.html#Comparison-Functions
//
//----------
//
// Arguments:
//	const void* v1:	One value (really const valtype*).
//	const void* v2:	Another   (really const valtype*).
//
// returns:
//	> 0 => v1 is greater than v2.
//	= 0 => v1 and v2 are the same.
//	< 0 => v1 is less than v2.
//
//----------

int valtype_ascending
   (const void*	_v1,
	const void*	_v2)
	{
	const valtype* v1 = (const valtype*) _v1;
	const valtype* v2 = (const valtype*) _v2;

	return (*v1 > *v2) - (*v1 < *v2);
	}


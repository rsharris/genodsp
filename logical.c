// logical.c-- genodsp operators performing logical operations

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include "utilities.h"
#include "genodsp_interface.h"
#include "logical.h"

//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_binarize--
//	Apply a threshold filter, converting values to 1 or 0 depending on whether
//	they are above (or at) a threshold or below it.
//
//----------

// private dspop subtype

typedef struct dspop_binarize
	{
	dspop		common;			// common elements shared with all operators
	char*		thresholdVarName;
	valtype		threshold;
	int			tiesAbove;
	valtype		oneVal;
	valtype		zeroVal;
	} dspop_binarize;


// op_binarize_short--

void op_binarize_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "apply a threshold filter\n");
	}


// op_binarize_usage--

void op_binarize_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sApply a threshold filter. Values above or below a threshold are converted to\n",  indent);
	fprintf (f, "%s1 or 0, respectively.\n",                                                         indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [<threshold>] [options]\n", indent, name);
	fprintf (f, "%s  <threshold>              the threshold\n",                                      indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	fprintf (f, "%s  --threshold=<variable>   (T=) get threshold from named variable\n",             indent);
	fprintf (f, "%s  --ties:below             values equal to the threshold are converted to\n",     indent);
	fprintf (f, "%s                           zero\n",                                               indent);
	fprintf (f, "%s                           (this is the default)\n",                              indent);
	fprintf (f, "%s  --ties:above             values equal to the threshold are converted to one\n", indent);
	fprintf (f, "%s  --one=<value>            (O=) value for locations above the threshold\n",       indent);
	fprintf (f, "%s                           (default is 1.0)\n",                                   indent);
	fprintf (f, "%s  --zero=<value>           (Z=) value for locations below the threshold\n",       indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	}


// op_binarize_parse--

dspop* op_binarize_parse (char* name, int _argc, char** _argv)
	{
	dspop_binarize*	op;
	int				argc = _argc;
	char**			argv = _argv;
	char*			arg, *argVal;
	valtype			tempVal;
	int				haveThreshold;

	// allocate and initialize our control record

	op = (dspop_binarize*) malloc (sizeof(dspop_binarize));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->thresholdVarName = NULL;
	op->threshold        = 0.0;
	op->tiesAbove        = false;
	op->oneVal           = 1.0;
	op->zeroVal          = 0.0;

	// parse arguments

	haveThreshold = false;

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --threshold=<variable> or T=<variable>

		if ((strcmp_prefix (arg, "--threshold=") == 0)
		 || (strcmp_prefix (arg, "T=")           == 0)
		 || (strcmp_prefix (arg, "--T=")         == 0))
			{
			if (haveThreshold) goto more_than_one_threshold;
			op->thresholdVarName = copy_string (argVal);
			haveThreshold = true;
			goto next_arg;
			}

		// --ties:below and --ties:above

		if ((strcmp (arg, "--ties:below") == 0)
		 || (strcmp (arg, "--ties=below") == 0))
			{
			op->tiesAbove = false;
			goto next_arg;
			}

		if ((strcmp (arg, "--ties:above") == 0)
		 || (strcmp (arg, "--ties=above") == 0))
			{
			op->tiesAbove = true;
			goto next_arg;
			}

		// --one=<value> or O=<value>

		if ((strcmp_prefix (arg, "--one=") == 0)
		 || (strcmp_prefix (arg, "O=")      == 0)
		 || (strcmp_prefix (arg, "--O=")    == 0))
			{
			tempVal = string_to_valtype (argVal);
			op->oneVal = tempVal;
			goto next_arg;
			}

		// --zero=<value> or Z=<value>

		if ((strcmp_prefix (arg, "--zero=") == 0)
		 || (strcmp_prefix (arg, "Z=")      == 0)
		 || (strcmp_prefix (arg, "--Z=")    == 0))
			{
			tempVal = string_to_valtype (argVal);
			op->zeroVal = tempVal;
			goto next_arg;
			}

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <threshold>

		if (!haveThreshold)
			{
			op->threshold = string_to_valtype (arg);
			haveThreshold = true;
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_binarize));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_threshold:
	fprintf (stderr, "[%s] threshold specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_binarize_free--

void op_binarize_free (dspop* _op)
	{
	dspop_binarize*	op = (dspop_binarize*) _op;

	if (op->thresholdVarName != NULL) free (op->thresholdVarName);
	free (op);
	}


// op_binarize_apply--

void op_binarize_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_binarize*	op = (dspop_binarize*) _op;
	valtype			cutoffThresh = op->threshold;
	int				tiesAbove    = op->tiesAbove;
	valtype			oneVal       = op->oneVal;
	valtype			zeroVal      = op->zeroVal;
	int				ok;
	u32				ix;

	// if the threshold is a named variable, fetch it now;  note that we copy
	// the value from the named variable, then destroy our reference to the
	// named variable

	if (op->thresholdVarName != NULL)
		{
		ok = named_global_exists (op->thresholdVarName, &cutoffThresh);
		if (!ok) goto no_threshold;
		op->threshold = cutoffThresh;
		fprintf (stderr, "[%s] using %s = " valtypeFmt " as threshold\n",
		                 _op->name, op->thresholdVarName, cutoffThresh);
		free (op->thresholdVarName);
		op->thresholdVarName = NULL;
		}

	// process the vector

	if (tiesAbove)
		{
		for (ix=0 ; ix<vLen ; ix++)
			{ v[ix] = (v[ix] >= cutoffThresh)? oneVal : zeroVal; }
		}
	else
		{
		for (ix=0 ; ix<vLen ; ix++)
			{ v[ix] = (v[ix] > cutoffThresh)? oneVal : zeroVal; }
		}

	// success

	return;

	// failure

no_threshold:
	fprintf (stderr, "[%s] attempt to use %s as threshold failed (no such variable)\n",
	                 _op->name, op->thresholdVarName);
	exit(EXIT_FAILURE);
	}

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_or--
//	"Logical OR" an incoming set of interval values (read from a file) with the
//	current set.
//
//----------

// private dspop subtype

typedef struct dspop_or
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			originOne;
	} dspop_or;


// op_or_short--

void op_or_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "logical OR an incoming set of interval values with the current set\n");
	}


// op_or_usage--

void op_or_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sLogical OR an incoming set of interval values with the current set. The\n",       indent);
	fprintf (f, "%sresult is 1.0 if the entry is non-zero in either set. Otherwise the result is\n", indent);
	fprintf (f, "%s0.0.\n",                                                                          indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --value=<col>            input intervals contain a value in the specified\n",   indent);
	fprintf (f, "%s                           column\n",                                             indent);
	fprintf (f, "%s  --novalue                input intervals have no value (value given is 1)\n",   indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	}


// op_or_parse--

dspop* op_or_parse (char* name, int _argc, char** _argv)
	{
	dspop_or*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_or*) malloc (sizeof(dspop_or));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename  = NULL;
	op->valColumn = (int) get_named_global ("valColumn", 4-1);
	op->originOne = (int) get_named_global ("originOne", false);

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --value=<col>

		if ((strcmp (arg, "--novalue") == 0)
		 || (strcmp (arg, "--novalues") == 0)
		 || (strcmp (arg, "--value=none") == 0))
			{ op->valColumn = -1;  goto next_arg; }

		if (strcmp_prefix (arg, "--value=") == 0)
			{
			tempInt = string_to_int (argVal) - 1;
			if (tempInt == -1)
				chastise ("[%s] value column can't be 0 (\"%s\")\n",         name, arg);
			if (tempInt < 0)
				chastise ("[%s] value column can't be negative (\"%s\")\n",  name, arg);
			if (tempInt < 3)
				chastise ("[%s] value column can't be 1, 2 or 3 (\"%s\")\n", name, arg);
			op->valColumn = tempInt;
			goto next_arg;
			}

		// --origin=one, --origin=zero

		if ((strcmp (arg, "--origin=one") == 0)
		 || (strcmp (arg, "--origin=1")   == 0))
			{ op->originOne = true;  goto next_arg; }

		if ((strcmp (arg, "--origin=zero") == 0)
		 || (strcmp (arg, "--origin=0")    == 0))
			{ op->originOne = false;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <filename>

		if (op->filename == NULL)
			{
			op->filename = copy_string (arg);
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if (op->filename == NULL) goto filename_missing;

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_or));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_or_free--

void op_or_free (dspop* _op)
	{
	dspop_or*	op = (dspop_or*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_or_apply--

void op_or_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_or*	op = (dspop_or*) _op;
	char*			filename = op->filename;
	FILE*			f;
	char			lineBuffer[1001];
	char			prevChrom[1001];
	valtype*		v = NULL;
	char*			chrom;
	spec*			chromSpec;
	u32				start, end, o, adjStart, adjEnd;
	valtype			val;
	u32				ix, chromIx;
	int				ok;

	f = fopen (filename, "rt");
	if (f == NULL) goto cant_open_file;

	if (trackOperations)
		{
		for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
			{
			chromSpec = chromsSorted[chromIx];
			chromSpec->flag = false;
			}
		}

	// convert existing intervals to true/false (one/zero)

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		v = chromSpec->valVector;
		for (ix=0 ; ix<chromSpec->length ; ix++)
			{ if (v[ix] != 0.0) v[ix] = 1.0; }
		}

	// read intervals and values and "or" them

	if (op->originOne) o = 1;
	              else o = 0;

	prevChrom[0] = 0;
	chromSpec    = NULL;

	v = NULL;
	while (true)
		{
		ok = read_interval (f, lineBuffer, sizeof(lineBuffer), op->valColumn,
	                        &chrom, &start, &end, &val);
		if (!ok) break;
		if (val == 0.0) continue;

		if (strcmp (chrom, prevChrom) != 0)
			{
			v = NULL;
			chromSpec = find_chromosome_spec (chrom);
			if (chromSpec != NULL) v = chromSpec->valVector;
			safe_strncpy (prevChrom, chrom, sizeof(prevChrom)-1);
			}

		if (chromSpec == NULL) continue;

		if ((trackOperations) && (!chromSpec->flag))
			{
			fprintf (stderr, "%s(%s)\n", op->common.name, chrom);
			chromSpec->flag = true;
			}

		start -= o;
		adjStart = start;
		adjEnd   = end;

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

		// 'or' over this interval

		for (ix=adjStart ; ix<adjEnd ; ix++)
			v[ix] = 1.0;
		}

	// success

	fclose (f);
	return;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "[%s] can't open \"%s\" for reading\n",
	                 op->common.name, filename);
	exit (EXIT_FAILURE);

chrom_too_short:
	fprintf (stderr, "[%s] in \"%s\", %s %d %d is beyond the end of the chromosome (L=%d)\n",
	                 op->common.name, filename, chrom, start, end, chromSpec->length);
	exit (EXIT_FAILURE);
	}

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_and--
//	"Logical AND" an incoming set of interval values (read from a file) with the
//	current set.
//
// This amounts to clearing any intervals absent from the incoming set.
//
//----------

// private dspop subtype

typedef struct dspop_and
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			originOne;
	int			debug;
	} dspop_and;


// op_and_short--

void op_and_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "logical AND an incoming set of interval values with the current set\n");
	}


// op_and_usage--

void op_and_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sLogical AND an incoming set of interval values with the current set. The\n",      indent);
	fprintf (f, "%sresult is 1.0 if the entry is non-zero in both sets. Otherwise the result\n",     indent);
	fprintf (f, "%sis 0.0. This amounts to clearing any intervals absent from the incoming set.\n",  indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --value=<col>            input intervals contain a value in the specified\n",   indent);
	fprintf (f, "%s                           column\n",                                             indent);
	fprintf (f, "%s  --novalue                input intervals have no value (value given is 1)\n",   indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe input intervals are required to be sorted along each chromosome, and\n",      indent);
	fprintf (f, "%snon-overlapping)\n",                                                              indent);
	}


// op_and_parse--

dspop* op_and_parse (char* name, int _argc, char** _argv)
	{
	dspop_and*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_and*) malloc (sizeof(dspop_and));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename  = NULL;
	op->valColumn = (int) get_named_global ("valColumn", 4-1);
	op->originOne = (int) get_named_global ("originOne", false);

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --value=<col>

		if ((strcmp (arg, "--novalue") == 0)
		 || (strcmp (arg, "--novalues") == 0)
		 || (strcmp (arg, "--value=none") == 0))
			{ op->valColumn = -1;  goto next_arg; }

		if (strcmp_prefix (arg, "--value=") == 0)
			{
			tempInt = string_to_int (argVal) - 1;
			if (tempInt == -1)
				chastise ("[%s] value column can't be 0 (\"%s\")\n",         name, arg);
			if (tempInt < 0)
				chastise ("[%s] value column can't be negative (\"%s\")\n",  name, arg);
			if (tempInt < 3)
				chastise ("[%s] value column can't be 1, 2 or 3 (\"%s\")\n", name, arg);
			op->valColumn = tempInt;
			goto next_arg;
			}

		// --origin=one, --origin=zero

		if ((strcmp (arg, "--origin=one") == 0)
		 || (strcmp (arg, "--origin=1")   == 0))
			{ op->originOne = true;  goto next_arg; }

		if ((strcmp (arg, "--origin=zero") == 0)
		 || (strcmp (arg, "--origin=0")    == 0))
			{ op->originOne = false;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <filename>

		if (op->filename == NULL)
			{
			op->filename = copy_string (arg);
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if (op->filename == NULL) goto filename_missing;

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_and));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_and_free--

void op_and_free (dspop* _op)
	{
	dspop_and*	op = (dspop_and*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_and_apply--

void op_and_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_and*	op = (dspop_and*) _op;
	char*		filename = op->filename;
	FILE*		f;
	char		lineBuffer[1001];
	char		prevChrom[1001];
	valtype*	v = NULL;
	char*		chrom;
	spec*		chromSpec;
	u32			start, end, o, prevEnd, adjStart, adjEnd;
	valtype		val;
	u32			ix, chromIx;
	int			ok;

	f = fopen (filename, "rt");
	if (f == NULL) goto cant_open_file;

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		chromSpec->flag = false;
		}

	// convert existing intervals to true/false (one/zero)

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		v = chromSpec->valVector;
		for (ix=0 ; ix<chromSpec->length ; ix++)
			{ if (v[ix] != 0.0) v[ix] = 1.0; }
		}

	// read intervals and values and multiply by them

	if (op->originOne) o = 1;
	              else o = 0;

	prevChrom[0] = 0;
	chromSpec    = NULL;
	prevEnd      = 0;

	v = NULL;
	while (true)
		{
		ok = read_interval (f, lineBuffer, sizeof(lineBuffer), op->valColumn,
	                        &chrom, &start, &end, &val);
		if (!ok) break;
		if (val == 0.0) continue; // treat zero as a missing interval

		// if this interval is on a different chromosome, finish the previous
		// chromosome and then switch to this one

		if (strcmp (chrom, prevChrom) != 0)
			{
			// clear any gap at the end of the previous chromosome

			if (v != NULL)
				{
				for (ix=prevEnd ; ix<chromSpec->length ; ix++)
					v[ix] = 0.0;
				}

			v = NULL;
			chromSpec = find_chromosome_spec (chrom);
			if (chromSpec != NULL)
				{
				v = chromSpec->valVector;
				prevEnd = 0;
				if (chromSpec->flag) goto chrom_not_together;
				}

			safe_strncpy (prevChrom, chrom, sizeof(prevChrom)-1);
			}

		// if the current chromosome is not of any interest, ignore this
		// interval

		if (chromSpec == NULL) continue;

		// if this is the first interval on this chromosome, 'mark' the
		// chromosome

		if (!chromSpec->flag)
			{
			if (trackOperations) fprintf (stderr, "%s(%s)\n", op->common.name, chrom);
			chromSpec->flag = true;
			}

		// validate the interval, and (if necessary) shift it onto the vector

		start -= o;
		adjStart = start;
		adjEnd   = end;

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

		// make sure intervals are in sorted order along the chromosome

		if (adjStart < prevEnd)
			goto intervals_out_of_order;

		// clear any gap before this interval

		for (ix=prevEnd ; ix<adjStart ; ix++)
			v[ix] = 0.0;

		// "and" over this interval;  actually, there is nothing to do-- the
		// surrounding code will clear all intervals *not* present in the input

		; // (do nothing)

		prevEnd = adjEnd;
		}

	// clear any gap at the end of the last chromosome observed

	if (v != NULL)
		{
		for (ix=prevEnd ; ix<chromSpec->length ; ix++)
			v[ix] = 0.0;
		}

	// clear any chromosomes that weren't observed in the incoming set

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		if (chromSpec->flag) continue;

		chrom = chromSpec->chrom;
		v     = chromSpec->valVector;
		start = 0;
		end   = chromSpec->length;

		if (trackOperations)
			fprintf (stderr, "%s(%s,absent)\n", op->common.name, chrom);

		for (ix=start ; ix<end ; ix++)
			v[ix] = 0.0;
		}

	// success

	fclose (f);
	return;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "[%s] can't open \"%s\" for reading\n",
	                 op->common.name, filename);
	exit (EXIT_FAILURE);

chrom_too_short:
	fprintf (stderr, "[%s] in \"%s\", %s %d %d is beyond the end of the chromosome (L=%d)\n",
	                 op->common.name, filename, chrom, start, end, chromSpec->length);
	exit (EXIT_FAILURE);

chrom_not_together:
	fprintf (stderr, "[%s] in \"%s\", not all intervals on %s are together (%d..%d begins new group)\n",
	                 op->common.name, filename, chrom, start, end);
	exit (EXIT_FAILURE);

intervals_out_of_order:
	fprintf (stderr, "[%s] in \"%s\", intervals on %s are not sorted (%d..%d after %d)\n",
	                 op->common.name, filename, chrom, start, end, chromSpec->start+prevEnd);
	exit (EXIT_FAILURE);
	}

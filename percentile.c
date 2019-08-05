// percentile.c-- genodsp operators finding percentiles

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
#include "percentile.h"

#define percentileStepUnits 1000		// resolution for percentiles is in
										// .. thousandths of a percent

static void set_percentile_name    (char* varName, u32 percentile);
static void combine_sorted_vectors (valtype* v, u32 vLen, valtype* w, u32 wLen);
//static void sort_two_increasing_blocks (valtype* v, u32 vSplit, u32 vLen);

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_percentile--
//	Determine specified percentiles in the data and report them;  if min and/or
//	max values are specified, only values within this range are considered.
//
//	Note that the genome vectors *are* altered;  the only guarantee made about
//	the state the vectors are left in is that the set of values is preserved.
//
//----------

// private dspop subtype
// $$$ I should change percentileLo, Hi, Step to an array of values

typedef struct dspop_percentile
	{
	dspop		common;			// common elements shared with all operators
	char*		preseveFilename;
	char*		mapFilename;
	u32			percentileLo;	// percentile range (in thousandths of a
	u32			percentileHi;	// .. percent)
	u32			percentileStep;
	u32			windowSize;
	valtype		minAllowed;
	valtype		maxAllowed;
	int			valPrecision;
	int			reportForBash;
	int			quiet;
	int			debug;
	int			debugShowIndex;
	int			debugStage;
	} dspop_percentile;


// op_percentile_short--

void op_percentile_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "identify percentiles in the data\n");
	}


// op_percentile_usage--

void op_percentile_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sIdentify percentiles in the current data. This is a destructive operation;\n",    indent);
	fprintf (f, "%safterwards the data is in an undependable state unless --preserve is used.\n",    indent);
	fprintf (f, "%s(but see the special case noted below)\n",                                        indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <low>[..<high>] [options]\n", indent, name);
	fprintf (f, "%s  <low>[..<high>]          range of percentiles to report (e.g. 95..100)\n",      indent);
	fprintf (f, "%s  --step=<value>           step size of percentiles reported (e.g 0.1)\n",        indent);
	fprintf (f, "%s                           (default is 1.0)\n",                                   indent);
	fprintf (f, "%s  --window=<length>        (W=) size of window;  only the first entry in each\n", indent);
	fprintf (f, "%s                           window is considered;  this can greatly improve\n",    indent);
	fprintf (f, "%s                           speed;  but see caveat below\n",                       indent);
	fprintf (f, "%s  --min=<value>            minimum data value to consider;  lower values are\n",  indent);
	fprintf (f, "%s                           ignored\n",                                            indent);
	fprintf (f, "%s                           (default is -inf)\n",                                  indent);
	fprintf (f, "%s  --max=<value>            maximum data value to consider;  higher values are\n", indent);
	fprintf (f, "%s                           ignored\n",                                            indent);
	fprintf (f, "%s                           (default is +inf)\n",                                  indent);
	fprintf (f, "%s  --precision=<number>     number of digits to round percentile values to\n",     indent);
	fprintf (f, "%s  --preserve=<filename>    preserve the current data and restore it upon\n",      indent);
	fprintf (f, "%s                           completion;  note that it is often faster to just\n",  indent);
	fprintf (f, "%s                           recreate the data instead, since preservation adds\n", indent);
	fprintf (f, "%s                           a write/read to/from a file\n",                        indent);
	fprintf (f, "%s                           (by default, data is left in an erratic state)\n",     indent);
	fprintf (f, "%s  --map=<filename>         write percentile values to a file, suitable for\n",    indent);
	fprintf (f, "%s                           use as a mapping file\n",                              indent);
	fprintf (f, "%s  --quiet                  don't report percentile values to the console\n",      indent);
	fprintf (f, "%s  --report:bash            report percentile values to stdout in a form\n",       indent);
	fprintf (f, "%s                           suitable for defining a bash environment varirable\n", indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sIn addition to reporting percentiles to stderr, this function sets a named\n",    indent);
	fprintf (f, "%svariable for each percentile reported. The name is of the form percentile<p>\n",  indent);
	fprintf (f, "%swhere <p> has only whatever precision is necessary (e.g. 95 or 95.72).\n",        indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sA special case which does NOT leave the data in an undependable state is when\n", indent);
	fprintf (f, "%sthe range of percentiles consists only of 0 and/or 100.  In this case a\n",       indent);
	fprintf (f, "%sdifferent algorithm is used, simply scanning the signal for smallest or\n",       indent);
	fprintf (f, "%slargest value.\n",                                                                indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sIf --window is used, successive uses of this operation without preserving or\n",  indent);
	fprintf (f, "%srecreating the input signal can give unexpectedly different results.  This is\n", indent);
	fprintf (f, "%sbecause the operation reorders the signal (bringing the selected values to\n",    indent);
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sthe front), then sorts them, and a subsequent call 'sees' a different sample\n",  indent);
	fprintf (f, "%sof the signal (with potential for bias).\n",                                      indent);
	}


// op_percentile_parse--

dspop* op_percentile_parse (char* name, int _argc, char** _argv)
	{
	dspop_percentile*	op;
	int					argc = _argc;
	char**				argv = _argv;
	char*				arg, *argVal, *argVal2;
	int					tempInt;
	valtype				tempVal, tempLo, tempHi;
	int					haveRange;

	// allocate and initialize our control record

	op = (dspop_percentile*) malloc (sizeof(dspop_percentile));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->preseveFilename = NULL;
	op->mapFilename     = NULL;
	op->percentileLo    = 0.0; // (inital value is never used)
	op->percentileHi    = 0.0; // (inital value is never used)
	op->percentileStep  = percentileStepUnits;
	op->windowSize      = (u32) get_named_global ("windowSize", 1);
	op->minAllowed      = -valtypeMax;
	op->maxAllowed      =  valtypeMax;
	op->valPrecision    = (int) get_named_global ("valPrecision",  0);
	op->reportForBash   = false;
	op->quiet           = false;
	op->debug           = false;
	op->debugShowIndex  = false;
	op->debugStage      = 0;

	// parse arguments

	haveRange = false;

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --step=<value>

		if (strcmp_prefix (arg, "--step=") == 0)
			{
		set_step:
			tempVal = string_to_valtype (argVal);
			if (tempVal == 0)
				chastise ("[%s] step can't be zero (\"%s\")\n", name, arg);
			if (tempVal < 0)
				chastise ("[%s] step can't be negative (\"%s\")\n", name, arg);
			if (tempVal < .001) tempVal = .001;
			op->percentileStep = (u32) (percentileStepUnits*tempVal + .5);
			goto next_arg;
			}

		// --window=<length> or W=<length>

		if ((strcmp_prefix (arg, "--window=") == 0)
		 || (strcmp_prefix (arg, "W=")        == 0)
		 || (strcmp_prefix (arg, "--W=")      == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt == 0)
				tempInt = 1;
			if (tempInt < 0)
				chastise ("[%s] window size can't be negative (\"%s\")\n", name, arg);
			op->windowSize = (u32) tempInt;
			goto next_arg;
			}

		// --min=<value> or --max=<value>

		if (strcmp_prefix (arg, "--min=") == 0)
			{ op->minAllowed = string_to_valtype (argVal);  goto next_arg; }

		if (strcmp_prefix (arg, "--max=") == 0)
			{ op->maxAllowed = string_to_valtype (argVal);  goto next_arg; }

		// --precision=<col>

		if (strcmp_prefix (arg, "--precision=") == 0)
			{
			tempInt = string_to_int (argVal);
			if (tempInt < 0)
				chastise ("[%s] precision can't be negative (\"%s\")\n", name, arg);
			op->valPrecision = tempInt;
			goto next_arg;
			}

		// --preserve=<filename>

		if (strcmp_prefix (arg, "--preserve=") == 0)
			{
			if ((op->preseveFilename != NULL)
			 && (strcmp (argVal, op->preseveFilename) != 0))
				chastise ("[%s] can't specify two files for data preservation\n"
				          "(\"%s\" and \"%s\")", name, op->preseveFilename, argVal);
			op->preseveFilename = copy_string (argVal);
			goto next_arg;
			}

		// --map=<filename>

		if ((strcmp_prefix (arg, "--map=")     == 0)
		 || (strcmp_prefix (arg, "--mapping=") == 0))
			{
			if ((op->mapFilename != NULL)
			 && (strcmp (argVal, op->mapFilename) != 0))
				chastise ("[%s] can't specify two files for mapping\n"
				          "(\"%s\" and \"%s\")", name, op->mapFilename, argVal);
			op->mapFilename = copy_string (argVal);
			goto next_arg;
			}

		// --report:bash

		if ((strcmp (arg, "--report:bash") == 0)
		 || (strcmp (arg, "--bash")        == 0))
			{ op->reportForBash = true;  goto next_arg; }

		// --quiet

		if ((strcmp (arg, "--quiet")  == 0)
		 || (strcmp (arg, "--silent") == 0))
			{ op->quiet = true;  goto next_arg; }

		// --debug arguments

		if (strcmp (arg, "--debug") == 0)
			{ op->debug = true;  goto next_arg; }

		if (strcmp (arg, "--debug=index") == 0)
			{ op->debugShowIndex = true;  goto next_arg; }

		if (strcmp (arg, "--debug=collect") == 0)
			{ op->debugStage = 1;  goto next_arg; }

		if (strcmp (arg, "--debug=sort1") == 0)
			{ op->debugStage = 2;  goto next_arg; }

		if (strcmp (arg, "--debug=sort2") == 0)
			{ op->debugStage = 3;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <low>,<high>

		if ((!haveRange) && (strchr(arg,',') != NULL))
			{
			argVal = strchr(arg,',');
			*(argVal++) = 0;

			tempLo = string_to_valtype (arg);
			tempHi = string_to_valtype (argVal);

			if (tempLo > tempHi)
				{ tempVal = tempHi;  tempHi = tempLo;  tempLo = tempVal; }

			if      (tempLo <   0.0) op->percentileLo = 0;
			else if (tempLo > 100.0) op->percentileLo = 100*percentileStepUnits;
			else                     op->percentileLo = (int) (percentileStepUnits*tempLo + .5);

			if      (tempHi <   0.0) op->percentileHi = 0;
			else if (tempHi > 100.0) op->percentileHi = 100*percentileStepUnits;
			else                     op->percentileHi = (int) (percentileStepUnits*tempHi + .5);

			haveRange = true;
			if (op->percentileLo < op->percentileHi)
				op->percentileStep = op->percentileHi - op->percentileLo;
			else
				op->percentileStep = 1;
			goto next_arg;
			}

		// <low>[..<high>][by<step>]

		if (!haveRange)
			{
			argVal = strstr(arg,"..");
			if (argVal != NULL) { *argVal = 0;  argVal += 2; }

			argVal2 = NULL;
			if (argVal != NULL)
				{
				argVal2 = strstr(argVal,"by");
				if (argVal2 != NULL) { *argVal2 = 0;  argVal2 += 2; }
				}

			if (argVal == NULL)
				tempLo = tempHi = string_to_valtype (arg);
			else
				{
				tempLo = string_to_valtype (arg);
				tempHi = string_to_valtype (argVal);
				}

			if (tempLo > tempHi)
				{ tempVal = tempHi;  tempHi = tempLo;  tempLo = tempVal; }

			if      (tempLo <   0.0) op->percentileLo = 0;
			else if (tempLo > 100.0) op->percentileLo = 100*percentileStepUnits;
			else                     op->percentileLo = (int) (percentileStepUnits*tempLo + .5);

			if      (tempHi <   0.0) op->percentileHi = 0;
			else if (tempHi > 100.0) op->percentileHi = 100*percentileStepUnits;
			else                     op->percentileHi = (int) (percentileStepUnits*tempHi + .5);

			haveRange = true;
			if (argVal2 != NULL)
				{ argVal = argVal2;  goto set_step; }
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if (!haveRange) goto range_missing;

	if ((op->reportForBash) && (op->quiet))
		chastise ("[%s] Can't use both --report:bash and --quiet\n", name);

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_percentile));
	exit(EXIT_FAILURE);

range_missing:
	fprintf (stderr, "[%s] no range of percentiles was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_percentile_free--

void op_percentile_free (dspop* _op)
	{
	dspop_percentile*	op = (dspop_percentile*) _op;

	if (op->preseveFilename != NULL) free (op->preseveFilename);
	if (op->mapFilename     != NULL) free (op->mapFilename);
	free (op);
	}


// op_percentile_apply--

void op_percentile_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_percentile*	op = (dspop_percentile*) _op;
	FILE*		mapF = NULL;
	char		varName[100];
	u32			percentileLo   = op->percentileLo;
	u32			percentileHi   = op->percentileHi;
	u32			percentileStep = op->percentileStep;
	valtype		minAllowed     = op->minAllowed;
	valtype		maxAllowed     = op->maxAllowed;
	u32			windowSize     = op->windowSize;
	int			valPrecision   = op->valPrecision;
	spec*		chromSpec, *dstChromSpec;
	valtype*	v, *dstV;
	u32			numValues, numValuesInLastChrom;
	u32			percentileNext, percentileIx, numValuesPassed;
	valtype		val;
	u32			ix, chromIx, dstChromIx, lastChromIx, dstIx, dstVLen;
	float		pPct;
	u32			pIx;
	valtype		pVal, minVal, maxVal;

	if (op->debug)
		{
		fprintf (stderr, "percentile(" valtypeFmt, percentileLo  /(float)percentileStepUnits);
		fprintf (stderr, "," valtypeFmt,           percentileHi  /(float)percentileStepUnits);
		fprintf (stderr, "," valtypeFmt,           percentileStep/(float)percentileStepUnits);
		if      (minAllowed == -valtypeMax)  fprintf (stderr, ",-inf");
		else if (minAllowed ==  valtypePuny) fprintf (stderr, ",1/inf");
		else                                 fprintf (stderr, "," valtypeFmt, minAllowed);
		if      (maxAllowed ==  valtypeMax)  fprintf (stderr, ",+inf");
		else if (minAllowed == -valtypePuny) fprintf (stderr, ",-1/inf");
		else                                 fprintf (stderr, "," valtypeFmt, maxAllowed);
		fprintf (stderr, ")\n");
		}

	// handle special case of percentile 0

	if ((op->percentileLo == 0)
	 && (op->percentileHi == 0)
	 && (op->mapFilename == NULL))
		{
		minVal = valtypeMax;
		numValues = 0;

		for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
			{
			chromSpec = chromsSorted[chromIx];
			v = chromSpec->valVector;

			for (ix=0 ; ix<chromSpec->length ; ix+=windowSize)
				{
				if (v[ix] < minAllowed) continue;
				if (v[ix] > maxAllowed) continue;

				if (v[ix] < minVal) minVal = v[ix];
				numValues++;
				}
			}

		if (numValues == 0) goto no_values;

		set_percentile_name (varName, op->percentileLo);
		set_named_global    (varName, minVal);

		return;
		}

	// handle special case of percentile 100

	if ((op->percentileLo == 100*percentileStepUnits)
	 && (op->percentileHi == 100*percentileStepUnits)
	 && (op->mapFilename == NULL))
		{
		maxVal = -valtypeMax;
		numValues = 0;

		for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
			{
			chromSpec = chromsSorted[chromIx];
			v = chromSpec->valVector;

			for (ix=0 ; ix<chromSpec->length ; ix+=windowSize)
				{
				if (v[ix] < minAllowed) continue;
				if (v[ix] > maxAllowed) continue;

				if (v[ix] > maxVal) maxVal = v[ix];
				numValues++;
				}
			}

		if (numValues == 0) goto no_values;

		set_percentile_name (varName, op->percentileHi);
		set_named_global    (varName, maxVal);

		return;
		}

	// handle special case of percentiles 0 and 100

	if ((op->percentileLo == 0)
	 && (op->percentileHi == 100*percentileStepUnits)
	 && (op->mapFilename == NULL))
		{
		maxVal = -valtypeMax;
		minVal = valtypeMax;
		numValues = 0;

		for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
			{
			chromSpec = chromsSorted[chromIx];
			v = chromSpec->valVector;

			for (ix=0 ; ix<chromSpec->length ; ix+=windowSize)
				{
				if (v[ix] < minAllowed) continue;
				if (v[ix] > maxAllowed) continue;

				if (v[ix] > maxVal) maxVal = v[ix];
				if (v[ix] < minVal) minVal = v[ix];
				numValues++;
				}
			}

		if (numValues == 0) goto no_values;

		set_percentile_name (varName, op->percentileHi);
		set_named_global    (varName, maxVal);
		set_percentile_name (varName, op->percentileLo);
		set_named_global    (varName, minVal);

		return;
		}

	// preserve the input data

	if (op->preseveFilename != NULL)
		write_all_chromosomes (op->preseveFilename);

	// if we're to write percentiles to a mapping file, open the file

	mapF = NULL;
	if (op->mapFilename != NULL)
		mapF = fopen (op->mapFilename, "wt");

	// take an initial pass thru the genome (looking only at the first value in
	// each window) and move all qualifying values to the "front";  we preserve
	// the overall content of the genome, effectively only shuffling the values

	numValues = 0;

	dstChromIx   = 0;
	dstChromSpec = chromsSorted[dstChromIx];
	dstV         = dstChromSpec->valVector;
	dstVLen      = dstChromSpec->length;
	dstIx        = 0;

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		v = chromSpec->valVector;

		for (ix=0 ; ix<chromSpec->length ; ix+=windowSize)
			{
			if (v[ix] < minAllowed) continue;
			if (v[ix] > maxAllowed) continue;
			numValues++;

			if (dstIx >= dstVLen)
				{
				// hop over to the next destination vector
				dstChromSpec = chromsSorted[++dstChromIx];
				dstV         = dstChromSpec->valVector;
				dstVLen      = dstChromSpec->length;
				dstIx        = 0;
				}

			// move the qualifying value to the front by swapping it with the
			// first remaining non-qualifying value
			val = dstV[dstIx];  dstV[dstIx] = v[ix];  v[ix] = val;
			dstIx++;
			}
		}

	lastChromIx          = dstChromIx;
	numValuesInLastChrom = dstIx;

	if (numValues == 0) goto no_values;

	percentileNext  = percentileLo;
	percentileIx    = (u32) (((u64) numValues) * percentileNext / (100.0*percentileStepUnits));
	numValuesPassed = 0;

	if (op->debug)
		fprintf (stderr, "  numValues=%u (%u in %s), percentileIx=%u\n",
		                 numValues,
		                 numValuesInLastChrom, chromsSorted[lastChromIx]->chrom,
		                 percentileIx);

	if (op->debugStage == 1)
		return;

	// sort the qualifying values, from low to high (only those at the front,
	// where the qualifying values are);  first we make a pass through the
	// vectors, sorting each individually;  then we make an outer pass through
	// the vectors, filling that vector with the lowest values;  the fill is
	// accomplished by an inner pass along the vectors, trading low values in
	// from the inner vector with high values in the outer vector;  at each of
	// these inner steps, the values within each vector reamin sorted
	//
	// once we have sorted enough values to know the desired percentile, we
	// stop

	for (chromIx=0 ; chromIx<=lastChromIx ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];

		v = chromSpec->valVector;
		if (chromIx < lastChromIx) vLen = chromSpec->length;
		                      else vLen = numValuesInLastChrom;

		if (op->debug) fprintf (stderr, "sorting %s (%u)\n", chromSpec->chrom, vLen);
		qsort (v, vLen, sizeof(valtype), valtype_ascending);
		}

	for (chromIx=0 ; chromIx<=lastChromIx ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];

		if ((op->debugStage == 2) && (chromIx == 1)) goto done;
		if ((op->debugStage == 3) && (chromIx == 2)) goto done;

		v = chromSpec->valVector;
		if (chromIx < lastChromIx) vLen = chromSpec->length;
		                      else vLen = numValuesInLastChrom;

		if (op->debug) fprintf (stderr, "resolving %s (%u)\n", chromSpec->chrom, vLen);

		if (chromIx < lastChromIx)
			{
			for (dstChromIx=chromIx+1 ; dstChromIx<=lastChromIx ; dstChromIx++)
				{
				dstChromSpec = chromsSorted[dstChromIx];

				dstV = dstChromSpec->valVector;
				if (dstChromIx < lastChromIx) dstVLen = dstChromSpec->length;
				                         else dstVLen = numValuesInLastChrom;

				if (op->debug)
					fprintf (stderr, "  merging %s and %s (%u)\n",
					                 chromSpec->chrom, dstChromSpec->chrom, dstVLen);
				combine_sorted_vectors(v, vLen, dstV, dstVLen);
				}
			}

		// report whatever desired percentiles are in this vector;  if the
		// final percentile is in this vector, we're done;  otherwise, reduce
		// the next-percentile index and continue on to the next outer vector

		while (percentileIx < vLen)
			{
			pPct = percentileNext / ((float)percentileStepUnits);
			pIx  = numValuesPassed + percentileIx;
			pVal = v[percentileIx];

			set_percentile_name (varName, percentileNext);
			set_named_global    (varName, pVal);

			if (op->reportForBash)
				fprintf (stdout, "%s=" valtypeFmtPrec " # bash command\n",
				                 varName, valPrecision, pVal);
			else if (!op->quiet)
				{
				fprintf (stderr, "percentile %.3f is ", pPct);
				if (op->debugShowIndex) fprintf (stderr, "[%u] ", pIx);
				fprintf (stderr, valtypeFmtPrec "\n", valPrecision, pVal);
				}

			if (mapF != NULL)
				fprintf (mapF, valtypeFmtPrec " %.3f\n", valPrecision, pVal, pPct);

			percentileNext += percentileStep;
			if (percentileNext > percentileHi) goto done;
			percentileIx =  (u32) (((u64) numValues) * percentileNext / (100.0*percentileStepUnits));
			percentileIx -= numValuesPassed;
			}

		percentileIx    -= vLen;
		numValuesPassed += vLen;

		if (chromIx == lastChromIx)
			{
			if (percentileIx == 0)  // special case for 100th percentile
				{
				pPct = percentileNext / ((float)percentileStepUnits);
				pIx  = numValuesPassed + percentileIx - 1;
				pVal = v[vLen-1];

				set_percentile_name (varName, percentileNext);
				set_named_global    (varName, pVal);

				if (!op->quiet)
					{
					fprintf (stderr, "percentile %.3f is ", pPct);
					if (op->debugShowIndex) fprintf (stderr, "[%u] ", pIx);
					fprintf (stderr, valtypeFmtPrec "\n", valPrecision, pVal);
					}
				if (mapF != NULL)
					fprintf (mapF, valtypeFmtPrec " %.3f\n", valPrecision, pVal, pPct);
				}
			else
				goto internal_error;
			}
		}
done:

	// restore the original data;  note that we detroy (but don't delete) the
	// scratch file once we are done with it

	if (op->preseveFilename != NULL)
		{
		FILE* f;

		read_all_chromosomes (op->preseveFilename);

		f = fopen (op->preseveFilename, "wb");  // destroy the scratch file
		if (f != NULL) fclose (f);
		}

	// success

	if (mapF != NULL) fclose (mapF);
	return;

	//////////
	// warning exits
	//////////

no_values:
	fprintf (stderr, "[%s] percentile can't be computed;  no input values meet the criteria\n",
	                 op->common.name);
	if (mapF != NULL) fclose (mapF);
	return;

	//////////
	// failure exits
	//////////

internal_error:
	fprintf (stderr, "[%s]internal error;  "
	                 "got to end of genome without finding the next percentile index\n",
	                 op->common.name);
	exit (EXIT_FAILURE);
	}


// set_percentile_name--

static void set_percentile_name (char* varName, u32 percentile)
	{
	float	pPct;
	u32		denom;
	int		precision;

	if (percentile % percentileStepUnits == 0)
		{
		sprintf (varName, "percentile%d", percentile / percentileStepUnits);
		return;
		}

	pPct = percentile / ((float)percentileStepUnits);

	for (precision=1,denom=percentileStepUnits/10 ; denom>=1 ; precision++,denom/=10)
		{
		if (percentile % denom == 0)
			{
			sprintf (varName, "percentile%.*f", precision, pPct);
			return;
			}
		}

	sprintf (varName, "percentile%f", pPct);
	}

//----------
//
// combine_sorted_vectors--
//	Combine vectors v and w so that v[i] <= w[j] for all i,j.  This is
//	accomplished by exchanging the N lowest elements from the head of w with
//	the N highest elements from the tail of v.  N is determined by a pre-scan
//	of the data.
//
//	That leaves us with each block consisting of two increasing partitions.
//	We then sort each block in place by merging its two partitions.
//
//	Diagramatically, we do this:
//
//		   +----+    +----+         +----+    +----+         +-------+    +-------+
//		v: | V1 | w: | W1 | ---> v: | V1 | w: | V3 | ---> v: | V1,W1 | w: | V3,W3 |
//		   | .. |    | .. |         | .. |    | .. |         |  ..   |    |  ..   |
//		   | .. |    | W2 |         | .. |    | V4 |         |  ..   |    |  ..   |
//		   | .. |    | W3 |         | .. |    | W3 |         |  ..   |    |  ..   |
//		   | .. |    | .. |         | .. |    | .. |         |  ..   |    |  ..   |
//		   | V2 |    | .. |         | V2 |    | .. |         |  ..   |    |  ..   |
//		   | V3 |    | .. |         | W1 |    | .. |         |  ..   |    |  ..   |
//		   | .. |    | W4 |         | .. |    | W4 |         |  ..   |    | V4,W4 |
//		   | V4 |    +----+         | W2 |    +----+         | V2,W2 |    +-------+
//		   +----+                   +----+                   +-------+
//
//----------
//
// Arguments:
//	valtype*	v:		The first block.
//	u32			vLen:	The length of the first block.
//	valtype*	w:		The second block.
//	u32			wLen:	The length of the second block.
//
// Returns:
//	(nothing)
//
//----------

static void combine_sorted_vectors
   (valtype*	v,
	u32			vLen,
	valtype*	w,
	u32			wLen)
	{
	u32			vIx, wIx, vSplit, wSplit;
	valtype		val;

	if ((vLen == 0) || (wLen == 0)) return;		// nothing to combine
	if (v[vLen-1] <= w[0]) return;				// already sorted

	// identify the biggest N such that N elements from the head of w are all
	// lower than N elements from the tail of v

	wIx = 0;
	vIx = vLen-1;
	while ((wIx+1 < wLen) && (vIx > 0))
		{
		wIx++;  vIx--;
		if (v[vIx] > w[wIx]) continue;
		wIx--;  vIx++;
		break;
		}

	wSplit = wIx;		// need to swap w[0..wSplit] with v[vSplit..vLen-1]
	vSplit = vIx;

	// exchange the N elements from the head of w with the N elements from the
	// tail of v, maintaining the order within each of those sublists

	for (wIx=0,vIx=vSplit ; wIx<=wSplit ; wIx++,vIx++)
		{ val = w[wIx];  w[wIx] = v[vIx];  v[vIx] = val; }

	// perform an in-place sort on each block
	// $$$ note that we might be able to take advantage of the fact that each
	// $$$ .. of the blocks is partitioned into sort increasing sub-blocks

	qsort (v, vLen, sizeof(valtype), valtype_ascending);
	qsort (w, wLen, sizeof(valtype), valtype_ascending);
	}

//----------
//
// sort_two_increasing_blocks--
//	Perform an in-place sort on a block that consists of two sorted partitions.
//	The algorithm is something like a merge sort, combined with a series of
//	insertion sorts.  We iteratively pull the ith smallest item to correct
//	position in the top partition.  If this item comes from the bottom partition
//	the ith item in the top partition is displaced and moved into the bottom
//	partition, but keeping that partition ordered.
//
//----------
//
// Arguments:
//	valtype*	v:		The block containing the two partitions.
//	u32			vSplit:	The length of the first partition in the block.
//	u32			vLen:	The total length of the block.
//
// Returns:
//	(nothing)
//
//----------

// This has note been debugged, and is known to have a problem
//
//static void sort_two_increasing_blocks
//   (valtype*	v,
//	u32			vSplit,
//	u32			vLen)
//	{
//	u32			topIx, botIx, ix;
//	valtype		topVal, botVal;
//
//	if (vSplit == 0)    return;
//	if (vSplit == vLen) return;
//
////…… fprintf (stderr, "vSplit=%u vLen=%u\n", vSplit, vLen);
////…… fprintf (stderr, "  v[%u]=" valtypeFmt "\n", vSplit-1, v[vSplit-1]);
////…… fprintf (stderr, "  v[%u]=" valtypeFmt "\n", vSplit,   v[vSplit]);
//
//	botIx = vSplit;
//	for (topIx=0 ; topIx<botIx ; topIx++)
//		{
//		// compare the first items in both blocks;  if the item in the top
//		// block is lower (or same) we keep it and loop back for the next item
//
////…… fprintf (stderr, "  ---\n");
////…… fprintf (stderr, "  top v[%u]=" valtypeFmt "\n", topIx, v[topIx]);
////…… fprintf (stderr, "  bot v[%u]=" valtypeFmt "\n", botIx, v[botIx]);
//
//		topVal = v[topIx];
//		botVal = v[botIx];
//		if (topVal <= botVal) continue;
//
//		// otherwise the item in the bottom block is lower;  we move it to the
//		// top block, then insert the item from the top block into the bottom
//		// block
//
//		v[topIx] = botVal;
//
//		for (ix=botIx ; ix+1<vLen ; ix++)
//			{
//			if (topVal <= v[ix+1]) break;
//			v[ix] = v[ix+1];
//			}
////…… fprintf (stderr, "  --> v[%u]=" valtypeFmt "\n", ix, v[ix]);
//		v[ix] = topVal;
//		}
//
//	}

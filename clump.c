// clump.c-- genodsp operators performing clump search

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
#include "clump.h"

// private dspop subtype, used by both operators in this module--

typedef struct dspop_clump
	{
	dspop		common;			// common elements shared with all operators
	char*		averageVarName;
	valtype		average;
	u32			minLength;
	valtype		oneVal;
	valtype		zeroVal;
	int			debug;
	int			debugDetail;
	int			progress;
	} dspop_clump;

// prototypes

static void clump_search (dspop* _op, char* vName, u32 vLen, valtype* v,
                          int aboveThresh);

//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_clump--
//	Apply a clump-finder, identifying intervals with an average above (or for
//	op_skimp, below) a specified threshold.
//
// We strive to find all intervals that (a) have an average value above or
// equal to the specified threshold, and (b) are at least as long as the
// specified length.  For op_skimp, we strive to find all intervals that have
// an average value below or equal to the threshold.
//
// Note that, mathematically, such intervals may overlap.  We effectively merge
// all overlaps, and then we trim non-contributing values from the ends of the
// intervals (i.e. values that are less than the specified average).  Because
// of the overlap, it is possible that a reported interval has an average on
// the wrong side of the threshold.  But all bases in that interval will be in
// *some* interval that satisfies the criteria (they just might not all be in
// the *same* interval).
//
// Because of the trimming, we may report intervals shorter than the minimum
// length.
//
// The algorithm here is inspired by [1], though there are many differences
// from the algorithm presented there.
//
// References:
//	[1]	 Allison, Lloyd. "Longest biased interval and longest non-negative sum
//		 interval." Bioinformatics 19.10 (2003): 1294-1295.
//	[1b] www.csse.monash.edu.au/~lloyd/tildeProgLang/Java2/Biased/Biased.java

//----------

// op_clump_short--

void op_clump_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "find intervals with an average above some threshold\n");
	}


// op_clump_usage--

static void op_clump_usage_args (char* name, FILE* f, char* indent);

void op_clump_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sFind intervals with an average above some threshold. Bases in such\n",       indent);
	fprintf (f, "%sintervals are replaced with 1s; other bases with 0s.\n",                          indent);
	fprintf (f, "%s\n",                                                                              indent);
	fprintf (f, "%sGiven length L and threshold T, we identify intervals with length at least L\n",  indent);
	fprintf (f, "%ssuch that the average over that interval is at least T.  After identifying\n",    indent);
	fprintf (f, "%ssuch intervals, the ends are then trimmed to remove signal lower than T\n",       indent);
	fprintf (f, "%s(since these are only reducing the average).  Thus the resulting intervals\n",    indent);
	fprintf (f, "%s(actually runs of 1 rather than intervals) may be shorter than L.\n",             indent);
	op_clump_usage_args (name, f, indent);
	}

static void op_clump_usage_args (char* name, FILE* f, char* indent)
	{
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [<average>] [options]\n", indent, name);
	fprintf (f, "%s  <average>                the threshold\n",                                      indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	fprintf (f, "%s  --average=<variable>     (T=) get average from named variable\n",                                      indent);
	fprintf (f, "%s  --length=<length>        (L=) minimum length of qualifying interval\n",         indent);
	fprintf (f, "%s                           (default is 100)\n",                                   indent);
	fprintf (f, "%s  --one=<value>            (O=) value to fill qualifying intervals\n",            indent);
	fprintf (f, "%s                           (default is 1.0)\n",                                   indent);
	fprintf (f, "%s  --zero=<value>           (Z=) zero value to fill non-qualifying intervals\n",   indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	}

// op_clump_parse--

dspop* op_clump_parse (char* name, int _argc, char** _argv)
	{
	dspop_clump*	op;
	int				argc = _argc;
	char**			argv = _argv;
	char*			arg, *argVal;
	int				tempInt;
	valtype			tempVal;
	int				haveAverage;

	// allocate and initialize our control record

	op = (dspop_clump*) malloc (sizeof(dspop_clump));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->averageVarName = NULL;
	op->average        = 0.0;
	op->minLength      = 100;
	op->oneVal         = 1.0;
	op->zeroVal        = 0.0;
	op->debug          = false;
	op->debugDetail    = false;
	op->progress       = 0;

	// parse arguments

	haveAverage = false;

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --average=<variable> or T=<variable>

		if ((strcmp_prefix (arg, "--average=") == 0)
		 || (strcmp_prefix (arg, "T=")         == 0)
		 || (strcmp_prefix (arg, "--T=")       == 0))
			{
			if (haveAverage) goto more_than_one_average;
			op->averageVarName = copy_string (argVal);
			haveAverage = true;
			goto next_arg;
			}

		// --length=<length> or L=<length>

		if ((strcmp_prefix (arg, "--length=") == 0)
		 || (strcmp_prefix (arg, "L=")        == 0)
		 || (strcmp_prefix (arg, "--L=")      == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt == 0)
				chastise ("[%s] minimum length can't be zero (\"%s\")\n", name, arg);
			if (tempInt < 0)
				chastise ("[%s] minimum length can't be negative (\"%s\")\n", name, arg);
			op->minLength = (u32) tempInt;
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

		// --debug arguments

		if (strcmp (arg, "--debug") == 0)
			{ op->debug = true;  goto next_arg; }

		if (strcmp (arg, "--debug=detail") == 0)
			{ op->debugDetail = true;  goto next_arg; }

		if (strcmp_prefix (arg, "--progress=") == 0)
			{
			op->progress = string_to_unitized_int (argVal, true);
			goto next_arg;
			}

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <average>

		if (!haveAverage)
			{
			op->average = string_to_valtype (arg);
			haveAverage = true;
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
	                 name, (int) sizeof(dspop_clump));
	exit(EXIT_FAILURE);

more_than_one_average:
	fprintf (stderr, "[%s] average threshold specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_clump_free--

void op_clump_free (dspop* _op)
	{
	dspop_clump*	op = (dspop_clump*) _op;

	if (op->averageVarName != NULL) free (op->averageVarName);
	free (op);
	}


// op_clump_apply--

void op_clump_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{ clump_search (_op,vName,vLen,v,/*above thresh*/ true); }


//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_skimp--
//	Apply a clump-finder, identifying intervals with an average below a
//	specified threshold.
//
// For more details, see the header for op_clump.
//
//----------

// op_skimp_short--

void op_skimp_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "find intervals with an average below some threshold\n");
	}


// op_skimp_usage--

void op_skimp_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sFind intervals with an average below some threshold. Bases in such\n",       indent);
	fprintf (f, "%sintervals are replaced with 1s; other bases with 0s.\n",                          indent);
	fprintf (f, "%s\n",                                                                              indent);
	fprintf (f, "%sGiven length L and threshold T, we identify intervals with length at least L\n",  indent);
	fprintf (f, "%ssuch that the average over that interval is no more than T.  After\n",            indent);
	fprintf (f, "%sidentifying such intervals, the ends are then trimmed to remove signal higher\n", indent);
	fprintf (f, "%sthan T (since these are only increasing the average).  Thus the resulting\n",     indent);
	fprintf (f, "%sintervals (actually runs of 1 rather than intervals) may be shorter than L.\n",   indent);
	op_clump_usage_args (name, f, indent);
	}


// op_skimp_parse--

dspop* op_skimp_parse (char* name, int _argc, char** _argv)
	{ return op_clump_parse (name, _argc, _argv); }


// op_skimp_free--

void op_skimp_free (dspop* op)
	{
	free (op);
	}


// op_skimp_apply--

void op_skimp_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{ clump_search (_op,vName,vLen,v,/*above thresh*/ false); }


//----------
//
// clump_search--
//	Perform the clump-finding.
//
//----------

static void clump_search
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v),
	arg_dont_complain(int		aboveThresh))
	{
	dspop_clump*	op = (dspop_clump*) _op;
	valtype			targetAvg = op->average;
	u32				minLength = op->minLength;
	valtype			oneVal    = op->oneVal;
	valtype			zeroVal   = op->zeroVal;
	valtype*		s         = get_scratch_vector();
	valtype*		minSums   = get_scratch_vector();
	u32*			minWhere  = (u32*) get_scratch_ints();
	int				ok;
	u32				ix, iy, scanIx;
	valtype			minSum, valSum, val;
	u32				numMinSums, minScan, minIx;
	u32				start, end, prevStart, prevEnd;
	int				allMonotonic;

	// if the threshold is a named variable, fetch it now;  note that we copy
	// the value from the named variable, then destroy our reference to the
	// named variable

	if (op->averageVarName != NULL)
		{
		ok = named_global_exists (op->averageVarName, &targetAvg);
		if (!ok) goto no_threshold;
		op->average = targetAvg;
		fprintf (stderr, "[%s] using %s = " valtypeFmt " as threshold\n",
		                 _op->name, op->averageVarName, targetAvg);
		free (op->averageVarName);
		op->averageVarName = NULL;
		}

	// perform a pre-scan do determine whether the sum in the algorithm would
	// be strictly decreasing over the entire vector;  this would cause the
	// minSums and minWhere arrays to overflow;  instead, if this case happens
	// we know there are no clump intervals and we can just erase the entire
	// vector and quit

	allMonotonic = true;
	for (ix=0 ; ix<vLen ; ix++)
		{
		if (aboveThresh) val = v[ix] - targetAvg;
		            else val = targetAvg - v[ix];
		if (val >= 0.0)
			{ allMonotonic = false;  break; }
		}

	if (allMonotonic)
		{
		if (op->debug)
			{
			if (aboveThresh) fprintf (stderr, " all decreasing\n");
			            else fprintf (stderr, " all increasing\n");
			}

		for (ix=0 ; ix<vLen ; ix++)
			v[ix] = zeroVal;
		goto release_memory;
		}

	// search for clumps, intervals for which the average is above (or at) the
	// threshold;  this is equivalent to intervals in which the sum, minus the
	// length times the average, is positive or zero
	//
	// nota bene: this algorithm can suffer from round off error in the sums,
	//            which can make the reported intervals less precise than we'd
	//            like

	minSum = valSum = 0.0;

	minSums [0] = valSum;
	minWhere[0] = (u32) -1;
	numMinSums  = 1;
	minScan     = 0;

	prevStart   = (u32) -1;
	prevEnd     = (u32) -1;

	for (ix=0 ; ix<vLen ; ix++)
		{
		if ((op->progress != 0)
		 && (ix % op->progress == 0))
			fprintf (stderr, " progress %s %s/%s (%.1f%%)\n",
			                 vName, ucommatize(ix), ucommatize(vLen),
			                 (100.0 * ix)/vLen);

		// invariant: minSums[minScan] <= valSum <  minSums[minScan-1]
		// with the implied assumption that minSums[-1] == infinity

		if (aboveThresh) val = v[ix] - targetAvg;
		            else val = targetAvg - v[ix];
		s[ix] = zeroVal;

		valSum += val;
		if (valSum < minSum)
			{
			minSum = valSum;
			minSums [numMinSums] = valSum;
			minWhere[numMinSums] = ix;
			numMinSums++;
			}

		// (re-establish the invariant)
		if (val < 0)								// the sum has decreased
			{ while (minSums[minScan] > valSum) minScan++; }
		else if (val > 0)							// the sum has increased
			{ while ((minScan > 0) && (minSums[minScan-1] <= valSum)) minScan--; }

		if (op->debugDetail)
			fprintf (stderr, " [%u] " valtypeFmt " " valtypeFmt " %u..\n",
			                 ix,val,valSum,minWhere[minScan]);

		// minScan points at the earliest index with minSums[minScan] <= valSum,
		// so the interval (minWhere[minScan]+1 to ix) has sum >= 0

		minIx = minWhere[minScan];
		if (ix - minIx < minLength) continue;

		start = minIx+1;
		end   = ix;

		if (op->debug)
			fprintf (stderr, " setting %u..%u\n", start, end+1);

		if ((prevStart == (u32) -1)
		 || (start > prevEnd+1))
			{
			// no overlap with previous interval (or no previous interval)
			for (iy=start ; iy<=end ; iy++) s[iy] = oneVal;
			prevStart = start;
			prevEnd   = end;
			}
		else if (start >= prevStart)
			{
			// new interval extends previous interval on right
			for (iy=prevEnd+1 ; iy<=end ; iy++) s[iy] = oneVal;
			prevEnd = end;
			}
		else
			{
			// new interval extends previous interval on left and right
			for (iy=start     ; iy<prevStart ; iy++) s[iy] = oneVal;
			for (iy=prevEnd+1 ; iy<=end      ; iy++) s[iy] = oneVal;
			prevStart = start;
			prevEnd   = end;
			}

		}

	// copy scratch array to vector, trimming any below-average ends off of
	// the intervals

	scanIx = 0;
	while (true)
		{
		// find the start of the next non-zero scratch interval, clearing the
		// vector as we go

		for (start=scanIx ; start<vLen ; start++)
			{
			if (s[start] != zeroVal) break;
			v[start] = zeroVal;
			}
		if (start >= vLen) break;

		// scan over the first part of the interval, clearing any below-average
		// values in the vector

		if (op->debug)
			{
			if (aboveThresh) fprintf (stderr, " clump start %u\n", start);
			            else fprintf (stderr, " skimp start %u\n", start);
			}

		for (ix=start ; ix<vLen ; ix++)
			{
			if (s[ix] == zeroVal) break;
			if (aboveThresh) { if (v[ix] >= targetAvg) break; }
			            else { if (v[ix] <= targetAvg) break; }
			v[ix] = zeroVal;
			}

		if ((ix >= vLen) || (s[ix] == zeroVal))
			{ scanIx = ix;  continue; }

		// scan over the middle and end parts of the interval, locating the
		// last average-or-better value in the vector

		if (op->debug)
			fprintf (stderr, "   start trimmed to %u\n", ix);

		start = ix;
		end   = ix++;
		for ( ; ix<vLen ; ix++)
			{
			if (s[ix] == zeroVal) break;
			if (aboveThresh) { if (v[ix] >= targetAvg) end = ix; }
			            else { if (v[ix] <= targetAvg) end = ix; }
			}

		if (op->debug)
			{
			fprintf (stderr, "   end %u\n", ix);
			fprintf (stderr, "   end trimmed to %u\n", end+1);
			}

		for (iy=start ; iy<=end ; iy++)
			v[iy] = oneVal;
		for (iy=end+1 ; iy<ix ; iy++)
			v[iy] = zeroVal;

		scanIx = ix;
		}

release_memory:
	release_scratch_vector(s);
	release_scratch_vector(minSums);
	release_scratch_ints  ((s32*) minWhere);

	// success

	return;

	// failure

no_threshold:
	fprintf (stderr, "[%s] attempt to use %s as threshold failed (no such variable)\n",
	                 _op->name, op->averageVarName);
	exit(EXIT_FAILURE);
	}

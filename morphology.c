// morphology.c-- genodsp operators performing Minkowkski morphology operations

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
#include "morphology.h"

//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_close--
//	Perform a "close" operation, in the sense of Minkowski/morphology set
//	operations on covered intervals.
//
//	The incoming signal is first binarized (conceptually), and the result is
//	1 or 0 depending depending on whether a position is within the closure or
//	not.
//
//----------

// private dspop subtype

typedef struct dspop_close
	{
	dspop		common;			// common elements shared with all operators
	u32			closingLength;
	int			haveThreshold;
	char*		thresholdVarName;
	valtype		threshold;
	valtype		oneVal;
	valtype		zeroVal;
	int			debug;
	} dspop_close;


// op_close_short--

void op_close_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

 	fprintf (f, "apply interval closure (fill small gaps between intervals)\n");
	}


// op_close_usage--

void op_close_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sCompute the morphological closure of the current set of intervals. This is\n",    indent);
	fprintf (f, "%sequivalent to \"filling\" short gaps between intervals.\n",                       indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe current signal is first binarized according to the given threshold; any\n",   indent);
	fprintf (f, "%slocations above the threshold are considered to be \"in\" the set. Closure is\n", indent);
	fprintf (f, "%sthen performed, adding any short gaps into the set. The result is one for\n",     indent);
	fprintf (f, "%slocations in the resulting zet, zero otherwise. Note that short gaps at the\n",   indent);
	fprintf (f, "%send of chromosomes are NOT filled.\n",                                            indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sStrictly speaking, the length parameter is twice the true closure length. For\n", indent);
	fprintf (f, "%smore information on set morphology, see\n",                                       indent);
	fprintf (f, "%s  en.wikipedia.org/wiki/Closing_(morphology).\n",                                 indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <length> [options]\n", indent, name);
	fprintf (f, "%s  <length>                 the \"length\" of closure; gaps of this length\n",     indent);
	fprintf (f, "%s                           or shorter will be filled\n",                          indent);
	fprintf (f, "%s  --threshold=<value>      (T=) a binarization threshold; values above this\n",   indent);
	fprintf (f, "%s                           are considered to be \"in\" the incoming intervals\n", indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                    indent);
	fprintf (f, "%s  --one=<value>            (O=) value for locations in the closure\n",            indent);
	fprintf (f, "%s                           (default is 1.0)\n",                                   indent);
	fprintf (f, "%s  --zero=<value>           (Z=) value for locations not in the closure\n",        indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	}

// op_close_parse--

dspop* op_close_parse (char* name, int _argc, char** _argv)
	{
	dspop_close*	op;
	int				argc = _argc;
	char**			argv = _argv;
	char*			arg, *argVal;
	valtype			tempVal;
	int				haveLength;

	// allocate and initialize our control record

	op = (dspop_close*) malloc (sizeof(dspop_close));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->closingLength    = 0;      // not used, user is required to set it
	op->haveThreshold    = false;
	op->thresholdVarName = NULL;
	op->threshold        = 0.0;
	op->oneVal           = 1.0;
	op->zeroVal          = 0.0;
	op->debug            = false;

	// parse arguments

	haveLength = false;

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
			if (op->haveThreshold) goto more_than_one_threshold;
			if (!try_string_to_valtype (argVal, &op->threshold))
				op->thresholdVarName = copy_string (argVal);
			else
				op->threshold = string_to_valtype (argVal);
			op->haveThreshold = true;
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

		// --debug argument

		if (strcmp (arg, "--debug") == 0)
			{ op->debug = true;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <length>

		if (!haveLength)
			{
			op->closingLength = (u32) string_to_unitized_int (arg, /*thousands*/ true);
			haveLength = true;
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if (!haveLength) goto length_missing;

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_close));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_threshold:
	fprintf (stderr, "[%s] threshold specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

length_missing:
	fprintf (stderr, "[%s] closing length was not provided\n", name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_close_free--

void op_close_free (dspop* _op)
	{
	dspop_close*	op = (dspop_close*) _op;

	if (op->thresholdVarName != NULL) free (op->thresholdVarName);
	free (op);
	}


// op_close_apply--

void op_close_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_close*	op = (dspop_close*) _op;
	valtype			closingLength = op->closingLength;
	valtype			cutoffThresh  = op->threshold;
	valtype			oneVal        = op->oneVal;
	valtype			zeroVal       = op->zeroVal;
	int				ok;
	int				inGap;
	u32				gapStartIx = 0; // placate compiler
	u32				ix, scanIx;

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

	// process the vector, filling short gaps (and binarizing);  note that we
	// never fill the first gap because its "true" length is infinite (it
	// includes everything on the number line before the start of the
	// chromosome)

	if (v[0] <= cutoffThresh) { gapStartIx = 0;  inGap = true; }
	                     else { v[0] = oneVal;   inGap = false;  }

	for (ix=1 ; ix<vLen ; ix++)
		{
		if (v[ix] <= cutoffThresh)			// new location is in a gap
			{
			if (inGap) continue;
			gapStartIx = ix;
			inGap = true;
			}
		else								// new location is in an interval
			{
			v[ix] = oneVal;
			if (!inGap) continue;
			if ((gapStartIx == 0) || (ix - gapStartIx > closingLength))
				{
				if (op->debug) fprintf (stderr, "clear %u bp gap from %u..%u\n",
				                                ix-gapStartIx, gapStartIx, ix);
				for (scanIx=gapStartIx ; scanIx<ix ; scanIx++) v[scanIx] = zeroVal;
				}
			else
				{
				if (op->debug) fprintf (stderr, "fill %u bp gap from %u..%u\n",
				                                ix-gapStartIx, gapStartIx, ix);
				for (scanIx=gapStartIx ; scanIx<ix ; scanIx++) v[scanIx] = oneVal;
				}
			inGap = false;
			}
		}

	// clear the final gap;  we never fill it because its "true" length is
	// infinite (it includes everything on the number line beyond the end of
	// the chromosome)

	if (inGap)
		{
		if (op->debug) fprintf (stderr, "clear %u bp gap from %u..%u\n",
				                        vLen-gapStartIx, gapStartIx, vLen);
		for (scanIx=gapStartIx ; scanIx<vLen ; scanIx++) v[scanIx] = zeroVal;
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
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_open--
//	Perform a "open" operation, in the sense of Minkowski/morphology set
//	operations on covered intervals.
//
//	The incoming signal is first binarized (conceptually), and the result is
//	1 or 0 depending depending on whether a position is within the opened set
//	or not.
//
//----------

// private dspop subtype

typedef struct dspop_open
	{
	dspop		common;			// common elements shared with all operators
	u32			openingLength;
	int			haveThreshold;
	char*		thresholdVarName;
	valtype		threshold;
	valtype		oneVal;
	valtype		zeroVal;
	} dspop_open;


// op_open_short--

void op_open_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "apply interval opening (remove small intervals)\n");
	}


// op_open_usage--

void op_open_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sCompute the morphological opening of the current set of intervals. This is\n",    indent);
	fprintf (f, "%sequivalent to \"removing\" short intervals.\n",                                   indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe current signal is first binarized according to the given threshold; any\n",   indent);
	fprintf (f, "%slocations above the threshold are considered to be \"in\" the set. Opening is\n", indent);
	fprintf (f, "%sthen performed, removing any short intervals from the set. The result is one\n",  indent);
	fprintf (f, "%sfor locations in the resulting zet, zero otherwise.\n",                           indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sStrictly speaking, the length parameter is twice the true opening length. For\n", indent);
	fprintf (f, "%smore information on set morphology, see\n",                                       indent);
	fprintf (f, "%s  en.wikipedia.org/wiki/Opening_(morphology).\n",                                 indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <length> [options]\n", indent, name);
	fprintf (f, "%s  <length>                 the \"length\" of opening; intervals of this\n",       indent);
	fprintf (f, "%s                           length or shorter will be removed\n",                  indent);
	fprintf (f, "%s  --threshold=<value>      (T=) a binarization threshold; values above this\n",   indent);
	fprintf (f, "%s                           are considered to be \"in\" the incoming intervals\n", indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                    indent);
	fprintf (f, "%s  --one=<value>            (O=) value for locations in the opened set\n",         indent);
	fprintf (f, "%s                           (default is 1.0)\n",                                   indent);
	fprintf (f, "%s  --zero=<value>           (Z=) value for locations not in the opened set\n",     indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	}

// op_open_parse--

dspop* op_open_parse (char* name, int _argc, char** _argv)
	{
	dspop_open*	op;
	int				argc = _argc;
	char**			argv = _argv;
	char*			arg, *argVal;
	valtype			tempVal;
	int				haveLength;

	// allocate and initialize our control record

	op = (dspop_open*) malloc (sizeof(dspop_open));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->openingLength    = 0;      // not used, user is required to set it
	op->haveThreshold    = false;
	op->thresholdVarName = NULL;
	op->threshold        = 0.0;
	op->oneVal           = 1.0;
	op->zeroVal          = 0.0;

	// parse arguments

	haveLength = false;

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
			if (op->haveThreshold) goto more_than_one_threshold;
			if (!try_string_to_valtype (argVal, &op->threshold))
				op->thresholdVarName = copy_string (argVal);
			else
				op->threshold = string_to_valtype (argVal);
			op->haveThreshold = true;
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

		// <length>

		if (!haveLength)
			{
			op->openingLength = (u32) string_to_unitized_int (arg, /*thousands*/ true);
			haveLength = true;
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if (!haveLength) goto length_missing;

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_open));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_threshold:
	fprintf (stderr, "[%s] threshold specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

length_missing:
	fprintf (stderr, "[%s] opening length was not provided\n", name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_open_free--

void op_open_free (dspop* _op)
	{
	dspop_open*	op = (dspop_open*) _op;

	if (op->thresholdVarName != NULL) free (op->thresholdVarName);
	free (op);
	}


// op_open_apply--

void op_open_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_open*	op = (dspop_open*) _op;
	valtype			openingLength = op->openingLength;
	valtype			cutoffThresh  = op->threshold;
	valtype			oneVal        = op->oneVal;
	valtype			zeroVal       = op->zeroVal;
	int				ok;
	int				inInterval;
	u32				intervalStartIx = 0; // placate compiler
	u32				ix, scanIx;

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

	// process the vector, clearing short intervals (and binarizing)

	if (v[0] > cutoffThresh) { intervalStartIx = 0;  inInterval = true; }
	                    else { v[0] = zeroVal;       inInterval = false;  }

	for (ix=1 ; ix<vLen ; ix++)
		{
		if (v[ix] > cutoffThresh)			// new location is in an interval
			{
			if (inInterval) continue;
			intervalStartIx = ix;
			inInterval      = true;
			}
		else								// new location is in a gap
			{
			v[ix] = zeroVal;
			if (!inInterval) continue;
			if (ix - intervalStartIx > openingLength)
				{ for (scanIx=intervalStartIx ; scanIx<ix ; scanIx++) v[scanIx] = oneVal; }
			else
				{ for (scanIx=intervalStartIx ; scanIx<ix ; scanIx++) v[scanIx] = zeroVal; }
			inInterval = false;
			}
		}

	// fill or clear the final interval

	if (inInterval)
		{
		if (vLen - intervalStartIx > openingLength)
			{ for (scanIx=intervalStartIx ; scanIx<vLen ; scanIx++) v[scanIx] = oneVal; }
		else
			{ for (scanIx=intervalStartIx ; scanIx<vLen ; scanIx++) v[scanIx] = zeroVal; }
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
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_dilate--
//	Perform a "dilation" operation, in the sense of Minkowski/morphology set
//	operations on covered intervals.
//
//	The incoming signal is first binarized (conceptually), and the result is
//	1 or 0 depending depending on whether a position is within the dilation set
//	or not.
//
//----------

// private dspop subtype

typedef struct dspop_dilate
	{
	dspop		common;			// common elements shared with all operators
	u32			dilationLength;
	u32			leftDilation;
	u32			rightDilation;
	int			haveThreshold;
	char*		thresholdVarName;
	valtype		threshold;
	valtype		oneVal;
	valtype		zeroVal;
	int			debug;
	} dspop_dilate;


// op_dilate_short--

void op_dilate_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "apply interval dilation (widen intervals)\n");
	}


// op_dilate_usage--

void op_dilate_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sCompute the morphological dilation of the current set of intervals. This is\n",    indent);
	fprintf (f, "%sequivalent to widening intervals.\n",                                              indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe current signal is first binarized according to the given threshold; any\n",    indent);
	fprintf (f, "%slocations above the threshold are considered to be \"in\" the set. Dilation is\n", indent);
	fprintf (f, "%sthen performed, widening any intervals in the set. The result is one for\n",       indent);
	fprintf (f, "%slocations in the resulting zet, zero otherwise.\n",                                indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sWe assume that the signal is equivalent to a zero beyond the left and right\n",    indent);
	fprintf (f, "%sends of the vector.\n",                                                            indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sStrictly speaking, the length parameter is twice the true dilation length.\n",     indent);
	fprintf (f, "%sFor more information on set morphology, see\n",                                    indent);
	fprintf (f, "%s  en.wikipedia.org/wiki/Dilation_(morphology).\n",                                 indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <length> [options]\n", indent, name);
	fprintf (f, "%s  <length>                 the \"length\" of dilation; intervals will be\n",       indent);
	fprintf (f, "%s                           widened by this amount, half on each side\n",           indent);
	fprintf (f, "%s  --threshold=<value>      (T=) a binarization threshold; values above this\n",    indent);
	fprintf (f, "%s                           are considered to be \"in\" the incoming intervals\n",  indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                    indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                     indent);
	fprintf (f, "%s  --one=<value>            (O=) value for locations in the dilation set\n",        indent);
	fprintf (f, "%s                           (default is 1.0)\n",                                    indent);
	fprintf (f, "%s  --zero=<value>           (Z=) value for locations not in the dilation set\n",    indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                    indent);
	fprintf (f, "%s  --left=<length>          intervals will be widened by this amount, only on\n",   indent);
	fprintf (f, "%s                           the left side\n",                                       indent);
	fprintf (f, "%s  --right=<length>         intervals will be widened by this amount, only on\n",   indent);
	fprintf (f, "%s                           the right side\n",                                      indent);
	}

// op_dilate_parse--

dspop* op_dilate_parse (char* name, int _argc, char** _argv)
	{
	dspop_dilate*	op;
	int				argc = _argc;
	char**			argv = _argv;
	char*			arg, *argVal;
	valtype			tempVal;
	int				haveLength;

	// allocate and initialize our control record

	op = (dspop_dilate*) malloc (sizeof(dspop_dilate));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->dilationLength   = 0;      // not used, user is required to set it
	op->leftDilation     = 0;
	op->rightDilation    = 0;
	op->haveThreshold    = false;
	op->thresholdVarName = NULL;
	op->threshold        = 0.0;
	op->oneVal           = 1.0;
	op->zeroVal          = 0.0;
	op->debug            = false;

	// parse arguments

	haveLength = false;

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
			if (op->haveThreshold) goto more_than_one_threshold;
			if (!try_string_to_valtype (argVal, &op->threshold))
				op->thresholdVarName = copy_string (argVal);
			else
				op->threshold = string_to_valtype (argVal);
			op->haveThreshold = true;
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

		// --left=<length> or --right=<length>

		if (strcmp_prefix (arg, "--left=") == 0)
			{
			tempVal = string_to_valtype (argVal);
			op->leftDilation = tempVal;
			goto next_arg;
			}

		if (strcmp_prefix (arg, "--right=") == 0)
			{
			tempVal = string_to_valtype (argVal);
			op->rightDilation = tempVal;
			goto next_arg;
			}

		// --debug argument

		if (strcmp (arg, "--debug") == 0)
			{ op->debug = true;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <length>

		if (!haveLength)
			{
			op->dilationLength = (u32) string_to_unitized_int (arg, /*thousands*/ true);
			haveLength = true;
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	// user has to specify the length, OR specify a left-only or right-only
	// length

	if (haveLength)
		{
		if ((op->leftDilation != 0) || (op->rightDilation != 0)) goto more_than_one_length;
		}
	else // if (!haveLength)
		{
		if ((op->leftDilation == 0) && (op->rightDilation == 0)) goto length_missing;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_dilate));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_threshold:
	fprintf (stderr, "[%s] threshold specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

length_missing:
	fprintf (stderr, "[%s] dilation length was not provided\n", name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_length:
	fprintf (stderr, "[%s] dilation length was provided in more than one way\n", name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_dilate_free--

void op_dilate_free (dspop* _op)
	{
	dspop_dilate*	op = (dspop_dilate*) _op;

	if (op->thresholdVarName != NULL) free (op->thresholdVarName);
	free (op);
	}


// op_dilate_apply--

void op_dilate_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_dilate*	op = (dspop_dilate*) _op;
	valtype			dilationLength = op->dilationLength;
	valtype			cutoffThresh   = op->threshold;
	valtype			oneVal         = op->oneVal;
	valtype			zeroVal        = op->zeroVal;
	int				ok;
	int				inInterval;
	u32				leftDilation, rightDilation;
	u32				gapStartIx, leftIx, rightIx, ix, scanIx;

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

	// process the vector, widening intervals (and binarizing);  note that
	// locations in intervals are binarized immediately;  locations in gaps
	// are not binarized until we have determined the start and end of the gap

	if ((op->leftDilation == 0) && (op->rightDilation == 0))
		{
		leftDilation  = dilationLength / 2;
		rightDilation = dilationLength - leftDilation;
		}
	else
		{
		leftDilation  = op->leftDilation;
		rightDilation = op->rightDilation;
		}

	gapStartIx = 0;
	inInterval = (v[0] > cutoffThresh);
	if (inInterval) v[0] = oneVal;

	for (ix=1 ; ix<vLen ; ix++)
		{
		if (v[ix] <= cutoffThresh)			// new location is in a gap
			{
			if (!inInterval) continue;

			// we've encountered the start of a new gap;  record the start

			gapStartIx = ix;
			inInterval = false;
			}
		else								// new location is in an interval
			{
			v[ix] = oneVal;
			if (inInterval) continue;

			// we've encountered the start of a new interval;  set and/or clear
			// items to narrow the gap between this interval and the last

			if (op->debug) fprintf (stderr, "%u bp gap from %u..%u\n",
			                                ix-gapStartIx, gapStartIx, ix);

			if (gapStartIx == 0)
				{
				// gap is at left edge, so there's no fill on the gap's left

				if (ix <= leftDilation)
					{
					// left-extension of new interval reaches left edge
					if (op->debug)
						{
						fprintf (stderr, "  narrowed gap is -%u..%u\n", leftDilation-ix, 0);
						fprintf (stderr, "  set   %u..%u\n", 0, ix);
						}
					for (scanIx=0 ; scanIx<ix ; scanIx++) v[scanIx] = oneVal;
					}
				else
					{
					// left-extension of new interval leaves a gap at left edge
					leftIx = ix - leftDilation;

					if (op->debug)
						{
						fprintf (stderr, "  narrowed gap is %u..%u\n", 0,      leftIx);
						fprintf (stderr, "  clear %u..%u\n",           0,      leftIx);
						fprintf (stderr, "  set   %u..%u\n",           leftIx, ix);
						}
					for (scanIx=0 ; scanIx<leftIx ; scanIx++) v[scanIx] = zeroVal;
					for (         ; scanIx<ix     ; scanIx++) v[scanIx] = oneVal;
					}
				}
			else
				{
				rightIx = gapStartIx + rightDilation;
				leftIx  = (ix <= leftDilation)? 0 : ix - leftDilation;

				if (op->debug) fprintf (stderr, "  narrowed gap is %u..%u\n", rightIx, leftIx);

				if (rightIx >= leftIx)
					{
					// right-extension of previous interval and left-extension
					// of new interval overlap (or abut), so gap disappears

					if (op->debug)
						fprintf (stderr, "  set   %u..%u\n", gapStartIx, ix);
					for (scanIx=gapStartIx ; scanIx<ix ; scanIx++) v[scanIx] = oneVal;
					}
				else
					{
					// right-extension of previous interval and left-extension
					// of new interval leave a narrowed gap between intervals

					if (op->debug)
						{
						fprintf (stderr, "  set   %u..%u\n", gapStartIx, rightIx);
						fprintf (stderr, "  clear %u..%u\n", rightIx,    leftIx);
						fprintf (stderr, "  set   %u..%u\n", leftIx,     ix);
						}
					for (scanIx=gapStartIx ; scanIx<rightIx ; scanIx++) v[scanIx] = oneVal;
					for (                  ; scanIx<leftIx  ; scanIx++) v[scanIx] = zeroVal;
					for (                  ; scanIx<ix      ; scanIx++) v[scanIx] = oneVal;
					}
				}

			inInterval = true;
			}
		}

	// fill the final gap

	if (!inInterval)
		{
		rightIx = gapStartIx + rightDilation;

		if (op->debug)
			{
			fprintf (stderr, "%u bp gap from %u..%u\n", vLen-gapStartIx, gapStartIx, vLen);
			fprintf (stderr, "  narrowed gap is %u..%u\n", rightIx, vLen);
			}

		if (gapStartIx == 0)
			{
			// gap is entire vector, so there's no fill on either edge

			if (op->debug)
				fprintf (stderr, "  clear %u..%u\n", 0, vLen);
			for (scanIx=0 ; scanIx<vLen ; scanIx++) v[scanIx] = zeroVal;
			}
		else if (rightIx >= vLen)
			{
			// right-extension of previous interval reaches right edge

			if (op->debug)
				fprintf (stderr, "  set   %u..%u\n", gapStartIx, vLen);
			for (scanIx=gapStartIx ; scanIx<vLen ; scanIx++) v[scanIx] = oneVal;
			}
		else
			{
			// right-extension of previous interval leaves a gap at right edge

			if (op->debug)
				{
				fprintf (stderr, "  set   %u..%u\n", gapStartIx, rightIx);
				fprintf (stderr, "  clear %u..%u\n", rightIx,    vLen);
				}
			for (scanIx=gapStartIx ; scanIx<rightIx ; scanIx++) v[scanIx] = oneVal;
			for (                  ; scanIx<vLen    ; scanIx++) v[scanIx] = zeroVal;
			}
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
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_erode--
//	Perform a "erosion" operation, in the sense of Minkowski/morphology set
//	operations on covered intervals.
//
//	The incoming signal is first binarized (conceptually), and the result is
//	1 or 0 depending depending on whether a position is within the erosion set
//	or not.
//
//----------

// private dspop subtype

typedef struct dspop_erode
	{
	dspop		common;			// common elements shared with all operators
	u32			erosionLength;
	u32			leftErosion;
	u32			rightErosion;
	int			haveThreshold;
	char*		thresholdVarName;
	valtype		threshold;
	valtype		oneVal;
	valtype		zeroVal;
	int			debug;
	} dspop_erode;


// op_erode_short--

void op_erode_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "apply interval erosion (shrink intervals)\n");
	}


// op_erode_usage--

void op_erode_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sCompute the morphological erosion of the current set of intervals. This is\n",     indent);
	fprintf (f, "%sequivalent to shrinking intervals.\n",                                             indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe current signal is first binarized according to the given threshold; any\n",    indent);
	fprintf (f, "%slocations above the threshold are considered to be \"in\" the set. Erosion is\n",  indent);
	fprintf (f, "%sthen performed, shrinking any intervals in the set. The result is one for\n",      indent);
	fprintf (f, "%slocations in the resulting zet, zero otherwise.\n",                                indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sWe assume that the signal is equivalent to a zero beyond the left and right\n",    indent);
	fprintf (f, "%sends of the vector.\n",                                                            indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sStrictly speaking, the length parameter is twice the true erosion length.\n",      indent);
	fprintf (f, "%sFor more information on set morphology, see\n",                                    indent);
	fprintf (f, "%s  en.wikipedia.org/wiki/Erosion_(morphology).\n",                                  indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <length> [options]\n", indent, name);
	fprintf (f, "%s  <length>                 the \"length\" of erosion; intervals will be\n",        indent);
	fprintf (f, "%s                           shrunk by this amount, half on each side\n",            indent);
	fprintf (f, "%s  --threshold=<value>      (T=) a binarization threshold; values above this\n",    indent);
	fprintf (f, "%s                           are considered to be \"in\" the incoming intervals\n",  indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                    indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                     indent);
	fprintf (f, "%s  --one=<value>            (O=) value for locations in the erosion set\n",         indent);
	fprintf (f, "%s                           (default is 1.0)\n",                                    indent);
	fprintf (f, "%s  --zero=<value>           (Z=) value for locations not in the erosion set\n",     indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                    indent);
	fprintf (f, "%s  --left=<length>          intervals will be shrunk by this amount, only on\n",    indent);
	fprintf (f, "%s                           the left side\n",                                       indent);
	fprintf (f, "%s  --right=<length>         intervals will be shrunk by this amount, only on\n",    indent);
	fprintf (f, "%s                           the right side\n",                                      indent);
	}

// op_erode_parse--

dspop* op_erode_parse (char* name, int _argc, char** _argv)
	{
	dspop_erode*	op;
	int				argc = _argc;
	char**			argv = _argv;
	char*			arg, *argVal;
	valtype			tempVal;
	int				haveLength;

	// allocate and initialize our control record

	op = (dspop_erode*) malloc (sizeof(dspop_erode));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->erosionLength    = 0;      // not used, user is required to set it
	op->leftErosion      = 0;
	op->rightErosion     = 0;
	op->haveThreshold    = false;
	op->thresholdVarName = NULL;
	op->threshold        = 0.0;
	op->oneVal           = 1.0;
	op->zeroVal          = 0.0;
	op->debug            = false;

	// parse arguments

	haveLength = false;

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
			if (op->haveThreshold) goto more_than_one_threshold;
			if (!try_string_to_valtype (argVal, &op->threshold))
				op->thresholdVarName = copy_string (argVal);
			else
				op->threshold = string_to_valtype (argVal);
			op->haveThreshold = true;
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

		// --left=<length> or --right=<length>

		if (strcmp_prefix (arg, "--left=") == 0)
			{
			tempVal = string_to_valtype (argVal);
			op->leftErosion = tempVal;
			goto next_arg;
			}

		if (strcmp_prefix (arg, "--right=") == 0)
			{
			tempVal = string_to_valtype (argVal);
			op->rightErosion = tempVal;
			goto next_arg;
			}

		// --debug argument

		if (strcmp (arg, "--debug") == 0)
			{ op->debug = true;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <length>

		if (!haveLength)
			{
			op->erosionLength = (u32) string_to_unitized_int (arg, /*thousands*/ true);
			haveLength = true;
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	// user has to specify the length, OR specify a left-only or right-only
	// length

	if (haveLength)
		{
		if ((op->leftErosion != 0) || (op->rightErosion != 0)) goto more_than_one_length;
		}
	else // if (!haveLength)
		{
		if ((op->leftErosion == 0) && (op->rightErosion == 0)) goto length_missing;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_erode));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_threshold:
	fprintf (stderr, "[%s] threshold specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

length_missing:
	fprintf (stderr, "[%s] erosion length was not provided\n", name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_length:
	fprintf (stderr, "[%s] erosion length was provided in more than one way\n", name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_erode_free--

void op_erode_free (dspop* _op)
	{
	dspop_erode*	op = (dspop_erode*) _op;

	if (op->thresholdVarName != NULL) free (op->thresholdVarName);
	free (op);
	}


// op_erode_apply--

void op_erode_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_erode*	op = (dspop_erode*) _op;
	valtype			erosionLength = op->erosionLength;
	valtype			cutoffThresh  = op->threshold;
	valtype			oneVal        = op->oneVal;
	valtype			zeroVal       = op->zeroVal;
	int				ok;
	int				inInterval;
	u32				leftErosion, rightErosion;
	u32				startIx, leftIx, rightIx, ix, scanIx;

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

	// process the vector, widening intervals (and binarizing);  note that
	// locations in gaps are binarized immediately;  locations in intervals
	// are not binarized until we have determined the start and end of the
	// interval;  also note that the loop takes one an extra pass, to include
	// the location *beyond* the right end, which simplifies the handling of
	// the final interval

	if ((op->leftErosion == 0) && (op->rightErosion == 0))
		{
		leftErosion  = erosionLength / 2;
		rightErosion = erosionLength - leftErosion;
		}
	else
		{
		leftErosion  = op->leftErosion;
		rightErosion = op->rightErosion;
		}

	startIx = 0;
	inInterval = (v[0] > cutoffThresh);
	if (!inInterval) v[0] = zeroVal;

	for (ix=1 ; ix<=vLen ; ix++)			// (nota bene: loop steps beyond right end)
		{
		if (ix == vLen) goto in_gap;
		
		if (v[ix] > cutoffThresh)			// new location is in an interval
			{
			if (inInterval) continue;

			// we've encountered the start of a new interval;  record the start

			startIx = ix;
			inInterval = true;
			}
		else								// new location is in a gap
			{
			v[ix] = zeroVal;
		in_gap:
			if (!inInterval) continue;

			// we've encountered the start of a new gap;  set and/or clear
			// items to narrow the interval between this gap and the last

			rightIx = startIx + rightErosion;
			leftIx  = ix - leftErosion;

			if (op->debug)
				{
				fprintf (stderr, "%u bp interval from %u..%u\n", ix-startIx, startIx, ix);
				fprintf (stderr, "  narrowed interval is %u..%u\n", rightIx, leftIx);
				}

			if (rightIx >= leftIx)
				{
				// right-extension of previous interval and left-extension
				// of new interval overlap (or abut), so interval disappears

				if (op->debug)
					fprintf (stderr, "  clear %u..%u\n", startIx, ix);
				for (scanIx=startIx ; scanIx<ix ; scanIx++) v[scanIx] = zeroVal;
				}
			else
				{
				// right-extension of previous interval and left-extension
				// of new interval leave a narrowed interval between gaps

				if (op->debug)
					{
					fprintf (stderr, "  clear %u..%u\n", startIx, rightIx);
					fprintf (stderr, "  set   %u..%u\n", rightIx, leftIx);
					fprintf (stderr, "  clear %u..%u\n", leftIx,  ix);
					}
				for (scanIx=startIx ; scanIx<rightIx ; scanIx++) v[scanIx] = zeroVal;
				for (               ; scanIx<leftIx  ; scanIx++) v[scanIx] = oneVal;
				for (               ; scanIx<ix      ; scanIx++) v[scanIx] = zeroVal;
				}

			inInterval = false;
			}
		}

	// success

	return;

	// failure

no_threshold:
	fprintf (stderr, "[%s] attempt to use %s as threshold failed (no such variable)\n",
	                 _op->name, op->thresholdVarName);
	exit(EXIT_FAILURE);
	}


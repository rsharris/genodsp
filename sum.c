// sum.c-- genodsp operators performing windowed sum

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
#include "sum.h"

// private dspop subtype, shared by op_window_sum and op_sliding_sum--

typedef struct dspop_sum
	{
	dspop		common;			// common elements shared with all operators
	u32			windowSize;
	valtype		denominator;
	valtype		zeroVal;
	int			windowIsChromosome;	// (only for op_window_sum)
	int			useActualDenom;	    // (only for op_window_sum)
	int			denomIsWindowSize;  // (only for op_window_sum)
	} dspop_sum;

//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_window_sum--
//	Apply a windowed sum;  the signal is replaced by its sum over non-
//	overlapping windows.  The window's sum is written to the first entry in the
//	window, with the remaining entries set to zero.
//
//----------

// op_window_sum_short--

void op_window_sum_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "sum over non-overlapping windows\n");
	}


// op_window_sum_usage--

void op_window_sum_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sApply a windowed sum;  the signal is replaced by its sum over non-overlapping\n", indent);
	fprintf (f, "%swindows. The window's sum is written to the first entry in the window, with\n",   indent);
	fprintf (f, "%sthe remaining entries set to zero.\n",                                            indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --window=chromosome      entire chromosome is a single window\n",               indent);
	fprintf (f, "%s  --window=<length>        (W=) size of window\n",                                indent);
	fprintf (f, "%s  --denom=<value>          (D=) denominator; \"window\" or \"W\" can be used\n",  indent);
	fprintf (f, "%s                           as the value; the value \"actual\" means to use the\n",indent);
	fprintf (f, "%s                           actual size of the window, in case it's truncated\n",  indent);
	fprintf (f, "%s                           (by default, we do not use a denominator)\n",          indent);
	fprintf (f, "%s  --zero=<value>           (Z=) zero value to fill window body\n",                indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	}


// op_window_sum_parse--

dspop* op_window_sum_parse (char* name, int _argc, char** _argv)
	{
	dspop_sum*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_sum*) malloc (sizeof(dspop_sum));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->windowSize         = (u32) get_named_global ("windowSize", 100);
	op->denominator        = 1.0;
	op->zeroVal            = 0.0;
	op->windowIsChromosome = false;
	op->useActualDenom     = false;
	op->denomIsWindowSize  = false;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// -window=chromosome, --window=<length> or W=<length>

		if (strcmp (arg, "--window=chromosome") == 0)
			{ op->windowIsChromosome = true;  goto next_arg; }

		if ((strcmp_prefix (arg, "--window=") == 0)
		 || (strcmp_prefix (arg, "W=")        == 0)
		 || (strcmp_prefix (arg, "--W=")      == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt == 0)
				chastise ("[%s] window size can't be zero (\"%s\")\n", name, arg);
			if (tempInt < 0)
				chastise ("[%s] window size can't be negative (\"%s\")\n", name, arg);
			if (tempInt < 3)
				{
				fprintf (stderr, "[%s] WARNING: raising window size from %d to %d\n",
				                 name, tempInt, 3);
				tempInt = 3;
				}
			op->windowSize = (u32) tempInt;
			op->windowIsChromosome = false;
			goto next_arg;
			}

		// --denom[inator]=<value> or D=<value>

		if ((strcmp_prefix (arg, "--denom=")       == 0)
		 || (strcmp_prefix (arg, "--denominator=") == 0)
		 || (strcmp_prefix (arg, "D=")             == 0)
		 || (strcmp_prefix (arg, "--D=")           == 0))
			{
			op->denominator       = 1.0;
			op->useActualDenom    = false;
			op->denomIsWindowSize = false;
			if (strcmp (argVal, "actual") == 0)
				{ op->useActualDenom = true;  goto next_arg; }
			if ((strcmp (argVal, "window") == 0)
			 || (strcmp (argVal, "W")      == 0))
				{ op->denomIsWindowSize = true;  goto next_arg; }
			tempVal = string_to_valtype (argVal);
			if (tempVal == 0)
				chastise ("[%s] denominator can't be zero (\"%s\")\n", name, arg);
			op->denominator = tempVal;
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

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if ((!op->windowIsChromosome) && (op->windowSize < 3))
		{
		fprintf (stderr, "[%s] WARNING: raising window size from %d to %d\n",
		                 name, op->windowSize, 3);
		op->windowSize = 3;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_sum));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_window_sum_free--

void op_window_sum_free (dspop* op)
	{
	free (op);
	}


// op_window_sum_apply--

void op_window_sum_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_sum*	op = (dspop_sum*) _op;
	u32			windowSize         = op->windowSize;
	valtype		denominator        = op->denominator;
	valtype		zeroVal            = op->zeroVal;
	int         useActualDenom     = op->useActualDenom;
	valtype		sum;
	u32			startIx, endIx, ix, fillIx;

	if (op->windowIsChromosome) windowSize  = vLen;
	if (op->denomIsWindowSize)  denominator = windowSize;

	sum = 0.0;   // (placate compiler)
	startIx = 0; // (placate compiler)
	for (ix=0 ; ix<vLen ; ix++)
		{
		if (ix % windowSize == 0)
			{
			startIx = ix;
			sum = v[ix];
			}
		else
		    sum += v[ix];

		endIx = ix+1;
		if ((endIx % windowSize == 0) || (endIx == vLen))
			{
			if (useActualDenom)
				v[startIx] = sum / (valtype) (endIx-startIx);
			else
				v[startIx] = sum / denominator;
			for (fillIx=startIx+1 ; fillIx<endIx ; fillIx++)
				v[fillIx] = zeroVal;
			}
		}

	}

//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_sliding_sum--
//	Apply a sliding sum;  the signal is replaced at every entry by the sum over
//	the window centered at that entry.  Input values beyond the ends of the
//	vector are considered to be zero.
//
//----------

// op_sliding_sum_short--

void op_sliding_sum_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "continuous sum over overlapping windows\n");
	}


// op_sliding_sum_usage--

void op_sliding_sum_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sApply a sliding sum;  the signal is replaced at every entry by the sum over\n",   indent);
	fprintf (f, "%sthe window centered at that entry. Input values beyond the ends of the\n",        indent);
	fprintf (f, "%svector are considered to be zero.\n",                                             indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --window=<length>        (W=) size of window\n",                                indent);
	fprintf (f, "%s  --denom=<value>          (D=) denominator\n",                                   indent);
	fprintf (f, "%s                           by default, we do not use a denominator\n",            indent);
	}


// op_sliding_sum_parse--

dspop* op_sliding_sum_parse (char* name, int _argc, char** _argv)
	{
	dspop_sum*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_sum*) malloc (sizeof(dspop_sum));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->windowSize  = (u32) get_named_global ("windowSize", 100);
	op->denominator = 1.0;
	// (op->zeroVal is not used)

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --window=<length> or W=<length>

		if ((strcmp_prefix (arg, "--window=") == 0)
		 || (strcmp_prefix (arg, "W=")        == 0)
		 || (strcmp_prefix (arg, "--W=")      == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt == 0)
				chastise ("[%s] window size can't be zero (\"%s\")\n", name, arg);
			if (tempInt < 0)
				chastise ("[%s] window size can't be negative (\"%s\")\n", name, arg);
			if (tempInt < 3)
				{
				fprintf (stderr, "[%s] WARNING: raising window size from %d to %d\n",
				                 name, tempInt, 3);
				tempInt = 3;
				}
			op->windowSize = (u32) tempInt;
			goto next_arg;
			}

		// --denom[inator]=<value> or Z=<value>

		if ((strcmp_prefix (arg, "--denom=")       == 0)
		 || (strcmp_prefix (arg, "--denominator=") == 0)
		 || (strcmp_prefix (arg, "D=")             == 0)
		 || (strcmp_prefix (arg, "--D=")           == 0))
			{
			if ((strcmp (argVal, "window") == 0)
			 || (strcmp (argVal, "W")      == 0))
				{ op->denominator = op->windowSize;  goto next_arg; }
			tempVal = string_to_valtype (argVal);
			if (tempVal == 0)
				chastise ("[%s] denominator can't be zero (\"%s\")\n", name, arg);
			op->denominator = tempVal;
			goto next_arg;
			}

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if (op->windowSize < 3)
		{
		fprintf (stderr, "[%s] WARNING: raising window size from %d to %d\n",
		                 name, op->windowSize, 3);
		op->windowSize = 3;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_sum));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_sliding_sum_free--

void op_sliding_sum_free (dspop* op)
	{
	free (op);
	}


// op_sliding_sum_apply--
//
// worksheet for sliding sum
//
// vLen       = 20
// windowSize = 9
// hOff       = 4
//
//	v:             *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//	ix:            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
//	cIx:           -  -  -  -  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
//	items in sum:  -  -  -  -  5  6  7  8  9  9  9  9  9  9  9  9  9  9  9  9  8  7  6  5

void op_sliding_sum_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_sum*	op = (dspop_sum*) _op;
	u32			windowSize  = op->windowSize;
	valtype		denominator = op->denominator;
	valtype*	s = get_scratch_vector();
	u32			hOff;
	valtype		sum;
	u32			ix, cIx;

	// compute the sliding sum

	hOff = (windowSize - 1) / 2;

	sum = 0.0;
	for (ix=0 ; ix<vLen+hOff ; ix++)
		{
		//if (dbgSlidingSum)
		//	{
		//	if (ix < vLen)        fprintf (stderr, "r [%d]\n", ix);
		//	if (ix >= windowSize) fprintf (stderr, "r [%d]\n", ix-windowSize);
		//	}
		if (ix < vLen)        sum += v[ix];
		if (ix >= windowSize) sum -= v[ix-windowSize];

		if (ix < hOff) continue;
		cIx = ix - hOff;

		//if (dbgSlidingSum)
		//	fprintf (stderr, "w [%d]\n", cIx);
		s[cIx] = sum;
		}

	// copy scratch array to vector

	for (ix=0 ; ix<vLen ; ix++)
		v[ix] = s[ix] / denominator;

	release_scratch_vector(s);
	}

//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_smooth--
//	Apply a smoothing filter (Hann window).
//
//----------

#define maxWindowSize ((50*1000)+1)

// private dspop subtype

typedef struct dspop_smooth
	{
	dspop		common;			// common elements shared with all operators
	u32			windowSize;
	} dspop_smooth;


// op_smooth_short--

void op_smooth_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "apply a smoothing filter (Hann window)\n");
	}


// op_smooth_usage--

void op_smooth_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sApply a smoothing filter (Hann window);  the signal is replaced at every\n",      indent);
	fprintf (f, "%sentry by the weighted sum over the window centered at that entry. Input\n",       indent);
	fprintf (f, "%svalues beyond the ends of the vector are considered to be zero.\n",               indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --window=<length>        (W=) size of window\n",                                indent);
	fprintf (f, "%s                           (if this is not odd, it will be increased by 1)\n",    indent);
	}


// op_smooth_parse--

dspop* op_smooth_parse (char* name, int _argc, char** _argv)
	{
	dspop_smooth*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_smooth*) malloc (sizeof(dspop_smooth));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->windowSize = (u32) get_named_global ("windowSize", 101);

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --window=<length> or W=<length>

		if ((strcmp_prefix (arg, "--window=") == 0)
		 || (strcmp_prefix (arg, "W=")        == 0)
		 || (strcmp_prefix (arg, "--W=")      == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt == 0)
				chastise ("[%s] window size can't be zero (\"%s\")\n", name, arg);
			if (tempInt < 0)
				chastise ("[%s] window size can't be negative (\"%s\")\n", name, arg);
			if (tempInt > maxWindowSize)
				chastise ("[%s] window size exceeds %u (\"%s\")\n", name, maxWindowSize, arg);
			if (tempInt < 3)
				{
				fprintf (stderr, "[%s] WARNING: raising window size from %d to %d\n",
				                 name, tempInt, 3);
				tempInt = 3;
				}
			if ((tempInt & 1) == 0)
				{
				fprintf (stderr, "[%s] WARNING: raising window size from %d to %d\n",
				                 name, tempInt, tempInt+1);
				tempInt++;
				}
			op->windowSize = (u32) tempInt;
			goto next_arg;
			}

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if ((op->windowSize & 1) == 0)
		{
		fprintf (stderr, "[%s] WARNING: raising window size from %d to %d\n",
		                 name, op->windowSize, op->windowSize+1);
		op->windowSize++;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_smooth));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_smooth_free--

void op_smooth_free (dspop* op)
	{
	free (op);
	}


// op_smooth_apply--

void op_smooth_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_smooth*	op = (dspop_smooth*) _op;
	u32			windowSize  = op->windowSize;
	valtype*	s = get_scratch_vector();
	valtype		window[maxWindowSize];
	double		x;
	valtype		sum;
	u32			hOff, ix, wIx, wStart, wEnd;

	// create window (Hann convolution kernel)

	hOff = (windowSize - 1) / 2;

	for (wIx=0 ; wIx<=hOff ; wIx++)
		{
		x = (wIx+1) / (double) (windowSize+1);
		window[wIx] = window[windowSize-1-wIx] = (1 - cos(2*M_PI*x)) / 2;
		}

	sum = 0.0;
	for (wIx=0 ; wIx<windowSize ; wIx++)
		sum += window[wIx];

	for (wIx=0 ; wIx<windowSize ; wIx++)
		window[wIx] /= sum;

	//for (wIx=0 ; wIx<windowSize ; wIx++)
	//	fprintf (stderr, "window[%d] = %.3f (%.3f)\n",
	//	                 wIx, window[wIx], window[wIx] * sum);

	for (ix=0 ; ix<vLen ; ix++)
		{
		wStart = 0;
		if (ix < hOff) wStart = hOff-ix;

		wEnd = windowSize-1;
		if (ix > vLen + hOff - windowSize) wEnd = vLen-1 + hOff-ix;

		sum = 0.0;
		for (wIx=wStart ; wIx<=wEnd ; wIx++)
			sum += window[wIx] * v[ix-hOff+wIx];
		s[ix] = sum;

		//fprintf (stderr, "s[%d] <-", ix);
		//for (wIx=wStart ; wIx<=wEnd ; wIx++)
		//	fprintf (stderr, " %d*%d", wIx, ix-hOff+wIx);
		//fprintf (stderr, "\n");
		}

	// copy scratch array to vector

	for (ix=0 ; ix<vLen ; ix++)
		v[ix] = s[ix];

	release_scratch_vector(s);
	}

//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_cumulative_sum--
//	Compute the cumulative sum of the current set of interval values.
//
//----------

// op_cumulative_sum_short--

void op_cumulative_sum_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "compute the cumulative sum of the current set of interval values\n");
	}


// op_cumulative_sum_usage--

void op_cumulative_sum_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sCompute the cumulative sum of the current set of interval values,\n",             indent);
	fprintf (f, "%sseparately for each chromsome.\n",                                                indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s\n", indent, name);
	}


// op_cumulative_sum_parse--

dspop* op_cumulative_sum_parse (char* name, int _argc, char** _argv)
	{
	dspop*	op;
	int		argc = _argc;
	char**	argv = _argv;
	char*	arg, *argVal;

	// allocate and initialize our control record

	op = (dspop*) malloc (sizeof(dspop));
	if (op == NULL) goto cant_allocate;

	op->atRandom = false;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	//next_arg:
		argv++;  argc--;
		continue;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_cumulative_sum_free--

void op_cumulative_sum_free (dspop* op)
	{
	free (op);
	}


// op_cumulative_sum_apply--

void op_cumulative_sum_apply
   (arg_dont_complain(dspop*	op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	u32		ix;
	valtype	valSum;

	valSum = 0.0;
	for (ix=0 ; ix<vLen ; ix++)
		{
		valSum += v[ix];
		v[ix]  =  valSum;
		}

	}


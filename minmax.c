// minmax.c-- genodsp operators dealing with minima and maxima

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
#include "minmax.h"

// miscellany

#define min_of(a,b) ((a <= b)? a: b)
#define noIndex ((u32) -1)

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_min_in_interval--
//	Find the minimum value(s) in each of a set of intervals.
//
//	The intervals are replaced by infinity everywhere except at their single
//	trough position.  Ties are broken by the position nearest the center of the
//	interval.  Any bases outside the intervals are also replaced by zeros.
//
//----------

// private dspop subtype

typedef struct dspop_minover
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			originOne;
	valtype		infinityVal;
	int			debug;
	} dspop_minover;

// op_min_in_interval_short--

void op_min_in_interval_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "find the minimum value in each of a set of intervals\n");
	}


// op_min_in_interval_usage--

void op_min_in_interval_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sFind the minimum value in each of a set of intervals. The intervals are\n",       indent);
	fprintf (f, "%sreplaced by infinity everywhere except at their single trough position. Ties\n",  indent);
	fprintf (f, "%sare broken by the position nearest the center of the interval. Any bases\n",      indent);
	fprintf (f, "%soutside the intervals are also replaced by infinity.\n",                          indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --origin=one             input intervals are origin-one, closed\n",             indent);
	fprintf (f, "%s  --origin=zero            input intervals are origin-zero, half-open\n",         indent);
	fprintf (f, "%s  --infinity=<value>       value to fill non-minima\n",                           indent);
	fprintf (f, "%s                           (default is +inf)\n",                                  indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe input intervals are required to be sorted along each chromosome, and\n",      indent);
	fprintf (f, "%snon-overlapping)\n",                                                              indent);
	}


// op_min_in_interval_parse--

dspop* op_min_in_interval_parse (char* name, int _argc, char** _argv)
	{
	dspop_minover*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_minover*) malloc (sizeof(dspop_minover));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename    = NULL;
	op->originOne   = (int) get_named_global ("originOne", false);
	op->infinityVal = valtypeMax;
	op->debug       = false;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --origin=one, --origin=zero

		if ((strcmp (arg, "--origin=one") == 0)
		 || (strcmp (arg, "--origin=1")   == 0))
			{ op->originOne = true;  goto next_arg; }

		if ((strcmp (arg, "--origin=zero") == 0)
		 || (strcmp (arg, "--origin=0")    == 0))
			{ op->originOne = false;  goto next_arg; }

		// --infinity=<value>

		if (strcmp_prefix (arg, "--infinity=") == 0)
			{
			tempVal = string_to_valtype (argVal);
			op->infinityVal = tempVal;
			goto next_arg;
			}

		// --debug arguments

		if (strcmp (arg, "--debug") == 0)
			{ op->debug = true;  goto next_arg; }

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
	                 name, (int) sizeof(dspop_minover));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_min_in_interval_free--

// op_min_in_interval_free--

void op_min_in_interval_free (dspop* _op)
	{
	dspop_minover*	op = (dspop_minover*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_min_in_interval_apply--

void op_min_in_interval_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_minover*	op = (dspop_minover*) _op;
	char*			filename    = op->filename;
	valtype			infinityVal = op->infinityVal;
	FILE*			f;
	char			lineBuffer[1001];
	char			prevChrom[1001];
	valtype*		v = NULL;
	char*			chrom;
	spec*			chromSpec;
	u32				start, end, o, prevEnd, adjStart, adjEnd;
	valtype			minVal, val;
	u32				ix, chromIx, minIx, inset, maxInset;
	int				ok;

	f = fopen (filename, "rt");
	if (f == NULL) goto cant_open_file;

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		chromSpec->flag = false;
		}

	// read intervals and values and find the minima

	if (op->originOne) o = 1;
	              else o = 0;

	prevChrom[0] = 0;
	chromSpec    = NULL;
	prevEnd      = 0;

	v = NULL;
	while (true)
		{
		ok = read_interval (f, lineBuffer, sizeof(lineBuffer), -1,
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
				if ((op->debug) && (prevEnd < chromSpec->length))
					fprintf (stderr, "zeroing %s %u..%u\n", prevChrom, prevEnd, chromSpec->length-1);

				for (ix=prevEnd ; ix<chromSpec->length ; ix++)
					v[ix] = infinityVal;
				}

			v = NULL;
			chromSpec = find_chromosome_spec (chrom);
			if (chromSpec != NULL)
				{
				v = chromSpec->valVector;
				prevEnd = 0;
				if (chromSpec->flag) goto chrom_not_together;
				}

			strncpy (prevChrom, chrom, sizeof(prevChrom)-1);
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

		if ((op->debug) && (prevEnd < adjStart))
			fprintf (stderr, "zeroing %s %u..%u\n", chrom, prevEnd, adjStart-1);

		for (ix=prevEnd ; ix<adjStart ; ix++)
			v[ix] = infinityVal;

		// find the minimum over this interval

		minVal   = v[adjStart];
		minIx    = adjStart;
		maxInset = 0;

		for (ix=adjStart+1 ; ix<adjEnd ; ix++)
			{
			if (v[ix] > minVal) continue;
			if (v[ix] < minVal)
				{
				minVal   = v[ix];
				minIx    = ix;
				maxInset = min_of (minIx-adjStart,adjEnd-minIx);
				continue;
				}
			// we have a tie, which min is closest to center?
			inset = min_of (ix-adjStart,adjEnd-ix);
			if (inset > maxInset)
				{
				minIx    = ix;
				maxInset = inset;
				}
			}

		if (op->debug)
			fprintf (stderr, "minimum over %s %u..%u is " valtypeFmt " at %u\n",
			                 chrom, adjStart, adjEnd-1, minVal, minIx);

		// clear the interval, except for the min

		for (ix=adjStart ; ix<adjEnd ; ix++)
			{ if (ix != minIx) v[ix] = infinityVal; }

		prevEnd = adjEnd;
		}

	// clear any gap at the end of the last chromosome observed

	if (v != NULL)
		{
		if ((op->debug) && (prevEnd < chromSpec->length))
			fprintf (stderr, "zeroing %s %u..%u\n", prevChrom, prevEnd, chromSpec->length-1);

		for (ix=prevEnd ; ix<chromSpec->length ; ix++)
			v[ix] = infinityVal;
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

		if (op->debug)
			fprintf (stderr, "zeroing (all of) %s %u..%u\n",
			                 chrom, chromSpec->start+start, chromSpec->start+end-1);

		for (ix=start ; ix<end ; ix++)
			v[ix] = infinityVal;
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

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_max_in_interval--
//	Find the maximum value(s) in each of a set of intervals.
//
//	The intervals are replaced by zeros everywhere except at their single peak
//	position.  Ties are broken by the position nearest the center of the
//	interval.  Any bases outside the intervals are also replaced by zeros.
//
//----------

// private dspop subtype

typedef struct dspop_maxover
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			originOne;
	valtype		zeroVal;
	int			debug;
	} dspop_maxover;

// op_max_in_interval_short--

void op_max_in_interval_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "find the maximum value in each of a set of intervals\n");
	}


// op_max_in_interval_usage--

void op_max_in_interval_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sFind the maximum value in each of a set of intervals. The intervals are\n",       indent);
	fprintf (f, "%sreplaced by zeros everywhere except at their single peak position. Ties are\n",   indent);
	fprintf (f, "%sbroken by the position nearest the center of the interval. Any bases\n",          indent);
	fprintf (f, "%soutside the intervals are also replaced by zeros.\n",                             indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --origin=one             input intervals are origin-one, closed\n",             indent);
	fprintf (f, "%s  --origin=zero            input intervals are origin-zero, half-open\n",         indent);
	fprintf (f, "%s  --zero=<value>           (Z=) zero value to fill non-maxima\n",                 indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe input intervals are required to be sorted along each chromosome, and\n",      indent);
	fprintf (f, "%snon-overlapping)\n",                                                              indent);
	}


// op_max_in_interval_parse--

dspop* op_max_in_interval_parse (char* name, int _argc, char** _argv)
	{
	dspop_maxover*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_maxover*) malloc (sizeof(dspop_maxover));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename  = NULL;
	op->originOne = (int) get_named_global ("originOne", false);
	op->zeroVal   = 0.0;
	op->debug     = false;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --origin=one, --origin=zero

		if ((strcmp (arg, "--origin=one") == 0)
		 || (strcmp (arg, "--origin=1")   == 0))
			{ op->originOne = true;  goto next_arg; }

		if ((strcmp (arg, "--origin=zero") == 0)
		 || (strcmp (arg, "--origin=0")    == 0))
			{ op->originOne = false;  goto next_arg; }

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
	                 name, (int) sizeof(dspop_maxover));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_max_in_interval_free--

// op_max_in_interval_free--

void op_max_in_interval_free (dspop* _op)
	{
	dspop_maxover*	op = (dspop_maxover*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_max_in_interval_apply--

void op_max_in_interval_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_maxover*	op = (dspop_maxover*) _op;
	char*			filename = op->filename;
	valtype			zeroVal  = op->zeroVal;
	FILE*			f;
	char			lineBuffer[1001];
	char			prevChrom[1001];
	valtype*		v = NULL;
	char*			chrom;
	spec*			chromSpec;
	u32				start, end, o, prevEnd, adjStart, adjEnd;
	valtype			maxVal, val;
	u32				ix, chromIx, maxIx, inset, maxInset;
	int				ok;

	f = fopen (filename, "rt");
	if (f == NULL) goto cant_open_file;

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		chromSpec->flag = false;
		}

	// read intervals and values and find the maxima

	if (op->originOne) o = 1;
	              else o = 0;

	prevChrom[0] = 0;
	chromSpec    = NULL;
	prevEnd      = 0;

	v = NULL;
	while (true)
		{
		ok = read_interval (f, lineBuffer, sizeof(lineBuffer), -1,
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
				if ((op->debug) && (prevEnd < chromSpec->length))
					fprintf (stderr, "zeroing %s %u..%u\n", prevChrom, prevEnd, chromSpec->length-1);

				for (ix=prevEnd ; ix<chromSpec->length ; ix++)
					v[ix] = zeroVal;
				}

			v = NULL;
			chromSpec = find_chromosome_spec (chrom);
			if (chromSpec != NULL)
				{
				v = chromSpec->valVector;
				prevEnd = 0;
				if (chromSpec->flag) goto chrom_not_together;
				}

			strncpy (prevChrom, chrom, sizeof(prevChrom)-1);
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

		if ((op->debug) && (prevEnd < adjStart))
			fprintf (stderr, "zeroing %s %u..%u\n", chrom, prevEnd, adjStart-1);

		for (ix=prevEnd ; ix<adjStart ; ix++)
			v[ix] = zeroVal;

		// find the maximum over this interval

		maxVal   = v[adjStart];
		maxIx    = adjStart;
		maxInset = 0;

		for (ix=adjStart+1 ; ix<adjEnd ; ix++)
			{
			if (v[ix] < maxVal) continue;
			if (v[ix] > maxVal)
				{
				maxVal   = v[ix];
				maxIx    = ix;
				maxInset = min_of (maxIx-adjStart,adjEnd-maxIx);
				continue;
				}
			// we have a tie, which max is closest to center?
			inset = min_of (ix-adjStart,adjEnd-ix);
			if (inset > maxInset)
				{
				maxIx    = ix;
				maxInset = inset;
				}
			}

		if (op->debug)
			fprintf (stderr, "maximum over %s %u..%u is " valtypeFmt " at %u\n",
			                 chrom, adjStart, adjEnd-1, maxVal, maxIx);

		// clear the interval, except for the max

		for (ix=adjStart ; ix<adjEnd ; ix++)
			{ if (ix != maxIx) v[ix] = zeroVal; }

		prevEnd = adjEnd;
		}

	// clear any gap at the end of the last chromosome observed

	if (v != NULL)
		{
		if ((op->debug) && (prevEnd < chromSpec->length))
			fprintf (stderr, "zeroing %s %u..%u\n", prevChrom, prevEnd, chromSpec->length-1);

		for (ix=prevEnd ; ix<chromSpec->length ; ix++)
			v[ix] = zeroVal;
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

		if (op->debug)
			fprintf (stderr, "zeroing (all of) %s %u..%u\n",
			                 chrom, chromSpec->start+start, chromSpec->start+end-1);

		for (ix=start ; ix<end ; ix++)
			v[ix] = zeroVal;
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

//----------
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_local_minima--
//	Apply a trough-finder (more of a local minima finder).  Any entry equal to
//	to the minimum (in the window centered at that entry) is kept.  All entries
//	that aren't the minimum are cleared.
//
//----------

// private dspop subtype

typedef struct dspop_localmin
	{
	dspop		common;			// common elements shared with all operators
	u32			neighborhood;
	valtype		infinityVal;
	} dspop_localmin;


// op_local_minima_short--

void op_local_minima_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "find local minima\n");
	}


// op_local_minima_usage--

void op_local_minima_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sIdentify local minima. Any entry *not* equal to the minimum over the\n",          indent);
	fprintf (f, "%sneighborhood centered at that entry is set to infinity. Input values beyond\n",   indent);
	fprintf (f, "%sthe ends of the vector are not considered.\n",                                    indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --neighborhood=<length>  (N=) size of neighborhood\n",                          indent);
	fprintf (f, "%s                           (default is 3)\n",                                     indent);
	fprintf (f, "%s                           (if this is not odd, it will be increased by 1)\n",    indent);
	fprintf (f, "%s  --infinity=<value>       value to fill non-minima\n",                           indent);
	fprintf (f, "%s                           (default is +inf)\n",                                  indent);
	}


// op_local_minima_parse--

dspop* op_local_minima_parse (char* name, int _argc, char** _argv)
	{
	dspop_localmin*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_localmin*) malloc (sizeof(dspop_localmin));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->neighborhood = 3;
	op->infinityVal  = valtypeMax;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --neighborhood=<length> or N=<length>

		if ((strcmp_prefix (arg, "--neighborhood=") == 0)
		 || (strcmp_prefix (arg, "N=")              == 0)
		 || (strcmp_prefix (arg, "--N=")            == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt == 0)
				chastise ("[%s] neighborhood can't be zero (\"%s\")\n", name, arg);
			if (tempInt < 0)
				chastise ("[%s] neighborhood can't be negative (\"%s\")\n", name, arg);
			if (tempInt < 3)
				{
				fprintf (stderr, "[%s] WARNING: raising neighborhood from %d to %d\n",
				                 name, tempInt, 3);
				tempInt = 3;
				}
			if ((tempInt & 1) == 0)
				{
				fprintf (stderr, "[%s] WARNING: raising neighborhood from %d to %d\n",
				                 name, tempInt, tempInt+1);
				tempInt++;
				}
			op->neighborhood = (u32) tempInt;
			goto next_arg;
			}

		// --infinity=<value>

		if (strcmp_prefix (arg, "--infinity=") == 0)
			{
			tempVal = string_to_valtype (argVal);
			op->infinityVal = tempVal;
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

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_localmin));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_local_minima_free--

void op_local_minima_free (dspop* op)
	{
	free (op);
	}


// op_local_minima_apply--

void op_local_minima_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_localmin*	op = (dspop_localmin*) _op;
	u32			neighborhood = op->neighborhood;
	valtype		infinityVal  = op->infinityVal;
	valtype*	s = get_scratch_vector();
	valtype		val;
	u32			hOff, ix, wIx, wStart, wEnd;

	// find the local minima

	hOff = (neighborhood - 1) / 2;

	for (ix=0 ; ix<vLen ; ix++)
		{
		if (ix < hOff) wStart = 0;
		          else wStart = ix - hOff;  

		if (ix + hOff >= vLen) wEnd = vLen - 1;
		                  else wEnd = ix + hOff;

		val = v[ix];
		for (wIx=wStart ; wIx<=wEnd ; wIx++)
			{
			if ((wIx != ix) && (v[wIx] < val))
				{ val = infinityVal;  break; }
			}

		s[ix] = val;
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
// op_local_maxima--
//	Apply a peak-finder (more of a local maxima finder).  Any entry equal to
//	to the maximum (in the window centered at that entry) is kept.  All entries
//	that aren't the maximum are cleared.
//
//----------

// private dspop subtype

typedef struct dspop_localmax
	{
	dspop		common;			// common elements shared with all operators
	u32			neighborhood;
	valtype		zeroVal;
	} dspop_localmax;


// op_local_maxima_short--

void op_local_maxima_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "find local maxima\n");
	}


// op_local_maxima_usage--

void op_local_maxima_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sIdentify local maxima. Any entry *not* equal to the maximum over the\n",          indent);
	fprintf (f, "%sneighborhood centered at that entry is set to zero. Input values beyond the\n",   indent);
	fprintf (f, "%sends of the vector are not considered.\n",                                        indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --neighborhood=<length>  (N=) size of neighborhood\n",                          indent);
	fprintf (f, "%s                           (default is 3)\n",                                     indent);
	fprintf (f, "%s                           (if this is not odd, it will be increased by 1)\n",    indent);
	fprintf (f, "%s  --zero=<value>           (Z=) zero value to fill non-maxima\n",                 indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	}


// op_local_maxima_parse--

dspop* op_local_maxima_parse (char* name, int _argc, char** _argv)
	{
	dspop_localmax*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_localmax*) malloc (sizeof(dspop_localmax));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->neighborhood = 3;
	op->zeroVal      = 0.0;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --neighborhood=<length> or N=<length>

		if ((strcmp_prefix (arg, "--neighborhood=") == 0)
		 || (strcmp_prefix (arg, "N=")              == 0)
		 || (strcmp_prefix (arg, "--N=")            == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt == 0)
				chastise ("[%s] neighborhood can't be zero (\"%s\")\n", name, arg);
			if (tempInt < 0)
				chastise ("[%s] neighborhood can't be negative (\"%s\")\n", name, arg);
			if (tempInt < 3)
				{
				fprintf (stderr, "[%s] WARNING: raising neighborhood from %d to %d\n",
				                 name, tempInt, 3);
				tempInt = 3;
				}
			if ((tempInt & 1) == 0)
				{
				fprintf (stderr, "[%s] WARNING: raising neighborhood from %d to %d\n",
				                 name, tempInt, tempInt+1);
				tempInt++;
				}
			op->neighborhood = (u32) tempInt;
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

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_localmax));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_local_maxima_free--

void op_local_maxima_free (dspop* op)
	{
	free (op);
	}


// op_local_maxima_apply--

void op_local_maxima_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_localmax*	op = (dspop_localmax*) _op;
	u32			neighborhood = op->neighborhood;
	valtype		zeroVal      = op->zeroVal;
	valtype*	s = get_scratch_vector();
	valtype		val;
	u32			hOff, ix, wIx, wStart, wEnd;

	// find the local maxima
	// $$$ this assumes an odd-sized window, but I don't think anything
	// $$$ .. enforces that;  I need to check this an other operators for this
	// $$$ .. issue

	hOff = (neighborhood - 1) / 2;

	for (ix=0 ; ix<vLen ; ix++)
		{
		if (ix < hOff) wStart = 0;
		          else wStart = ix - hOff;  

		if (ix + hOff >= vLen) wEnd = vLen - 1;
		                  else wEnd = ix + hOff;

		val = v[ix];
		for (wIx=wStart ; wIx<=wEnd ; wIx++)
			{
			if ((wIx != ix) && (v[wIx] > val))
				{ val = zeroVal;  break; }
			}

		s[ix] = val;
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
// op_best_local_min--
//	Apply a best-local-min filter.  Each entry is replaced by the minimum in
//	the window centered at that entry.
//
//----------

// private dspop subtype

typedef struct dspop_bestmin
	{
	dspop		common;			// common elements shared with all operators
	u32			windowSize;
	int			debug;
	} dspop_bestmin;


// op_best_local_min_short--

void op_best_local_min_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "replace each entry with the minimum value in its local window\n");
	}


// op_best_local_min_usage--

void op_best_local_min_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sReplace each entry with the minimum value in the window centered at that\n",      indent);
	fprintf (f, "%sentry. Input values beyond the ends of the vector are not considered.\n",         indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --window=<length>  (W=) size of window\n",                                      indent);
	}


// op_best_local_min_parse--

dspop* op_best_local_min_parse (char* name, int _argc, char** _argv)
	{
	dspop_bestmin*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_bestmin*) malloc (sizeof(dspop_bestmin));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->windowSize = (u32) get_named_global ("windowSize", 100);
	op->debug      = false;

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

		// --debug arguments

		if (strcmp (arg, "--debug") == 0)
			{ op->debug = true;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_bestmin));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_best_local_min_free--

void op_best_local_min_free (dspop* op)
	{
	free (op);
	}


// op_best_local_min_apply--

void op_best_local_min_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_bestmin*	op = (dspop_bestmin*) _op;
	u32			windowSize = op->windowSize;
	valtype*	s = get_scratch_vector();
	valtype		minVal;
	u32			wLft, wRgt, ix, ixLft, ixRgt, wIx, bestIx;

	//////////
	// find the minimum in each window
	//////////

	wLft = (windowSize - 1) / 2;
	wRgt = (windowSize - 1) - wLft;

	// find min for window centered at 0

	minVal = v[0];
	for (wIx=1 ; wIx<=wRgt ; wIx++)
		{
		if (wIx >= vLen) break;
		if (v[wIx] < minVal) minVal = v[wIx];
		}

	s[0] = minVal;

	if (op->debug)
		{
		bestIx = noIndex;
		for (wIx=0 ; wIx<=wRgt ; wIx++)
			{
			if (wIx >= vLen) break;
			if (v[wIx] == minVal) bestIx = wIx;
			}
		fprintf (stderr, "min[%u]=" valtypeFmt ", from initial window [%u]\n",
		                 0, minVal, bestIx);
		}

	// for each position ix (except 0), find min for window centered at ix;  at
	// each ix we try to make use of our knowledge of the min over the window
	// centered at ix-1

	for (ix=1 ; ix<vLen ; ix++)
		{
		if (ix       < wLft+1) ixLft = noIndex;
		                  else ixLft = ix - (wLft+1);
		if (ix + wRgt >= vLen) ixRgt = noIndex;
		                  else ixRgt = ix + wRgt;

		// if the new value being added to the window is as small as any in the
		// old window, it's the min of the new window

		if ((ixRgt != noIndex) && (v[ixRgt] <= minVal))
			{
			minVal = v[ixRgt];
			if (op->debug)
				fprintf (stderr, "min[%u]=" valtypeFmt ", from new value [%u]\n",
				                 ix, minVal, ixRgt);
			}

		// otherwise (the new value is bigger than min of the old window), if
		// the old value being removed from the window was not the min of that
		// window, the min of the new window is the same as the min of the old 

		else if ((ixLft == noIndex) || (v[ixLft] > minVal))
			{
			if (op->debug)
				fprintf (stderr, "min[%u]=" valtypeFmt ", from old window\n",
				                 ix, minVal);
			}

		// otherwise, we need to search the new window to find the min

		else
			{
			if (ixLft == noIndex) ixLft = 0;
			                 else ixLft = ixLft+1;
			if (ixRgt == noIndex) ixRgt = vLen-1;
			minVal = v[ixLft];
			for (wIx=ixLft+1 ; wIx<=ixRgt ; wIx++)
				{ if (v[wIx] < minVal) minVal = v[wIx]; }

			if (op->debug)
				{
				bestIx = noIndex;
				for (wIx=ixLft ; wIx<=ixRgt ; wIx++)
					{ if (v[wIx] == minVal) bestIx = wIx; }
				fprintf (stderr, "min[%u]=" valtypeFmt ", from new window [%u]\n",
				                 ix, minVal, bestIx);
				}
			}

		s[ix] = minVal;
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
// op_best_local_max--
//	Apply a best-local-max filter.  Each entry is replaced by the maximum in
//	the window centered at that entry.
//
//----------

// private dspop subtype

typedef struct dspop_bestmax
	{
	dspop		common;			// common elements shared with all operators
	u32			windowSize;
	int			debug;
	} dspop_bestmax;


// op_best_local_max_short--

void op_best_local_max_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "replace each entry with the maximum value in its local window\n");
	}


// op_best_local_max_usage--

void op_best_local_max_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sReplace each entry with the maximum value in the window centered at that\n",      indent);
	fprintf (f, "%sentry. Input values beyond the ends of the vector are not considered.\n",         indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --window=<length>  (W=) size of window\n",                                      indent);
	}


// op_best_local_max_parse--

dspop* op_best_local_max_parse (char* name, int _argc, char** _argv)
	{
	dspop_bestmax*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_bestmax*) malloc (sizeof(dspop_bestmax));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->windowSize = (u32) get_named_global ("windowSize", 100);
	op->debug      = false;

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

		// --debug arguments

		if (strcmp (arg, "--debug") == 0)
			{ op->debug = true;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_bestmax));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_best_local_max_free--

void op_best_local_max_free (dspop* op)
	{
	free (op);
	}


// op_best_local_max_apply--

void op_best_local_max_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_bestmax*	op = (dspop_bestmax*) _op;
	u32			windowSize = op->windowSize;
	valtype*	s = get_scratch_vector();
	valtype		maxVal;
	u32			wLft, wRgt, ix, ixLft, ixRgt, wIx, bestIx;

	//////////
	// find the maximum in each window
	//////////

	wLft = (windowSize - 1) / 2;
	wRgt = (windowSize - 1) - wLft;

	// find max for window centered at 0

	maxVal = v[0];
	for (wIx=1 ; wIx<=wRgt ; wIx++)
		{
		if (wIx >= vLen) break;
		if (v[wIx] > maxVal) maxVal = v[wIx];
		}

	s[0] = maxVal;

	if (op->debug)
		{
		bestIx = noIndex;
		for (wIx=0 ; wIx<=wRgt ; wIx++)
			{
			if (wIx >= vLen) break;
			if (v[wIx] == maxVal) bestIx = wIx;
			}
		fprintf (stderr, "max[%u]=" valtypeFmt ", from initial window [%u]\n",
		                 0, maxVal, bestIx);
		}

	// for each position ix (except 0), find max for window centered at ix;  at
	// each ix we try to make use of our knowledge of the max over the window
	// centered at ix-1

	for (ix=1 ; ix<vLen ; ix++)
		{
		if (ix       < wLft+1) ixLft = noIndex;
		                  else ixLft = ix - (wLft+1);
		if (ix + wRgt >= vLen) ixRgt = noIndex;
		                  else ixRgt = ix + wRgt;

		// if the new value being added to the window is as big as any in the
		// old window, it's the max of the new window

		if ((ixRgt != noIndex) && (v[ixRgt] >= maxVal))
			{
			maxVal = v[ixRgt];
			if (op->debug)
				fprintf (stderr, "max[%u]=" valtypeFmt ", from new value [%u]\n",
				                 ix, maxVal, ixRgt);
			}

		// otherwise (the new value is smaller than max of the old window), if
		// the old value being removed from the window was not the max of that
		// window, the max of the new window is the same as the max of the old 

		else if ((ixLft == noIndex) || (v[ixLft] < maxVal))
			{
			if (op->debug)
				fprintf (stderr, "max[%u]=" valtypeFmt ", from old window\n",
				                 ix, maxVal);
			}

		// otherwise, we need to search the new window to find the max

		else
			{
			if (ixLft == noIndex) ixLft = 0;
			                 else ixLft = ixLft+1;
			if (ixRgt == noIndex) ixRgt = vLen-1;
			maxVal = v[ixLft];
			for (wIx=ixLft+1 ; wIx<=ixRgt ; wIx++)
				{ if (v[wIx] > maxVal) maxVal = v[wIx]; }

			if (op->debug)
				{
				bestIx = noIndex;
				for (wIx=ixLft ; wIx<=ixRgt ; wIx++)
					{ if (v[wIx] == maxVal) bestIx = wIx; }
				fprintf (stderr, "max[%u]=" valtypeFmt ", from new window [%u]\n",
				                 ix, maxVal, bestIx);
				}
			}

		s[ix] = maxVal;
		}

	// copy scratch array to vector

	for (ix=0 ; ix<vLen ; ix++)
		v[ix] = s[ix];

	release_scratch_vector(s);
	}

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_min_with--
//	Combine an incoming set of interval values (read from a file) with the
//	current set, keeping the minimum value at each location.
//
//----------

// private dspop subtype

typedef struct dspop_min_with
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			originOne;
	int			destroyFile;
	} dspop_min_with;


// op_min_with_short--

void op_min_with_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "take the minimum of an incoming set of interval values and the current set\n");
	}


// op_min_with_usage--

void op_min_with_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sTake the minimum of an incoming set of interval values and the current set.\n",   indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --value=<col>            input intervals contain a value in the specified\n",   indent);
	fprintf (f, "%s                           column\n",                                             indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	fprintf (f, "%s  --destroy                destroy the file after reading it\n",                  indent);
	fprintf (f, "%s                           (BE SURE THAT'S WHAT YOU WANT, IT CAN'T BE UNDONE)\n", indent);
	}


// op_min_with_parse--

dspop* op_min_with_parse (char* name, int _argc, char** _argv)
	{
	dspop_min_with*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_min_with*) malloc (sizeof(dspop_min_with));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename    = NULL;
	op->valColumn   = (int) get_named_global ("valColumn", 4-1);
	op->originOne   = (int) get_named_global ("originOne", false);
	op->destroyFile = false;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --value=<col>

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

		// --destroy

		if (strcmp (arg, "--destroy") == 0)
			{ op->destroyFile = true;  goto next_arg; }

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
	                 name, (int) sizeof(dspop_min_with));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_min_with_free--

void op_min_with_free (dspop* _op)
	{
	dspop_min_with*	op = (dspop_min_with*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_min_with_apply--

void op_min_with_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_min_with*	op = (dspop_min_with*) _op;
	char*		filename = op->filename;
	FILE*		f;
	char		lineBuffer[1001];
	char		prevChrom[1001];
	valtype*	v = NULL;
	char*		chrom;
	spec*		chromSpec;
	u32			start, end, o, adjStart, adjEnd;
	valtype		val;
	u32			ix, chromIx;
	int			ok;

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

	// read intervals and values and add them

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

		// add over this interval

		for (ix=adjStart ; ix<adjEnd ; ix++)
			{ if (val < v[ix]) v[ix] = val; }
		}

	// success

	fclose (f);

	if (op->destroyFile)
		remove (filename);

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
// op_max_with--
//	Combine an incoming set of interval values (read from a file) with the
//	current set, keeping the maximum value at each location.
//
//----------

// private dspop subtype

typedef struct dspop_max_with
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			originOne;
	int			destroyFile;
	} dspop_max_with;


// op_max_with_short--

void op_max_with_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "take the maximum of an incoming set of interval values and the current set\n");
	}


// op_max_with_usage--

void op_max_with_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sTake the maximum of an incoming set of interval values and the current set.\n",   indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --value=<col>            input intervals contain a value in the specified\n",   indent);
	fprintf (f, "%s                           column\n",                                             indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	fprintf (f, "%s  --destroy                destroy the file after reading it\n",                  indent);
	fprintf (f, "%s                           (BE SURE THAT'S WHAT YOU WANT, IT CAN'T BE UNDONE)\n", indent);
	}


// op_max_with_parse--

dspop* op_max_with_parse (char* name, int _argc, char** _argv)
	{
	dspop_max_with*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_max_with*) malloc (sizeof(dspop_max_with));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename    = NULL;
	op->valColumn   = (int) get_named_global ("valColumn", 4-1);
	op->originOne   = (int) get_named_global ("originOne", false);
	op->destroyFile = false;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --value=<col>

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

		// --destroy

		if (strcmp (arg, "--destroy") == 0)
			{ op->destroyFile = true;  goto next_arg; }

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
	                 name, (int) sizeof(dspop_max_with));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_max_with_free--

void op_max_with_free (dspop* _op)
	{
	dspop_max_with*	op = (dspop_max_with*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_max_with_apply--

void op_max_with_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_max_with*	op = (dspop_max_with*) _op;
	char*		filename = op->filename;
	FILE*		f;
	char		lineBuffer[1001];
	char		prevChrom[1001];
	valtype*	v = NULL;
	char*		chrom;
	spec*		chromSpec;
	u32			start, end, o, adjStart, adjEnd;
	valtype		val;
	u32			ix, chromIx;
	int			ok;

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

	// read intervals and values and add them

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

		// add over this interval

		for (ix=adjStart ; ix<adjEnd ; ix++)
			{ if (val > v[ix]) v[ix] = val; }
		}

	// success

	fclose (f);

	if (op->destroyFile)
		remove (filename);

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


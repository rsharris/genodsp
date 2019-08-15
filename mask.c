// mask.c-- genodsp operators performing interval masking

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
#include "mask.h"

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_mask--
//	Apply a list of masking intervals (read from a file).
//
//----------

// private dspop subtype

typedef struct dspop_mask
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			haveMaskVal;
	char*		maskValVarName;
	valtype		maskVal;
	int			originOne;
	} dspop_mask;


// op_mask_short--

void op_mask_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "apply a list of masking intervals (read from a file)\n");
	}


// op_mask_usage--

void op_mask_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sApply a list of masking intervals (read from a file).\n",                         indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --mask=<value>           (M=) value to fill masked intervals with\n",           indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                    indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	}


// op_mask_parse--

dspop* op_mask_parse (char* name, int _argc, char** _argv)
	{
	dspop_mask*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;

	// allocate and initialize our control record

	op = (dspop_mask*) malloc (sizeof(dspop_mask));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename       = NULL;
	op->haveMaskVal    = false;
	op->maskValVarName = NULL;
	op->maskVal        = 0.0;
	op->originOne      = (int) get_named_global ("originOne", false);

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --mask=<value> or M=<value>

		if ((strcmp_prefix (arg, "--mask=") == 0)
		 || (strcmp_prefix (arg, "M=")      == 0)
		 || (strcmp_prefix (arg, "--M=")    == 0))
			{
			if (op->haveMaskVal) goto more_than_one_mask_val;
			if (!try_string_to_valtype (argVal, &op->maskVal))
				op->maskValVarName = copy_string (argVal);
			else
				op->maskVal = string_to_valtype (argVal);
			op->haveMaskVal = true;
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
	                 name, (int) sizeof(dspop_mask));
	exit(EXIT_FAILURE);

more_than_one_mask_val:
	fprintf (stderr, "[%s] mask value specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)


filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_mask_free--

void op_mask_free (dspop* _op)
	{
	dspop_mask*	op = (dspop_mask*) _op;

	if (op->filename       != NULL) free (op->filename);
	if (op->maskValVarName != NULL) free (op->maskValVarName);
	free (op);
	}


// op_mask_apply--

void op_mask_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_mask*	op = (dspop_mask*) _op;
	char*		filename = op->filename;
	valtype		maskVal  = op->maskVal;
	FILE*		f;
	char		lineBuffer[1000];
	char		prevChrom[1000];
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

	// if the mask value is a named variable, fetch it now;  note that we copy
	// the value from the named variable, then destroy our reference to the
	// named variable

	if (op->maskValVarName != NULL)
		{
		ok = named_global_exists (op->maskValVarName, &maskVal);
		if (!ok) goto no_mask_val;
		op->maskVal = maskVal;
		fprintf (stderr, "[%s] using %s = " valtypeFmt " as mask value\n",
		                 _op->name, op->maskValVarName, maskVal);
		free (op->maskValVarName);
		op->maskValVarName = NULL;
		}

	// read intervals and apply them as masks

	if (op->originOne) o = 1;
	              else o = 0;

	prevChrom[0] = 0;
	chromSpec    = NULL;

	v = NULL;
	while (true)
		{
		ok = read_interval (f, lineBuffer, sizeof(lineBuffer), -1,
	                        &chrom, &start, &end, &val);
		if (!ok) break;

		//fprintf (stderr, "%s %u %u %f\n", chrom, start, end, val);

		if (strcmp (chrom, prevChrom) != 0)
			{
			v = NULL;
			chromSpec = find_chromosome_spec (chrom);
			if (chromSpec != NULL) v = chromSpec->valVector;
			strncpy (prevChrom, chrom, sizeof(prevChrom));
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

		for (ix=adjStart ; ix<adjEnd ; ix++)
			v[ix] = maskVal;
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

no_mask_val:
	fprintf (stderr, "[%s] attempt to use %s as mask value failed (no such variable)\n",
	                 _op->name, op->maskValVarName);
	exit(EXIT_FAILURE);

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
// op_mask_not--
//	Apply the complement of a list of masking intervals (read from a file).
//
// This amounts to clearing any intervals absent from the incoming set.
//
//----------

// private dspop subtype

typedef struct dspop_masknot
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	valtype		maskVal;
	int			originOne;
	int			debug;
	} dspop_masknot;


// op_mask_not_short--

void op_mask_not_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "apply the complement of a list of masking intervals (read from a file)\n");
	}


// op_mask_not_usage--

void op_mask_not_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sApply the complement of a list of masking intervals (read from a file). This\n",  indent);
	fprintf (f, "%samounts to clearing any intervals absent from the incoming set.\n",               indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --mask=<value>           (M=) value to fill masked intervals with\n",           indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe input intervals are required to be sorted along each chromosome, and\n",      indent);
	fprintf (f, "%snon-overlapping)\n",                                                              indent);
	}


// op_mask_not_parse--

dspop* op_mask_not_parse (char* name, int _argc, char** _argv)
	{
	dspop_masknot*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_masknot*) malloc (sizeof(dspop_masknot));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename  = NULL;
	op->maskVal   = 0.0;
	op->originOne = (int) get_named_global ("originOne", false);

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --mask=<value> or M=<value>

		if ((strcmp_prefix (arg, "--mask=") == 0)
		 || (strcmp_prefix (arg, "M=")      == 0)
		 || (strcmp_prefix (arg, "--M=")    == 0))
			{
			tempVal = string_to_valtype (argVal);
			op->maskVal = tempVal;
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
	                 name, (int) sizeof(dspop_masknot));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_mask_not_free--

void op_mask_not_free (dspop* _op)
	{
	dspop_masknot*	op = (dspop_masknot*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_mask_not_apply--

void op_mask_not_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_masknot*	op = (dspop_masknot*) _op;
	char*		filename = op->filename;
	valtype		maskVal  = op->maskVal;
	FILE*		f;
	char		lineBuffer[1000];
	char		prevChrom[1000];
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

	// read intervals and values and multiply by them

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
				for (ix=prevEnd ; ix<chromSpec->length ; ix++)
					v[ix] = maskVal;
				}

			v = NULL;
			chromSpec = find_chromosome_spec (chrom);
			if (chromSpec != NULL)
				{
				v = chromSpec->valVector;
				prevEnd = 0;
				if (chromSpec->flag) goto chrom_not_together;
				}

			strncpy (prevChrom, chrom, sizeof(prevChrom));
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
			v[ix] = maskVal;

		// "fail to mask" over this interval;  actually, there is nothing to
		// do-- the surrounding code will clear all intervals *not* present in
		// the input

		; // (do nothing)

		prevEnd = adjEnd;
		}

	// clear any gap at the end of the last chromosome observed

	if (v != NULL)
		{
		for (ix=prevEnd ; ix<chromSpec->length ; ix++)
			v[ix] = maskVal;
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
			v[ix] = maskVal;
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
// op_clip--
//	Clip the current set of interval values to a specified minimum and/or
//	maximum.
//
//----------

// private dspop subtype

typedef struct dspop_clip
	{
	dspop		common;			// common elements shared with all operators
	int			haveMinVal;
	char*		minValVarName;
	valtype		minVal;
	int			haveMaxVal;
	char*		maxValVarName;
	valtype		maxVal;
	} dspop_clip;

// op_clip_short--

void op_clip_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "clip the current set of interval values to a specified min and/or max\n");
	}


// op_clip_usage--

void op_clip_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sClip the current set of interval values to a specified minimum and/or.\n",        indent);
	fprintf (f, "%smaximum.\n",                                                                      indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --min=<value>            minimum value in resulting signal\n",                  indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                    indent);
	fprintf (f, "%s  --max=<value>            maximum value in resulting signal\n",                  indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                    indent);
	}


// op_clip_parse--

dspop* op_clip_parse (char* name, int _argc, char** _argv)
	{
	dspop_clip*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;

	// allocate and initialize our control record

	op = (dspop_clip*) malloc (sizeof(dspop_clip));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->haveMinVal    = false;
	op->minValVarName = NULL;
	op->minVal        = 0.0;
	op->haveMaxVal    = false;
	op->maxValVarName = NULL;
	op->maxVal        = 0.0;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --min=<value>

		if ((strcmp_prefix (arg, "--minimum=") == 0)
		 || (strcmp_prefix (arg, "--min=")     == 0))
			{
			if (op->haveMinVal) goto more_than_one_minimum;
			if (!try_string_to_valtype (argVal, &op->minVal))
				op->minValVarName = copy_string (argVal);
			op->haveMinVal = true;
			goto next_arg;
			}

		// --max=<value>

		if ((strcmp_prefix (arg, "--maximum=") == 0)
		 || (strcmp_prefix (arg, "--max=")     == 0))
			{
			if (op->haveMaxVal) goto more_than_one_maximum;
			if (!try_string_to_valtype (argVal, &op->maxVal))
				op->maxValVarName = copy_string (argVal);
			op->haveMaxVal = true;
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

	if ((!op->haveMaxVal) && (!op->haveMinVal)) goto limits_missing;

	if ((op->haveMaxVal) && (op->haveMinVal))
		{ if (op->minVal > op->maxVal) goto conflicting_limits; }

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_clip));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_minimum:
	fprintf (stderr, "[%s] minimum limit specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_maximum:
	fprintf (stderr, "[%s] maximum limit specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

limits_missing:
	fprintf (stderr, "[%s] neither minimum nor maximum limit was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

conflicting_limits:
	fprintf (stderr, "[%s] conflicting limits (" valtypeFmt ">" valtypeFmt ")\n",
	                 name, op->minVal, op->maxVal);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_clip_free--

void op_clip_free (dspop* _op)
	{
	dspop_clip*	op = (dspop_clip*) _op;

	if (op->minValVarName != NULL) free (op->minValVarName);
	if (op->maxValVarName != NULL) free (op->maxValVarName);
	free (op);
	}


// op_clip_apply--

void op_clip_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_clip*	op = (dspop_clip*) _op;
	valtype		minVal = op->minVal;
	valtype		maxVal = op->maxVal;
	int			ok;
	u32			ix;

	// if either limit is a named variable, fetch it now;  note that we copy
	// the value from the named variable, then destroy our reference to the
	// named variable

	if (op->minValVarName != NULL)
		{
		ok = named_global_exists (op->minValVarName, &minVal);
		if (!ok) goto no_minimum;
		op->minVal = minVal;
		fprintf (stderr, "[%s] using %s = " valtypeFmt " as minimum limit\n",
		                 _op->name, op->minValVarName, minVal);
		free (op->minValVarName);
		op->minValVarName = NULL;
		}

	if (op->maxValVarName != NULL)
		{
		ok = named_global_exists (op->maxValVarName, &maxVal);
		if (!ok) goto no_maximum;
		op->maxVal = maxVal;
		fprintf (stderr, "[%s] using %s = " valtypeFmt " as maximum limit\n",
		                 _op->name, op->maxValVarName, maxVal);
		free (op->maxValVarName);
		op->maxValVarName = NULL;
		}

	// apply limits over the vector

	if (!op->haveMaxVal)
		{ // clip to minimum only
		for (ix=0 ; ix<vLen ; ix++)
			{ if (v[ix] < minVal) v[ix] = minVal; }
		}
	else if (!op->haveMinVal)
		{ // clip to maximum only
		for (ix=0 ; ix<vLen ; ix++)
			{ if (v[ix] > maxVal) v[ix] = maxVal; }
		}
	else
		{ // clip to minimum and maximum
		for (ix=0 ; ix<vLen ; ix++)
			{
			if      (v[ix] < minVal) v[ix] = minVal;
			else if (v[ix] > maxVal) v[ix] = maxVal;
			}
		}

	// success

	return;

	// failure

no_minimum:
	fprintf (stderr, "[%s] attempt to use %s as minimum failed (no such variable)\n",
	                 _op->name, op->minValVarName);
	exit(EXIT_FAILURE);

no_maximum:
	fprintf (stderr, "[%s] attempt to use %s as maximum failed (no such variable)\n",
	                 _op->name, op->minValVarName);
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
// op_erase--
//	Erase any values in the current set of interval values that occur within
//	(or, alternatively, outside of) a specified minimum and/or maximum.
//
//----------

// private dspop subtype

typedef struct dspop_erase
	{
	dspop		common;			// common elements shared with all operators
	int			haveMinVal;
	char*		minValVarName;
	valtype		minVal;
	int			haveMaxVal;
	char*		maxValVarName;
	valtype		maxVal;
	int			keepInside;
	valtype		zeroVal;
	} dspop_erase;

// op_erase_short--

void op_erase_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "erase values within (or outside of) a specified min and/or max\n");
	}


// op_erase_usage--

void op_erase_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%serase any values in the current set of interval values that occur within a\n",    indent);
	fprintf (f, "%sspecified minimum and/or maximum.\n",                                             indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  --min=<value>            minimum value in resulting signal\n",                  indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                    indent);
	fprintf (f, "%s                           (default is negative infinity)\n",                     indent);
	fprintf (f, "%s  --max=<value>            maximum value in resulting signal\n",                  indent);
	fprintf (f, "%s                           <value> can be a named variable\n",                    indent);
	fprintf (f, "%s                           (default is infinity)\n",                              indent);
	fprintf (f, "%s  --keep:outside           keep values outside the min and max; erase values\n",  indent);
	fprintf (f, "%s                           inside\n",                                             indent);
	fprintf (f, "%s                           (this is the default)\n",                              indent);
	fprintf (f, "%s  --keep:inside            keep values within the min and max; erase values\n",   indent);
	fprintf (f, "%s                           outside\n",                                            indent);
	fprintf (f, "%s  --zero=<value>           (Z=) value for erased locations\n",                    indent);
	fprintf (f, "%s                           (default is 0.0)\n",                                   indent);
	}


// op_erase_parse--

dspop* op_erase_parse (char* name, int _argc, char** _argv)
	{
	dspop_erase*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_erase*) malloc (sizeof(dspop_erase));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->haveMinVal    = false;
	op->minValVarName = NULL;
	op->minVal        = 0.0;
	op->haveMaxVal    = false;
	op->maxValVarName = NULL;
	op->maxVal        = 0.0;
	op->keepInside    = false;
	op->zeroVal       = 0.0;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --min=<value>

		if ((strcmp_prefix (arg, "--minimum=") == 0)
		 || (strcmp_prefix (arg, "--min=")     == 0))
			{
			if (op->haveMinVal) goto more_than_one_minimum;
			if (!try_string_to_valtype (argVal, &op->minVal))
				op->minValVarName = copy_string (argVal);
			op->haveMinVal = true;
			goto next_arg;
			}

		// --max=<value>

		if ((strcmp_prefix (arg, "--maximum=") == 0)
		 || (strcmp_prefix (arg, "--max=")     == 0))
			{
			if (op->haveMaxVal) goto more_than_one_maximum;
			if (!try_string_to_valtype (argVal, &op->maxVal))
				op->maxValVarName = copy_string (argVal);
			op->haveMaxVal = true;
			goto next_arg;
			}

		// --keep:outside and --keep:inside

		if ((strcmp (arg, "--keep:outside") == 0)
		 || (strcmp (arg, "--keep=outside") == 0))
			{
			op->keepInside = false;
			goto next_arg;
			}

		if ((strcmp (arg, "--keep:inside") == 0)
		 || (strcmp (arg, "--keep=inside") == 0))
			{
			op->keepInside = true;
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

	if ((!op->haveMaxVal) && (!op->haveMinVal)) goto limits_missing;

	if ((op->haveMaxVal) && (op->haveMinVal))
		{ if (op->minVal > op->maxVal) goto conflicting_limits; }

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_erase));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_minimum:
	fprintf (stderr, "[%s] minimum limit specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

more_than_one_maximum:
	fprintf (stderr, "[%s] maximum limit specified more than once (at \"%s\")\n",
	                 name, arg);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

limits_missing:
	fprintf (stderr, "[%s] neither minimum nor maximum limit was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

conflicting_limits:
	fprintf (stderr, "[%s] conflicting limits (" valtypeFmt ">" valtypeFmt ")\n",
	                 name, op->minVal, op->maxVal);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_erase_free--

void op_erase_free (dspop* _op)
	{
	dspop_erase*	op = (dspop_erase*) _op;

	if (op->minValVarName != NULL) free (op->minValVarName);
	if (op->maxValVarName != NULL) free (op->maxValVarName);
	free (op);
	}


// op_erase_apply--

void op_erase_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_erase* op = (dspop_erase*) _op;
	valtype		minVal = op->minVal;
	valtype		maxVal = op->maxVal;
	int			ok;
	u32			ix;

	// if either limit is a named variable, fetch it now;  note that we copy
	// the value from the named variable, then destroy our reference to the
	// named variable

	if (op->minValVarName != NULL)
		{
		ok = named_global_exists (op->minValVarName, &minVal);
		if (!ok) goto no_minimum;
		op->minVal = minVal;
		fprintf (stderr, "[%s] using %s = " valtypeFmt " as minimum limit\n",
		                 _op->name, op->minValVarName, minVal);
		free (op->minValVarName);
		op->minValVarName = NULL;
		}

	if (op->maxValVarName != NULL)
		{
		ok = named_global_exists (op->maxValVarName, &maxVal);
		if (!ok) goto no_maximum;
		op->maxVal = maxVal;
		fprintf (stderr, "[%s] using %s = " valtypeFmt " as maximum limit\n",
		                 _op->name, op->maxValVarName, maxVal);
		free (op->maxValVarName);
		op->maxValVarName = NULL;
		}

	// apply limits over the vector

	if (op->keepInside)
		{
		// keep the inside, erase the outside

		if (!op->haveMaxVal)
			{ // erase below minimum only
			for (ix=0 ; ix<vLen ; ix++)
				{ if (v[ix] < minVal) v[ix] = op->zeroVal; }
			}
		else if (!op->haveMinVal)
			{ // erase above maximum only
			for (ix=0 ; ix<vLen ; ix++)
				{ if (v[ix] > maxVal) v[ix] = op->zeroVal; }
			}
		else
			{ // erase outside of minimum and maximum
			for (ix=0 ; ix<vLen ; ix++)
				{ if ((v[ix] < minVal) || (v[ix] > maxVal)) v[ix] = op->zeroVal; }
			}
		}
	else
		{
		// keep the outside, erase the inside

		if (!op->haveMaxVal)
			{ // erase above minimum only
			for (ix=0 ; ix<vLen ; ix++)
				{ if (v[ix] >= minVal) v[ix] = op->zeroVal; }
			}
		else if (!op->haveMinVal)
			{ // erase below maximum only
			for (ix=0 ; ix<vLen ; ix++)
				{ if (v[ix] <= maxVal) v[ix] = op->zeroVal; }
			}
		else
			{ // erase inside of minimum and maximum
			for (ix=0 ; ix<vLen ; ix++)
				{ if ((v[ix] >= minVal) && (v[ix] <= maxVal)) v[ix] = op->zeroVal; }
			}
		}

	// success

	return;

	// failure

no_minimum:
	fprintf (stderr, "[%s] attempt to use %s as minimum failed (no such variable)\n",
	                 _op->name, op->minValVarName);
	exit(EXIT_FAILURE);

no_maximum:
	fprintf (stderr, "[%s] attempt to use %s as maximum failed (no such variable)\n",
	                 _op->name, op->minValVarName);
	exit(EXIT_FAILURE);
	}


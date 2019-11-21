// multiply.c-- genodsp operators performing interval multiplication

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
#include "add.h"

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_multiply--
//	Multiply an incoming set of interval values (read from a file) by the
//	current set.
//
//----------

// private dspop subtype

typedef struct dspop_multiply
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			originOne;
	int			debug;
	} dspop_multiply;


// op_multiply_short--

void op_multiply_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "multiply an incoming set of interval values by the current set\n");
	}


// op_multiply_usage--

void op_multiply_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sMultiply an incoming set of interval values by the current set. Note that\n",     indent);
	fprintf (f, "%sintervals absent in the incoming set are considered to be zeros.\n",              indent);
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


// op_multiply_parse--

dspop* op_multiply_parse (char* name, int _argc, char** _argv)
	{
	dspop_multiply*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_multiply*) malloc (sizeof(dspop_multiply));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename  = NULL;
	op->valColumn = (int) get_named_global ("valColumn", 4-1);
	op->originOne = (int) get_named_global ("originOne", false);
	op->debug     = false;

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
	                 name, (int) sizeof(dspop_multiply));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_multiply_free--

void op_multiply_free (dspop* _op)
	{
	dspop_multiply*	op = (dspop_multiply*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_multiply_apply--

void op_multiply_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_multiply*	op = (dspop_multiply*) _op;
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
				if ((op->debug) && (prevEnd < chromSpec->length))
					fprintf (stderr, "zeroing %s %u..%u\n", prevChrom, prevEnd, chromSpec->length-1);

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
			v[ix] = 0.0;

		// multiply over this interval

		if (op->debug)
			fprintf (stderr, "multiplying by " valtypeFmt " over %s %u..%u\n",
			                 val, chrom, adjStart, adjEnd-1);

		for (ix=adjStart ; ix<adjEnd ; ix++)
			v[ix] *= val;

		prevEnd = adjEnd;
		}

	// clear any gap at the end of the last chromosome observed

	if (v != NULL)
		{
		if ((op->debug) && (prevEnd < chromSpec->length))
			fprintf (stderr, "zeroing %s %u..%u\n", prevChrom, prevEnd, chromSpec->length-1);

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

		if (op->debug)
			fprintf (stderr, "zeroing (all of) %s %u..%u\n",
			                 chrom, chromSpec->start+start, chromSpec->start+end-1);

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

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_divide--
//	Divide the current set of intervals by an incoming set (read from a file).
//
//----------

// private dspop subtype

typedef struct dspop_divide
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			originOne;
	valtype		infinityVal;
	int			debug;
	} dspop_divide;


// op_divide_short--

void op_divide_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "divide the current set of intervals by an incoming set\n");
	}


// op_divide_usage--

void op_divide_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sDivide the current set of intervals by an incoming set. Note that\n",             indent);
	fprintf (f, "%sintervals absent in the incoming set are considered to be zeros.\n",              indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --value=<col>            input intervals contain a value in the specified\n",   indent);
	fprintf (f, "%s                           column\n",                                             indent);
	fprintf (f, "%s  --novalue                input intervals have no value (value given is 1)\n",   indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	fprintf (f, "%s  --infinity=<value>       value to fill where denominator is zero\n",            indent);
	fprintf (f, "%s                           (default is +inf)\n",                                  indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%sThe input intervals are required to be sorted along each chromosome, and\n",      indent);
	fprintf (f, "%snon-overlapping)\n",                                                              indent);
	}


// op_divide_parse--

dspop* op_divide_parse (char* name, int _argc, char** _argv)
	{
	dspop_divide*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;
	valtype		tempVal;

	// allocate and initialize our control record

	op = (dspop_divide*) malloc (sizeof(dspop_divide));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename    = NULL;
	op->valColumn   = (int) get_named_global ("valColumn", 4-1);
	op->originOne   = (int) get_named_global ("originOne", false);
	op->infinityVal = valtypeMax;
	op->debug       = false;

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
	                 name, (int) sizeof(dspop_divide));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_divide_free--

void op_divide_free (dspop* _op)
	{
	dspop_divide*	op = (dspop_divide*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_divide_apply--

void op_divide_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_divide*	op = (dspop_divide*) _op;
	char*			filename    = op->filename;
	valtype			infinityVal = op->infinityVal;
	FILE*			f;
	char			lineBuffer[1001];
	char			prevChrom[1001];
	valtype*		v = NULL;
	char*			chrom;
	spec*			chromSpec;
	u32				start, end, o, prevEnd, adjStart, adjEnd;
	valtype			val;
	u32				ix, chromIx;
	int				ok;

	f = fopen (filename, "rt");
	if (f == NULL) goto cant_open_file;

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		chromSpec->flag = false;
		}

	// read intervals and values and divide by them

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
				if ((op->debug) && (prevEnd < chromSpec->length))
					fprintf (stderr, "infinitizing %s %u..%u\n", prevChrom, prevEnd, chromSpec->length-1);

				for (ix=prevEnd ; ix<chromSpec->length ; ix++)
					v[ix] = (v[ix]>=0)? infinityVal : -infinityVal;
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

		// infinitize any gap before this interval

		if ((op->debug) && (prevEnd < adjStart))
			fprintf (stderr, "infinitizing %s %u..%u\n", chrom, prevEnd, adjStart-1);

		for (ix=prevEnd ; ix<adjStart ; ix++)
			v[ix] = (v[ix]>=0)? infinityVal : -infinityVal;

		// divide over this interval

		if (op->debug)
			fprintf (stderr, "dividing by " valtypeFmt " over %s %u..%u\n",
			                 val, chrom, adjStart, adjEnd-1);

		for (ix=adjStart ; ix<adjEnd ; ix++)
			v[ix] /= val;  // nota bene: val is never zero

		prevEnd = adjEnd;
		}

	// infinitize any gap at the end of the last chromosome observed

	if (v != NULL)
		{
		if ((op->debug) && (prevEnd < chromSpec->length))
			fprintf (stderr, "infinitizing %s %u..%u\n", prevChrom, prevEnd, chromSpec->length-1);

		for (ix=prevEnd ; ix<chromSpec->length ; ix++)
			v[ix] = (v[ix]>=0)? infinityVal : -infinityVal;
		}

	// infinitize any chromosomes that weren't observed in the incoming set

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
			fprintf (stderr, "infinitizing (all of) %s %u..%u\n",
			                 chrom, chromSpec->start+start, chromSpec->start+end-1);

		for (ix=start ; ix<end ; ix++)
			v[ix] = (v[ix]>=0)? infinityVal : -infinityVal;
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


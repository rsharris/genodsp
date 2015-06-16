// add.c-- genodsp operators performing interval addition

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
// op_add--
//	Add an incoming set of interval values (read from a file) to the current
//	set.
//
//----------

// private dspop subtype

typedef struct dspop_add
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			originOne;
	int			destroyFile;
	} dspop_add;


// op_add_short--

void op_add_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "add an incoming set of interval values to the current set\n");
	}


// op_add_usage--

void op_add_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sAdd an incoming set of interval values to the current set.\n",                    indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --value=<col>            input intervals contain a value in the specified\n",   indent);
	fprintf (f, "%s                           column\n",                                             indent);
	fprintf (f, "%s  --novalue                input intervals have no value (value given is 1)\n",   indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	fprintf (f, "%s  --destroy                destroy the file after reading it\n",                  indent);
	fprintf (f, "%s                           (BE SURE THAT'S WHAT YOU WANT, IT CAN'T BE UNDONE)\n", indent);
	}


// op_add_parse--

dspop* op_add_parse (char* name, int _argc, char** _argv)
	{
	dspop_add*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_add*) malloc (sizeof(dspop_add));
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
	                 name, (int) sizeof(dspop_add));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_add_free--

void op_add_free (dspop* _op)
	{
	dspop_add*	op = (dspop_add*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_add_apply--

void op_add_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_add*	op = (dspop_add*) _op;
	char*		filename = op->filename;
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
		if (val == 0.0) continue;

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

		// add over this interval

		for (ix=adjStart ; ix<adjEnd ; ix++)
			v[ix] += val;
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
// op_subtract--
//	Subtract an incoming set of interval values (read from a file) from the
//	current set.
//
//----------

// private dspop subtype

typedef struct dspop_subtract
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			originOne;
	int			destroyFile;
	} dspop_subtract;


// op_subtract_short--

void op_subtract_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "subtract an incoming set of interval values from the current set\n");
	}


// op_subtract_usage--

void op_subtract_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sSubtract an incoming set of interval values from the current set.\n",                    indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --value=<col>            input intervals contain a value in the specified\n",   indent);
	fprintf (f, "%s                           column\n",                                             indent);
	fprintf (f, "%s  --novalue                input intervals have no value (value given is 1)\n",   indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	fprintf (f, "%s  --destroy                destroy the file after reading it\n",                  indent);
	fprintf (f, "%s                           (BE SURE THAT'S WHAT YOU WANT, IT CAN'T BE UNDONE)\n", indent);
	}


// op_subtract_parse--

dspop* op_subtract_parse (char* name, int _argc, char** _argv)
	{
	dspop_subtract*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			tempInt;

	// allocate and initialize our control record

	op = (dspop_subtract*) malloc (sizeof(dspop_subtract));
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
	                 name, (int) sizeof(dspop_subtract));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_subtract_free--

void op_subtract_free (dspop* _op)
	{
	dspop_subtract*	op = (dspop_subtract*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_subtract_apply--

void op_subtract_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_subtract*	op = (dspop_subtract*) _op;
	char*			filename = op->filename;
	FILE*			f;
	char			lineBuffer[1000];
	char			prevChrom[1000];
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

	// read intervals and values and subtract them

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

		// subtract over this interval

		for (ix=adjStart ; ix<adjEnd ; ix++)
			v[ix] -= val;
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
// [[-- a dsp operation function group, operating on a single chromosome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_add_constant--
//	Add a constant value to the current set of interval values.
//
//----------

// private dspop subtype

typedef struct dspop_addconst
	{
	dspop		common;			// common elements shared with all operators
	valtype		val;
	} dspop_addconst;

// op_add_constant_short--

void op_add_constant_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "add a constant value to the current set of interval values\n");
	}


// op_add_constant_usage--

void op_add_constant_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sAdd a constant value to the current set of interval values.\n",                   indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <value>\n", indent, name);
	}


// op_add_constant_parse--

dspop* op_add_constant_parse (char* name, int _argc, char** _argv)
	{
	dspop_addconst*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			haveVal;

	// allocate and initialize our control record

	op = (dspop_addconst*) malloc (sizeof(dspop_addconst));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	// parse arguments

	haveVal = false;

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <value>

		if (!haveVal)
			{
			op->val = string_to_valtype (arg);
			haveVal = true;
			goto next_arg;
			}

		// unknown argument

		chastise ("[%s] Can't understand \"%s\"\n", name, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	if (!haveVal) goto constant_missing;

	return (dspop*) op;

cant_allocate:
	fprintf (stderr, "[%s] failed to allocate control record (%d bytes)\n",
	                 name, (int) sizeof(dspop_addconst));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

constant_missing:
	fprintf (stderr, "[%s] no constant value was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_add_constant_free--

void op_add_constant_free (dspop* op)
	{
	free (op);
	}


// op_add_constant_apply--

void op_add_constant_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_addconst*	op = (dspop_addconst*) _op;
	valtype			val = op->val;
	u32				ix;

	if (val == 0.0) return;

	for (ix=0 ; ix<vLen ; ix++)
		v[ix] += val;

	}

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_invert--
//	Apply an inversion filter;  arithmetic inversion preserving min and max.
//
//----------

// private dspop subtype

typedef struct dspop_invert
	{
	dspop		common;			// common elements shared with all operators
	int			haveMidVal;
	valtype		midVal;
	} dspop_invert;

// op_invert_short--

void op_invert_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "apply an arithmetic inversion filter, preserving min and max\n");
	}


// op_invert_usage--

void op_invert_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sApply an inversion filter;  arithmetic inversion preserving min and max.\n",      indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [<value>]\n", indent, name);
	fprintf (f, "%s  <value>                  'middle value' to reflect values about\n",             indent);
	fprintf (f, "%s                           (by default we use (min+max)/2)\n",                    indent);
	}


// op_invert_parse--

dspop* op_invert_parse (char* name, int _argc, char** _argv)
	{
	dspop_invert*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;
	int			haveMidVal;

	// allocate and initialize our control record

	op = (dspop_invert*) malloc (sizeof(dspop_invert));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->haveMidVal = false;
	op->midVal     = 0.0;

	// parse arguments

	haveMidVal = false;

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("[%s] Can't understand \"%s\"\n", name, arg);

		// <value> and variants

		if ((!haveMidVal)
		 && ((strcmp (arg, "zero") == 0) || (strcmp (arg, "negate") == 0)))
			{
			op->midVal = 0.0;
			op->haveMidVal = haveMidVal = true;
			goto next_arg;
			}

		if ((!haveMidVal)
		 && (strcmp (arg, "one") == 0))
			{
			op->midVal = 1.0;
			op->haveMidVal = haveMidVal = true;
			goto next_arg;
			}

		if ((!haveMidVal)
		 && ((strcmp (arg, "1/2") == 0) || (strcmp (arg, "binary") == 0)))
			{
			op->midVal = 0.5;
			op->haveMidVal = haveMidVal = true;
			goto next_arg;
			}

		if (!haveMidVal)
			{
			op->midVal = string_to_valtype (arg);
			op->haveMidVal = haveMidVal = true;
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
	                 name, (int) sizeof(dspop_invert));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_invert_free--

void op_invert_free (dspop* op)
	{
	free (op);
	}


// op_invert_apply--

void op_invert_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_invert*	op = (dspop_invert*) _op;
	valtype			midVal, minVal, maxVal;
	spec*			chromSpec;
	u32				ix, chromIx;

	// if the user didn't specify a middle value, derive one from the min and
	// max;  the derived value is such that the min and max will be preserved
	// in the output

	if (op->haveMidVal)
		midVal = op->midVal;
	else
		{
		chromSpec = chromsSorted[0];
		v = chromSpec->valVector;
		minVal = maxVal = v[0];

		for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
			{
			chromSpec = chromsSorted[chromIx];
			v = chromSpec->valVector;

			for (ix=0 ; ix<chromSpec->length ; ix++)
				{
				if (v[ix] < minVal) minVal = v[ix];
				if (v[ix] > maxVal) maxVal = v[ix];
				}
			}

		midVal = (minVal + maxVal) / 2.0;
		}

	// perform the inversion

	for (chromIx=0 ; chromsSorted[chromIx]!=NULL ; chromIx++)
		{
		chromSpec = chromsSorted[chromIx];
		v = chromSpec->valVector;

		for (ix=0 ; ix<chromSpec->length ; ix++)
			v[ix] = 2*midVal - v[ix];
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
// op_absolute_value--
//	Compute the absolute value of the current set of interval values.
//
//----------

// op_absolute_value_short--

void op_absolute_value_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "compute the absolute value of the current set of interval values\n");
	}


// op_absolute_value_usage--

void op_absolute_value_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sCompute the absolute value of the current set of interval values.\n",             indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s\n", indent, name);
	}


// op_absolute_value_parse--

dspop* op_absolute_value_parse (char* name, int _argc, char** _argv)
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


// op_absolute_value_free--

void op_absolute_value_free (dspop* op)
	{
	free (op);
	}


// op_absolute_value_apply--

void op_absolute_value_apply
   (arg_dont_complain(dspop*	op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	u32		ix;

	for (ix=0 ; ix<vLen ; ix++)
		{ if (v[ix] < 0) v[ix] = -v[ix]; }

	}


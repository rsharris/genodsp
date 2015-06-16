// opio.c-- genodsp operators performing input/output

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
#include "opio.h"

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_input--
//	Read intervals from a file.
//
//----------

// private dspop subtype

typedef struct dspop_input
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			valColumn;
	int			missingVal;
	int			overlapOp;
	int			originOne;
	int			destroyFile;
	} dspop_input;

// op_input_short--

void op_input_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "read intervals from a file (replacing the current set)\n");
	}


// op_input_usage--

void op_input_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sRead intervals from a file (replacing the current set).\n",                       indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --value=<col>            input intervals contain a value in the specified\n",   indent);
	fprintf (f, "%s                           column\n",                                             indent);
	fprintf (f, "%s  --novalue                input intervals have no value (value given is 1)\n",   indent);
	fprintf (f, "%s  --missing=<value>        value for intervals missing from input\n",             indent);
	fprintf (f, "%s                           (by default, this is zero)\n",                         indent);
	fprintf (f, "%s  --overlap=sum            when a position is in multiple intervals, sum the\n",  indent);
	fprintf (f, "%s                           values\n",                                             indent);
	fprintf (f, "%s                           (this is the default)\n",                              indent);
	fprintf (f, "%s  --overlap=minimum        when a position is in multiple intervals, keep the\n", indent);
	fprintf (f, "%s                           minimum value\n",                                      indent);
	fprintf (f, "%s  --overlap=maximum        when a position is in multiple intervals, keep the\n", indent);
	fprintf (f, "%s                           maximum value\n",                                      indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	fprintf (f, "%s  --destroy                destroy the file after reading it\n",                  indent);
	fprintf (f, "%s                           (BE SURE THAT'S WHAT YOU WANT, IT CAN'T BE UNDONE)\n", indent);
	}


// op_input_parse--

dspop* op_input_parse (char* name, int _argc, char** _argv)
	{
	dspop_input*	op;
	int				argc = _argc;
	char**			argv = _argv;
	char*			arg, *argVal;
	valtype			tempVal;
	int				tempInt;

	// allocate and initialize our control record

	op = (dspop_input*) malloc (sizeof(dspop_input));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename    = NULL;
	op->valColumn   = (int) get_named_global ("valColumn", 4-1);
	op->missingVal  = 0.0;
	op->overlapOp   = ri_overlapSum;
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

		// --missing=<value>

		if (strcmp_prefix (arg, "--missing=") == 0)
			{
			tempVal = string_to_valtype (argVal);
			op->missingVal = tempVal;
			goto next_arg;
			}

		// --overlap=sum, etc.

		if (strcmp (arg, "--overlap=sum") == 0)
			{ op->overlapOp = ri_overlapSum;  goto next_arg; }

		if ((strcmp (arg, "--overlap=minimum") == 0)
		 || (strcmp (arg, "--overlap=min"    ) == 0))
			{ op->overlapOp = ri_overlapMin;  goto next_arg; }

		if ((strcmp (arg, "--overlap=maximum") == 0)
		 || (strcmp (arg, "--overlap=max")     == 0))
			{ op->overlapOp = ri_overlapMax;  goto next_arg; }

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
	                 name, (int) sizeof(dspop_input));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_input_free--

void op_input_free (dspop* _op)
	{
	dspop_input*	op = (dspop_input*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_input_apply--

void op_input_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	_v))
	{
	dspop_input*	op = (dspop_input*) _op;
	FILE*			f;

	f = fopen (op->filename, "rt");
	if (f == NULL) goto cant_open_file;

	read_intervals (f, op->valColumn, op->originOne, op->overlapOp,
	                /*clear*/ true, op->missingVal);
	fclose (f);

	if (op->destroyFile)
		remove (op->filename);

	// success

	return;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "[%s] can't open \"%s\" for reading\n",
	                 op->common.name, op->filename);
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
// op_output--
//	Write intervals to a file.
//
//----------

// private dspop subtype

typedef struct dspop_output
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			noOutputValues;
	int			valPrecision;
	int			collapseRuns;
	int			showUncovered;
	int			originOne;
	} dspop_output;


// op_output_short--

void op_output_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "write the current set of intervals to a file\n");
	}


// op_output_usage--

void op_output_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sWrite the current set of intervals to a file.\n",                                 indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --nooutputvalue          don't write value with output intervals\n",            indent);
	fprintf (f, "%s  --precision=<number>     number of digits to round output values to\n",         indent);
	fprintf (f, "%s  --nocollapse             in output, don't collapse runs of identical values\n", indent);
	fprintf (f, "%s                           to intervals\n",                                       indent);
	fprintf (f, "%s  --uncovered:hide         don't output intervals that have no coverage\n",       indent);
	fprintf (f, "%s  --uncovered:show         in output, include intervals that have no coverage\n", indent);
	fprintf (f, "%s  --uncovered:NA           in output, mark uncovered intervals as NA\n",          indent);
	fprintf (f, "%s  --origin=one             input/output intervals are origin-one, closed\n",      indent);
	fprintf (f, "%s  --origin=zero            input/output intervals are origin-zero, half-open\n",  indent);
	}


// op_output_parse--

dspop* op_output_parse (char* name, int _argc, char** _argv)
	{
	dspop_output*	op;
	int				argc = _argc;
	char**			argv = _argv;
	char*			arg, *argVal;
	int				tempInt;

	// allocate and initialize our control record

	op = (dspop_output*) malloc (sizeof(dspop_output));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = true;

	op->filename       = NULL;
	op->noOutputValues = (int) get_named_global ("noOutputValues", false);
	op->valPrecision   = (int) get_named_global ("valPrecision",   0);
	op->collapseRuns   = (int) get_named_global ("collapseRuns",   true);
	op->showUncovered  = (int) get_named_global ("showUncovered",  uncovered_hide);
	op->originOne      = (int) get_named_global ("originOne",      false);

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --precision=<col>

		if (strcmp_prefix (arg, "--precision=") == 0)
			{
			tempInt = string_to_int (argVal);
			if (tempInt < 0)
				chastise ("[%s] precision can't be negative (\"%s\")\n", name, arg);
			op->valPrecision = tempInt;
			goto next_arg;
			}

		if ((strcmp (arg, "--nooutputvalue")  == 0)
		 || (strcmp (arg, "--nooutputvalues") == 0))
			{ op->noOutputValues = true; goto next_arg; }

		// --nocollapse

		if (strcmp (arg, "--nocollapse") == 0)
			{ op->collapseRuns = false;  goto next_arg; }

		// --uncovered:hide

		if ((strcmp (arg, "--uncovered:hide") == 0)
		 || (strcmp (arg, "--hide:uncovered") == 0))
			{ op->showUncovered = uncovered_hide;  goto next_arg; }

		// --uncovered:show

		if ((strcmp (arg, "--uncovered:show") == 0)
		 || (strcmp (arg, "--show:uncovered") == 0))
			{ op->showUncovered = uncovered_show;  goto next_arg; }

		// --uncovered:NA

		if ((strcmp (arg, "--uncovered:NA") == 0)
		 || (strcmp (arg, "--uncovered:mark") == 0)
		 || (strcmp (arg, "--mark:uncovered") == 0)
		 || (strcmp (arg, "--markgaps")       == 0))
			{ op->showUncovered = uncovered_NA;  goto next_arg; }

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
	                 name, (int) sizeof(dspop_output));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_output_free--

void op_output_free (dspop* _op)
	{
	dspop_output*	op = (dspop_output*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_output_apply--

void op_output_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_output*	op = (dspop_output*) _op;
	FILE*			f;

	f = fopen (op->filename, "wt");
	if (f == NULL) goto cant_open_file;

	report_intervals (f, op->valPrecision, op->noOutputValues,
	                  op->collapseRuns, op->showUncovered, op->originOne);

	// success

	fclose (f);
	return;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "[%s] can't open \"%s\" for writing\n",
	                 op->common.name, op->filename);
	exit (EXIT_FAILURE);
	}


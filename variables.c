// variables.c-- genodsp 'operator' performing operations on global variables

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
#include "variables.h"

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_show_variables--
//	This is not really a dsp operation, as it doesn't perform any operation
//	on the vectors.  Instead, it just provides some functions realting to the
//	global variable state.
//
//----------

// op_show_variables_short--

void op_show_variables_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "inspect the named variables state\n");
	}


// op_show_variables_usage--

void op_show_variables_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sInspect the named variables state\n",                                            indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s [options]\n", indent, name);
	fprintf (f, "%s  (no options yet)\n",                                                            indent);
	}

// op_show_variables_parse--

dspop* op_show_variables_parse (char* name, int _argc, char** _argv)
	{
	dspop*	op;
	int		argc = _argc;
	char**	argv = _argv;
	char*	arg, *argVal;

	// allocate and initialize our control record

	op = (dspop*) malloc (sizeof(dspop));
	if (op == NULL) goto cant_allocate;

	op->atRandom = true;

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


// op_show_variables_free--

void op_show_variables_free (dspop* op)
	{
	free (op);
	}


// op_show_variables_apply--

void op_show_variables_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	fprintf (stderr, "variables:\n");
	report_named_globals (stderr, "  ");
	}

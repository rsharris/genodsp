// map.c-- genodsp operator performing general signal mapping

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
#include "map.h"

// piecewise mapping 

typedef struct mappingel
	{
	valtype	vIn;				// function input value
	valtype	vOut;				// function output value
	} mappingel;

typedef struct mapping
	{
	// note that v is allocated as part of the same heap block
	u32		len;				// number of elements in vIn[] and vOut[]
	mappingel* v;				// function input-to-output values
	} mapping;

// prototypes for private functions

static mapping* read_mapping    (FILE* f);
static int      read_value_pair (FILE* f, char* buffer, int bufferLen,
                                 valtype* val1, valtype* val2);

// miscellany

#define noIndex ((u32) -1)

//----------
// [[-- a dsp operation function group, operating on the whole genome --]]
//
// See genodsp_interface.h, "headers for dsp operator function groups" for
// function descriptions and argument details.
//
//----------
//
// op_map--
//	Map signal values per a piecewise-linear function.
//
//----------

// private dspop subtype

typedef struct dspop_map
	{
	dspop		common;			// common elements shared with all operators
	char*		filename;
	int			destroyFile;
	int			debug;
	} dspop_map;


// op_map_short--

void op_map_short (char* name, int nameWidth, FILE* f, char* indent)
	{
	int nameFill = nameWidth-2 - strlen(name);
	if (indent == NULL) indent = "";

	if (nameFill > 0) fprintf (f, "%s%s:%*s", indent, name, nameFill+1, " ");
	             else fprintf (f, "%s%s: ", indent, name);

	fprintf (f, "map values according to a piecewise-linear function\n");
	}


// op_map_usage--

void op_map_usage (char* name, FILE* f, char* indent)
	{
	if (indent == NULL) indent = "";
	//             3456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (f, "%sMap values according to a piecewise-linear function.\n",                          indent);
	fprintf (f, "%s\n", indent);
	fprintf (f, "%susage: %s <filename> [options]\n", indent, name);
	fprintf (f, "%s  --destroy                destroy the file after reading it\n",                  indent);
	fprintf (f, "%s                           (BE SURE THAT'S WHAT YOU WANT, IT CAN'T BE UNDONE)\n", indent);
	fprintf (f, "%s\n",                                                                              indent);
	fprintf (f, "%sThe file describes a piecewise-linear function. It is a text file consisting\n",  indent);
	fprintf (f, "%sof two values per line. The first value is an input value, the second is the\n",  indent);
	fprintf (f, "%scorresponding output value. The order of lines is unimportant, but usually\n",    indent);
	fprintf (f, "%sthey are sorted by input values.\n",                                              indent);
	fprintf (f, "%s\n",                                                                              indent);
	fprintf (f, "%sThe function maps given input values to the corresponding output value. Input\n", indent);
	fprintf (f, "%svalues not appearing in the file are interpolated linearly using the two\n",      indent);
	fprintf (f, "%sneighboring input values. Input values outside the range in the file are\n",      indent);
	fprintf (f, "%smapped to the output value for the minimum or maximum input value.\n",            indent);
	}


// op_map_parse--

dspop* op_map_parse (char* name, int _argc, char** _argv)
	{
	dspop_map*	op;
	int			argc = _argc;
	char**		argv = _argv;
	char*		arg, *argVal;

	// allocate and initialize our control record

	op = (dspop_map*) malloc (sizeof(dspop_map));
	if (op == NULL) goto cant_allocate;

	op->common.atRandom = false;

	op->filename    = NULL;
	op->destroyFile = false;
	op->debug       = false;

	// parse arguments

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --value=<col>

		// --destroy

		if (strcmp (arg, "--destroy") == 0)
			{ op->destroyFile = true;  goto next_arg; }

		// --debug

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
	                 name, (int) sizeof(dspop_map));
	exit(EXIT_FAILURE);

filename_missing:
	fprintf (stderr, "[%s] no filename was provided\n",
	                 name);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


// op_map_free--

void op_map_free (dspop* _op)
	{
	dspop_map*	op = (dspop_map*) _op;

	if (op->filename != NULL) free (op->filename);
	free (op);
	}


// op_map_apply--

void op_map_apply
   (arg_dont_complain(dspop*	_op),
	arg_dont_complain(char*		vName),
	arg_dont_complain(u32		vLen),
	arg_dont_complain(valtype*	v))
	{
	dspop_map*	op = (dspop_map*) _op;
	mapping*	map = NULL;
	FILE*		f;
	valtype		minIn, maxIn, outForMin, outForMax;
	valtype		pieceLo, pieceHi, outForLo, outForHi, inDiff, outDiff;
	u32			maxIx, pieceIx, ixLo, ixHi, ixMid;
	u32			ix;
	valtype		inVal;

	// read the mapping file

	f = fopen (op->filename, "rt");
	if (f == NULL) goto cant_open_file;

	map = read_mapping (f);
	fclose (f);
	if (map == NULL) goto cat_read_mapping;

	if (op->destroyFile)
		remove (op->filename);

	if (op->debug)
		{
		fprintf (stderr, "mapping:\n");
		for (pieceIx=0 ; pieceIx<map->len ; pieceIx++)
			fprintf (stderr, "  [%u] " valtypeFmt " -> " valtypeFmt "\n",
			                 pieceIx, map->v[pieceIx].vIn, map->v[pieceIx].vOut);
		}

	maxIx     = map->len - 1;
	minIn     = map->v[0].vIn;
	maxIn     = map->v[maxIx].vIn;
	outForMin = map->v[0].vOut;
	outForMax = map->v[maxIx].vOut;

	// apply the mapping

	pieceIx  = noIndex;
	pieceLo  = pieceHi  = 0.0;	// (placate compiler)
	outForLo = outForHi = 0.0;	// (placate compiler)
	inDiff   = outDiff  = 0.0;	// (placate compiler)
	ixLo     = ixHi     = 0;	// (placate compiler)

	for (ix=0 ; ix<vLen ; ix++)
		{
		if (op->debug)
			{
			fprintf (stderr, "[%u] " valtypeFmt "\n", ix, v[ix]);
			if (v[ix] <= minIn)
				{
				fprintf (stderr, "  below minimum\n");
				fprintf (stderr, "  --> " valtypeFmt "\n", outForMin);
				}
			if (v[ix] >= maxIn)
				{
				fprintf (stderr, "  above minimum\n");
				fprintf (stderr, "  --> " valtypeFmt "\n", outForMax);
				}
			}

		// if the next unmapped value is beyond the range of the mapping,
		// assign it the edge value

		inVal = v[ix];
		if (inVal <= minIn) { v[ix] = outForMin;  continue; }
		if (inVal >= maxIn) { v[ix] = outForMax;  continue; }

		// otherwise, find the piece the unmapped value is on;  if it is on
		// the most recent piece we've used, or a piece next to it, we can
		// shortcut the search;  otherwise, we do a binary search over the
		// possible pieces (full range the first time, limited to above or
		// below the most recent piece aferwards);  the rational behind the
		// shortcut operations is that most signals will vary smoothly, and
		// we want to avoid doing the binary search if we can
		
		if (pieceIx == noIndex)
			{ ixLo = 0;  ixHi = maxIx;  goto find_piece; }

		if ((inVal >= pieceLo) && (inVal <= pieceHi))
			{
			if (op->debug)
				fprintf (stderr, "  on same piece\n");
			goto piece_known;
			}

		if (inVal < pieceLo)
			{
			if ((pieceIx == 0) || (inVal < map->v[pieceIx-1].vIn))
				{ ixLo = 0;  ixHi = pieceIx-1;  goto find_piece; }
			ixLo = pieceIx-1;
			if (op->debug)
				fprintf (stderr, "  on next piece below\n");
			goto piece_found;
			}

		if (inVal > pieceHi)
			{
			if ((pieceIx > maxIx-2) || (inVal > map->v[pieceIx+2].vIn))
				{ ixLo = pieceIx+2;  ixHi = maxIx;  goto find_piece; }
			ixLo = pieceIx+1;
			if (op->debug)
				fprintf (stderr, "  on next piece above\n");
			goto piece_found;
			}

	find_piece:
		if (op->debug)
			fprintf (stderr, "  on new piece\n");

		while (ixLo+1 < ixHi)
			{
			// invariant: v[ixLo].vIn <= inVal < v[ixHi].vIn

			if (op->debug)
				fprintf (stderr, "  %u..%u " valtypeFmt ".." valtypeFmt "\n",
								 ixLo,ixHi, map->v[ixLo].vIn, map->v[ixHi].vIn);

			ixMid = (ixLo + ixHi) / 2;
			if      (inVal < map->v[ixMid].vIn) ixHi = ixMid;
			else if (inVal > map->v[ixMid].vIn) ixLo = ixMid;
			else                              { ixLo = ixMid;  break; }
			}

		while ((ixLo < maxIx) && (map->v[ixLo].vIn == map->v[ixLo+1].vIn))
			ixLo++;

	piece_found:
		if (op->debug)
			fprintf (stderr, "  piece %u"
							 " maps " valtypeFmt ".." valtypeFmt
							 " to "   valtypeFmt ".." valtypeFmt "\n",
							 ixLo,
							 map->v[ixLo].vIn,  map->v[ixLo+1].vIn,
							 map->v[ixLo].vOut, map->v[ixLo+1].vOut);

		pieceIx  = ixLo;
		pieceLo  = map->v[pieceIx  ].vIn;
		pieceHi  = map->v[pieceIx+1].vIn;
		outForLo = map->v[pieceIx  ].vOut;
		outForHi = map->v[pieceIx+1].vOut;
		inDiff   = pieceHi  - pieceLo;
		outDiff  = outForHi - outForLo;

		// map the value on the piece

	piece_known:
		if (inVal == pieceLo)
			v[ix] = outForLo;
		else if (inVal == pieceHi)
			v[ix] = outForHi;
		else
			{
			v[ix] = outForLo + (inVal - pieceLo) * outDiff / inDiff;
			if (op->debug)
				fprintf (stderr, "  " valtypeFmt
				                 " maps to " valtypeFmt " + (" valtypeFmt " * " valtypeFmt ")\n",
				                 inVal, outForLo, inVal-pieceLo, outDiff/inDiff);
			}

		if (op->debug)
			fprintf (stderr, "  --> " valtypeFmt "\n", v[ix]);
		}

	// success

	free (map);
	return;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "[%s] can't open \"%s\" for reading\n",
	                 op->common.name, op->filename);
	exit (EXIT_FAILURE);

cat_read_mapping:
	fprintf (stderr, "[%s] problem with mapping file \"%s\"\n",
	                 op->common.name, op->filename);
	exit (EXIT_FAILURE);
	}

//----------
//
// read_mapping--
//	Read a piecewise-linear function description from a file.
//
// The file is a text file consisting of two values per line.  The first value
// is an input value, the second is the corresponding output value.  Order of
// lines is unimportant, but usually the lines are sorted by input values.
//
//----------
//
// Arguments:
//	FILE*	f:	The file to read from.  This should have been opened as "rt",
//				and we assume it is rewindable.
//
// Returns:
//	A pointer to the mapping structure;  NULL if there were any problems.
//
//----------

static int mapping_element_ascending (const void* _v1, const void* _v2);
static int mapping_element_ascending (const void* _v1, const void* _v2)
	{
	const mappingel* v1 = (const mappingel*) _v1;
	const mappingel* v2 = (const mappingel*) _v2;

	return (v1->vIn > v2->vIn) - (v1->vIn < v2->vIn);
	}

// read_mapping--

static mapping* read_mapping
   (FILE*		f)
	{
	char		lineBuffer[1000];
	mapping*	map = NULL;
	valtype		vIn, vOut;
	u32			numPieces, ix;
	u32			bytesNeeded, vOffset;
	int			ok;

	// make a first pass through the file to count the number of pieces

	numPieces = 0;

	while (true)
		{
		ok = read_value_pair (f, lineBuffer, sizeof(lineBuffer), NULL, NULL);
		if (!ok) break;
		numPieces++;
		}

	// allocate the mapping structure

	bytesNeeded =  round_up_16 (sizeof(mapping));
	vOffset     =  bytesNeeded;
	bytesNeeded += numPieces * sizeof(mappingel);

	map = (mapping*) malloc (bytesNeeded);
	if (map == NULL) return NULL;

	map->len = numPieces;
	map->v   = (mappingel*) (((char*) map) + vOffset);

	// re-read the file and fill the mapping arrays

	rewind (f);

	ix = 0;
	while (true)
		{
		ok = read_value_pair (f, lineBuffer, sizeof(lineBuffer), &vIn, &vOut);
		if (!ok) break;
		map->v[ix].vIn  = vIn;
		map->v[ix].vOut = vOut;
		ix++;
		}

	// sort the pieces by increasing input value

	qsort (map->v, numPieces, sizeof(mappingel), mapping_element_ascending);

	return map;
	}

//----------
//
// read_value_pair--
//	Read the next value pair from a file.
//
//----------
//
// Arguments:
//	FILE*		f:			File to read from.
//	char*		buffer:		Buffer to read the line into.  Note that the caller
//							.. should not expect anything about the contents of
//							.. this buffer upon return.
//	int			bufferLen:	Number of bytes allocated for the buffer.
//	valtype*	val1:		Place to return the first value.
//	valtype*	val2:		Place to return the second value.
//
// Returns:
//	true if we were successful;  false if there are no more lines in the file.
//	Failures result in program termination.
//
//----------

static int read_value_pair
   (FILE*		f,
	char*		buffer,
	int			bufferLen,
	valtype*	_val1,
	valtype*	_val2)
	{
	static u32	lineNumber = 0;
	static int	missingEol = false;
	int			lineLen;
	char*		scan, *mark, *field;
	valtype		val1, val2;

	// read the next line

try_again:

	if (fgets (buffer, bufferLen, f) == NULL)
		return false;

	lineNumber++;

	// check for lines getting split by fgets (the final line in the file might
	// not have a newline, but no internal lines can be that way)

	if (missingEol) goto missing_eol;

	lineLen = strlen(buffer);
	if (lineLen != 0)
		missingEol = (buffer[lineLen-1] != '\n');

	// parse the line

	scan = skip_whitespace(buffer);
	if (*scan == 0)   goto try_again;  // empty line
	if (*scan == '#') goto try_again;  // comment line

	scan = buffer;
	if (*scan == ' ') goto no_val1;
	field = scan;
	mark = skip_darkspace(scan);
	scan = skip_whitespace(mark);
	if (*mark != 0) *mark = 0;
	val1 = string_to_valtype (field);

	if (*scan == 0) goto no_val2;
	field = scan;
	mark = skip_darkspace(scan);
	scan = skip_whitespace(mark);
	if (*mark != 0) *mark = 0;
	val2 = string_to_valtype (field);

	//////////
	// success
	//////////

	if (_val1  != NULL) *_val1   = val1;
	if (_val2  != NULL) *_val2   = val2;

	return true;

	//////////
	// failure exits
	//////////

missing_eol:
	fprintf (stderr, "problem at line %u, line is longer than internal buffer\n",
			 lineNumber-1);
	exit (EXIT_FAILURE);

no_val1:
	fprintf (stderr, "problem at line %u, line contains no first value\n",
			 lineNumber);
	exit (EXIT_FAILURE);

no_val2:
	fprintf (stderr, "problem at line %u, line contains no second value\n",
			 lineNumber);
	exit (EXIT_FAILURE);
	}

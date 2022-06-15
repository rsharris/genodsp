#ifndef genodsp_interface_H		// (prevent multiple inclusion)
#define genodsp_interface_H

// establish ownership of global variables

#ifdef globals_owner
#define global
#else
#define global extern
#endif

//----------
//
// chromosome value vectors--
//
//----------

// item values (values along the chromosomes)

typedef double valtype;
#define string_to_valtype(s)       ((valtype) string_to_double(s))
#define try_string_to_valtype(s,v) try_string_to_double(s,(valtype*)v)
#define valtypeFmt     "%f"
#define valtypeFmtPrec "%.*f"
#define valtypeMax     DBL_MAX
#define valtypePuny    DBL_MIN

// chromosomes-of-interest
//
// chromsOfInterest is a linked list of the chromosomes, in the order they were
// input;  operators that write chromosomes to a file should use this order
//
// chromsSorted is an array of the chromosomes, sorted from longest to shortest;
// operators that need to access the list of chromosomes should use this array
// instead of the linked list;  the array is terminated by a NULL entry

typedef struct spec
	{
	struct spec* next;			// next spec in a linked list
	char*		chrom;			// chromosome name
	int			flag;			// (internal use)
	u32			start;			// number of uninteresting bases at the start
								// .. of the chromosome (these are not included
								// .. in the vector)
	u32			length;			// the length of v[], i.e. the number of
								// .. interesting bases in the chromosome;
								// .. this is guaranteed to be non-zero
	valtype*	valVector;		// vector of values
	} spec;

#ifdef globals_owner
global spec*  chromsOfInterest = NULL;
global spec** chromsSorted     = NULL;
#else
global spec*  chromsOfInterest;
global spec** chromsSorted;
#endif

//----------
//
// operators--
//
//----------

// dsp operator functions
//
// each operator consists of five functions
//	short: short descriptive text (one line)
//	usage: long descriptive text
//	parse: parse command line arguments and allocate control record
//	free:  de-allocate control record
//	apply: apply function to vector(s)
//
// headers for these functions are show later in this file

#define opfuncargs_short (char*,int,FILE*,char*)
#define opfuncargs_usage (char*,FILE*,char*)
#define opfuncargs_parse (char*,int,char**)
#define opfuncargs_free  (struct dspop*)
#define opfuncargs_apply (struct dspop*,char*,u32,valtype*)

typedef void          (*opfunc_short) opfuncargs_short;
typedef void          (*opfunc_usage) opfuncargs_usage;
typedef struct dspop* (*opfunc_parse) opfuncargs_parse;
typedef void          (*opfunc_free)  opfuncargs_free;
typedef void          (*opfunc_apply) opfuncargs_apply;

#define dspprototypes(funcName) \
void          funcName##_short opfuncargs_short; \
void          funcName##_usage opfuncargs_usage;  \
struct dspop* funcName##_parse opfuncargs_parse; \
void          funcName##_free  opfuncargs_free;  \
void          funcName##_apply opfuncargs_apply;

// linked list for dsp operators
//
// the list will actually contain a mixture of records for different operators;
// each operator should define its own type (essentially a subtype of dspop), a
// struct with dspop as the first element

typedef struct dspop
	{
	struct dspop*	next;		// next operation in a linked list
	char*			name;		// operation
	opfunc_apply	funcApply;	// functions that perform the operation
	opfunc_free		funcFree;
	int				atRandom;	// true => this function needs 'random' access
								//         .. to hop around the whole genome
	} dspop;

typedef struct dspinfo
	{
	char*			name;		// operation
	opfunc_short	funcShort;	// functions that perform the operation
	opfunc_usage	funcUsage;
	opfunc_parse	funcParse;
	opfunc_free		funcFree;
	opfunc_apply	funcApply;
	} dspinfo;

#define dspinforecord(name,funcName) \
	{ name, funcName##_short, funcName##_usage, funcName##_parse, funcName##_free, funcName##_apply }

#define dspinfoalias(name) \
	{ name, NULL, NULL, NULL, NULL, NULL }

//----------
//
// miscellany--
//
//----------

#ifndef M_PI
#define M_PI 3.14159265358979323846264
#endif

// global command line options 

#ifdef globals_owner
global int trackOperations  = false;
global int reportComments   = false;
global u32 reportInputProgress = 0;
#else
global int trackOperations;
global int reportComments;
global u32 reportInputProgress;
#endif

// values for showUncovered

#define uncovered_NA   -1
#define uncovered_show 1
#define uncovered_hide 0

// values for read_interval's overlapOp

#define ri_overlapSum 0
#define ri_overlapMin 1
#define ri_overlapMax 2

//----------
//
// protypes for entries into genodsp.c--
//
//----------

void     chastise               (const char* format, ...);
spec*    find_chromosome_spec   (char* chrom);
void     read_intervals         (FILE* f, int valCol, int originOne,
                                 int overlapOp, int clear, valtype missingVal);
int      read_interval          (FILE* f,
                                 char* buffer, int bufferLen, int valCol,
                                 char** chrom, u32* start, u32* end,
                                 valtype* val);
void     report_intervals       (FILE* f,
                                 int precision,
                                 int noOutputValues, int collapseRuns,
                                 int showUncovered, int originOne);
void     read_all_chromosomes   (char* filename);
void     write_all_chromosomes  (char* filename);
valtype* get_scratch_vector     (void);
s32*     get_scratch_ints       (void);
void     release_scratch_vector (valtype* v);
void     release_scratch_ints   (s32* v);
void     set_named_global       (char* name, valtype val);
valtype  get_named_global       (char* name, valtype defaultVal);
int      named_global_exists    (char* name, valtype* val);
void     report_named_globals   (FILE* f, char* indent);
void     tracking_report        (const char* format, ...);
int      valtype_ascending      (const void* v1, const void* v2);

////////////
//
// headers for dsp operator function groups--
//
////////////

//----------
//
// op_short--
//	Print a short description of this operation, usually on one line.
//
//----------
//
// Arguments:
//	char*	name:		The name of the operator (what the caller likes to call
//						.. it).
//	int		nameWidth:	The width of the name field.  This routine has to add
//						.. padding after the name field (and colon) so that the
//						.. the description starts after this many characters
//						.. (not counting indent)
//	FILE*	f:			The stream to print to.
//	char*	indent:		A prefix to print at the start of every line.  Usually
//						.. this is a string of spaces.  This may be NULL, which
//						.. should have the same effect as an empty string.
//
// Returns:
//	(nothing)
//
//----------

//----------
//
// op_usage--
//	Print a longer description of this operation.  This should look like the
//	typical usuage text for a command-line program, describing the arguments
//	to use, etc.
//
//----------
//
// Arguments:
//	char*	name:	The name of the operator (what the caller likes to call it).
//	FILE*	f:		The stream to print to.
//	char*	indent:	A prefix to print at the start of every line.  Usually this
//					.. is a string of spaces.  This may be NULL, which should
//					.. have the same effect as an empty string.
//
// Returns:
//	(nothing)
//
//----------

//----------
//
// op_parse--
//	Parse a command line for this operation, and allocate a control record.
//
//----------
//
// Arguments:
//	char*	name:	The name of the operator (what the caller likes to call it).
//					.. This is expected to be used only for presenting error
//					.. messages to the user.
//	int		argc:	Number of arguments in argv[].
//	char**	argv:	Command line arguments.  Note that, unlike the usual argv
//					.. passed to main(), argv[0] contains an argument, *not*
//					.. the program name.
//
// Returns:
//	The operator's control record, allocated from the heap.  This is compatible
//	with a dspop record, but usually will be a type private to the operator,
//	containing additional fields.
//
//	This function must fill in one field in the common section of the record,
//	atRandom, which tells the caller whether the operator performs on multiple
//	vectors (atRandom=true) or single vectors (atRandom=false).
//
//	The operator will be responsible for de-allocating this control record, in
//	op_free().
//
//----------

//----------
//
// op_free--
//	De-allocate memory allocated for this operation.
//
//----------
//
// Arguments:
//	dspop*		op:		Pointer to the operator's control record to be de-
//						.. allocated (the record was allocated by op_parse).
//						.. This is compatible with a dspop record, but often
//						.. this function will cast it to the type of its
//						.. private control record.  This routine must free the
//						.. record and anything else allocated in the non-common
//						.. part of the record.
//
// Returns:
//	(nothing)
//
//----------

//----------
//
// op_apply--
//	Apply operation to vector(s).
//
//----------
//
// Arguments:
//	dspop*		op:		Pointer to the operator's control record that was
//						.. allocated by op_parse.  This is compatible with
//						.. a dspop record, but usually this function will cast
//						.. it to the type of its private control record.
//	char*		vName:	The name of the vector.  Typically a chromosome name,
//						.. or "*" if the operator deals with multiple vectors.
//	u32			vLen:	Number of entries in v[].  If this is a "whole genome"
//						operator, vLen is the maximum number of entries in any
//						vector.
//	valtype*	v:		The vector to operate opon.
//
// Returns:
//	(nothing)
//
//----------

#endif // genodsp_interface_H

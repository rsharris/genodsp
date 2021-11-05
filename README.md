genodsp -- General workbench for processing signals along genomic intervals
==============================================
Jun/16/2015, Bob Harris (rsharris *at* bx *dot* psu *dot* edu)

GenoDSP (pronounced jen-OH-dee-us-PEE) is a workbench program for performing
DSP-like pipelines on genomic data.  DSP stands for digitial signal processing.

The basic premise is that we have a vector of values along the genome.  More
precisely a set of vectors, one along each chromosome.  Collectively these are
called "the signal".  We read the signal from a file, perform a series of
operations on it, and then output the resulting signal.  Input and output files
are typically four column chromosome-start-end-value files.

The program is a serious memory hog.  Values are double-precicison floating-
point.  On my machine these take 8 bytes, so for the human genome the program
will need about 24G bytes.

Each operation works by modifying the current signal, and has its own parameter
settings.  A pipeline is specified using the "=" character.  This is easier to
describe with an example.

    cat something.dat \
      | genodsp --chromosomes=my_genome.chroms \
          = smooth --window=101 \
          = localmax --neghborhood=11 \
      > peaks.dat
 
In the example, a file describing the signal is piped into genodsp, with the
lengths of the chromosomes given in my_genome.chroms.  We smooth the signal
over 101bp windows, then take local maxima in 11bp neghborhoods, and save the
result to a peaks file.

General usage syntax is shown if you run genodsp with no arguments.  This shows
arguments applicable to general operation (most relating to I/O), followed by
a one-line description of each operator.  To get detailed information on a
specific operator, do "genodsp --help=<operator>".


===== Build/Installation =====

The program is comprised of a few dozen .c and .h files.  Running the command
"make" in the directory with those files should produce the executable.  As of
this writing there are no unusual packages to download/build/install/swear at.
The only library linked with is the standard math library.

To install it, copy the executable somewhere into your shell's PATH.


===== Chromosome Lengths =====

The --chromosomes=<filename> option is mandatory.  This indicates how much
space should be allocated for each chromosome.

While it would be nice if the program automatically determined chromosome
lengths from the signal itself, that is problematic.  The first issue is that
some of the genome might be absent from the input signal, and thus we couldn't
know the true lengths.  Secondly, we'd have to make two passes over the input,
one to determine the lengths, then another pass after allocating the vectors.
(or we could dynamically grow the allocated vectors but that introduces other
issues.)

The chromosomes filename is a simple two column format, whitespace-separated.
The first column is the chromosome name.  The second is the chromosome's length.
The output of Jim Kent's twoBitInfo is compatible with this.

The order of chromosomes in this file is used to determine the order when we
output the signal.  So if you like your data sorted in a certain way, make sure
that's the way the chromosomes appear in this file.

The chromosomes in this file also determine the chromosomes that we're
interested in.  If the signal contains additional chromosomes, those will be
ignored.


===== Signal Input/Output =====

If your input contains overlapping intervals, the values at those positions are
summed.  Thus a very easy way to compute coverage depth is as follows:

    cat intervals.dat \
      | genodsp --chromosomes=my_genome.chroms --novalue \
      > depth.dat

In that command, --novalue tells the program to ignore any values in the input
signal file, and just assume every interval given has the value 1.  (So the
file is chromosome-start-end and any extra columns are ignored.)  As each
interval is read, 1 is added (to the vector) at every position in the interval.
So each position ends up with the count of intervals that covered that base.

By default the output is four columns.  Runs of the same signal value are
collapsed into one interval, one line of output.  Zeros are not output.  Values
are rounded to the nearest whole number.

Each of those default traits can be overridden by command line parameters.  See
the usage text for details.

Usually the program reads the signal, applies a series of operators, and then
outputs the result.  However, there are input and output operators, so it is
possible to write the result of any stage of the pipeline to a file.

Also, several operators read a second signal from a file and involve it in a
computation with the current signal.


===== Named Variables =====

In some cases an operator computes a value that another operator would like to
use.  An example is the percentile operator, with another operator using the
computed percentile as a threshold.  This is probably clearer with an example:

    cat signal.dat \
      | genodsp --chromosomes=my_genome.chroms \
          = sum --window=100 --denom=100 \
          = output signal.save \
          = percentile 99 --window=100 --min=1/inf --precision=3 \
          = input signal.save --destroy \
          = binarize --threshold=percentile99 \
      > high_signal.dat

In that example, the computed result of the percentile operation is a named
variable, "percentile99".  The binarize operation then uses that as its
threshold (binarize converts everything above the threshold to 1, everthing
below to 0).

That example also shows a use of input and output operators.  The percentile
operation is desctructive-- the state of the resulting signal is unusable.  So
we save the input signal to a file and then restore it afterwards (and delete
the temporary file) befor applying the threshold.  However, in many cases it
is faster to simply reread and recompute the signal instead of incurring the
cost of writing a file.

At present, only a few commands actually allow the use of named variables.


===== Finding Clumps =====

Here's an example that indentifies all intervals that average at the 99.5th
percentile or better over 100 bp.

    genodsp --chromosomes=my_genome.chroms \
      = input signal.dat \
      = sum --window=100 --denom=100 \
      = percentile 99.5 --window=100 --min=1/inf --quiet \
      = input signal.dat \
      = sum --window=100 --denom=100 \
      = clump --average=percentile99.5 --length=100 \
      > clumps.dat


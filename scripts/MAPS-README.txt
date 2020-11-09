Comai Lab, UC Davis Genome Center 
Meric Lieberman, Isabelle Henry, 2012 
This work is the property of UC Davis Genome Center - Comai Lab 

Use at your own risk. 
We cannot provide support. 
All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.


MAPS - Mutation and Polymorphism Survey 


This package uses the information present in a mapped coverage file, in this case a parsed samtools mpileup file, to either detect putative mutations or genotype information. This is performed in two steps. The first step uses basic cutoffs to select interesting positions with adequate coverage and library representation. The second step then allows fine tuning to carefully select a specific output profile, reducing the interesting position list to a smaller list of positions with the relevant mutation detection/genotyping data in an simple output format.

Please note that this can operate in either genotyping or mutation detection modes. The mode must be selected using command line parameters for both part one and part two. Mode cannot be switched between part one and part two, so be sure to use the same -m (--mode) parameter for both parts.

For multicore/threading purposes, the first program reads the entire file into memory before running, so it is recommended not to run it on systems with limited memory. It typically requires > 1.5 times the size of the parsed mpileup file in RAM to run. For example, if the mpileup.txt file is 20 gbps, > 30 Gbps of RAM are necessary to run the program on the whole file at once. Depending on the number of samples and the coverage, mpileup files can be very large and the required RAM will not be available. In this case, it is recommended to break the parsed mpileup file into smaller chunks to be processed separately, or leave them in chunks if they were already separated during the mpileup parsing step (i. e. by chromosome or scaffold, see ??? link to the bwa-doall page?).

Both programs use command line parameters to set options, a list of which can be found i) by using a simple -h flag after the program name while executing (./program.py -h), ii) here in the program descriptions, and iii) at the top of the programs themselves.

Note: To see run paramaters: ./maps-part1.py -h or ./maps-part1.py --help AND ./maps-part2.py -h or ./maps-part2.py --help

The program is run in two steps to allow the second step to be quickly run repeatedly, allowing fast optimization of the final cutoff parameters.


Run Example, for a parsed mpileup file parsed_test_mpileup.txt.

Step 1:

To run in mutation mode with 10 threads, minimum of 6 libraries, minimum coverage of 20, maximum coverage of 2000, and default heterozygous position parameters: ./maps-part1.py -f parsed_test_mpileup.txt -o test-maps-part1.txt -l 6 -t 10 -v 20 -c 2000 -m m OR ./maps-part1.py --file parsed_test_mpileup.txt --out test-maps-part1.txt --minlibs 6 --numthreads 10 --mincov 20 --maxcov 2000 --mode m

Step 2: Continue the run from above, with heterozygous mutation minimum coverage of 6, homozygous minimum coverage of 3, and heterozygous minimum mutant allele percentage of 20% ./maps-part2.py -f test-maps-part1.txt -o results-test-maps-part1.txt -c 20 -d 6 -p 20.0 -s 3 -l 6 -m m OR ./maps-part2.py --file test-maps-part1.txt --out results-test-maps-part1.txt --mincov 20 --hetmincov 6 --hetminper 20.0 --hommincov 3 --minlibs 6 --mode m

Below is the individual documentation for the two programs of the MAPS package.


MAPS - Part 1, maps-part1.py

This program takes a parsed mpileup file and generates a "interesting" position file that serves as input to the second part of the MAPS pipeline. This first step is typically run once while the second program can be run quickly many times to try different mutation cutoffs.

INPUT: This program takes a parsed mpileup.txt file as input.

OUTPUT: 
i. First, the user specifies an output file name for the positions of interest, this is used as input for the second MAPS program. 
ii. There is an assay-[input file name] file generated that is an assay of the coverage of all "good" positions that pass the basic cutoffs. 
iii. There is a snps-[input file name] file generated that records any "good" position that is a SNP position (a position that is polymporphic between the reference genome and the samples) across all libraries that have data for that position. 
iv. There is a types-[input file name] file generated that outputs the counts of the different types of positions found. Only "good" positions are counted

A position is defined as “good” if it passes the cutoffs for minimum position coverage, maximum position coverage, and minimum of libraries with data for that position. This is also what we mean by assayed.

What are the different types of positions: 
i. het- These are "good" positions that pass the het cutoff parameters specified. (explained below) 
ii. ref - These are "good" positions that match the reference. 
iii. snp - These are positions that across all lib data for a position are a polymorphic between the reference and the samples. 
iv. other - These are the "interesting" positions used in MAPS part 2. 
v. rejected - These are positions that do not pass the basic position cutoffs mentioned to be a "good" position above. 
vi. unknown - These are assayed positions for which the data that are too complex / strange to be useful in step 2.

NOTE: 
i. This first program reads the entire file into memory before running, so it is recommended not to run it on systems with limited memory. It typically requires > 1.5 times the size of the parsed mpileup file in RAM to run. For example, if the mpileup.txt file is 20 gbps, 30 Gbps of RAM are necessary to run the program on the whole file at once. Depending on the number of samples and the coverage, mpileup files can be very large and the required RAM will not be available. In this case, it is recommended to break the parsed mpileup file into smaller chunks to be processed separately, or leave them in chunks if they were already separated during the mpileup parsing step (i. e. by chromosome or scaffold. 
ii. This program is threaded, so it can be used with the -t flag to specify the number of cores to be used.

PARAMETERS, Default value in []: 
1. REQUIRED: 
i. -f or --file, The input parsed mpileup file. [required] 
ii. -o or --out, The output interesting position destination file. [required] 
2. OPTIONAL: 
i. -t or --thread, Number of cores to be used while processing [1] 
ii. -l or --MinLibs, Minimum number of libraries covered at least once to be considered a valid position. [3] 
iii. --minCov or -v, Minimum total coverage to be considered as a valid position. Total coverage is the sum of the coverage at that position in all invididual libraries. [10] 
iv. --maxCov or -c, Maximum total position coverage to be considered as a valid position. Total coverage is the sum of the coverage at that position in all invididual libraries. [2000] 
v. --hetBothMinPer or -b, Minimum percentage for consideration of heterozygous position for both bases combined [95.00], i.e. % SNP base 1 and %SNP base 2 percent must sum to above this cutoff. This is based on the total coverage at a position. For example, if a base has the following coverage: 20A (50%), 16T (40%) and 4G (10%). SNP1 is A and SNP2 is T. SNP1 + SNP2 = 90%. 
vi. --hetOneMinPer or -i, :Minimum percentage for consideration of heterozygous position for each of the two most common bases.[20.00] Each of the two SNP abse percentages must be at least this cutoff percent. 
vii. --minCovNonMutBase or -u, Minimum coverage of a non-mutant base [5] 
viii. --mode or -m, Output Mode: m == Mutation Detection, g == Genotyping [m]


MAPS - Part 2, maps-part2.py


This program uses the list of potentially interesting positions generated by part1 and filters each of them based on a second set of more or less stringent criteria. Typically, users will run this second step multiple times to optimize the criteria, in order to obtain a final output that best represents their input data.

INPUT: This program takes the list of potentially interesting position from the MAPS part 1 program as input.

OUTPUT: 
The user specifies an output file name for the output mutation / gentoyping file. 
There is an non-assay-[input file name] file generated that is a table of the assayed positions that are lost do to cutoff parameters.

PARAMETERS, Default value in []: 
1. REQUIRED: 
i. -f or --file, The input file is the output of the MAPS part 1 program [required] 
ii. -o or --out, The output mutation / gentotyping file [required] 
2. OPTIONAL: 
i. --minCov or -v, Minimum total position coverage to be considered as a valid position. Should be the same as MAPS part 1. [1]

The next four parameters analyze the basecalls at a particular position and for a particular library.

In genotyping mode: each library is dealt with independently. 
ii. --hetMinCov or -d, minimum coverage of each of the two observed major basecalls in order to be called heterozygous for that library. [5] 
iii. --hetMinPer or -p, minimum percentage of each of the two observed major basecalls in order to be called heterozygous for that library. [20.00] 
iv. --hetMinLibCov or -D, Minimum total coverage for consideration of heterozygous. This parameter is only used in genotyping mode. [5] 
v. --homMinCov or -s, Minimum coverage in order for a position to be called homozygous for that library. [2]

In mutation detection mode: a potential “mutant” allele is determined as present in only one library and not in any of the others (which all carry the WT allele). Whether or not it qualifies as a potential mutation depends on the following criteria: 
ii. --hetMinCov or -d, minimum coverage of mutant allele to be called heterozygous for that library [5] 
iii. --hetMinPer or -p, minimum percentage of the mutant and WT alleles in order to be called heterozygous for that library. [20.00] 
v. --homMinCov or -s, Minimum coverage in order for a position to be called homozygous for that library. [2]

In both modes again: 
vi. -l or --MinLibs, Minimum number of libraries covered at least once to be considered a valid position. Should be the same as MAPS part 1. [3] 
vii. --mode or -m, Output Mode: m == Mutation Detection, g == Genotyping. Should be the same as MAPS part 1. [m]


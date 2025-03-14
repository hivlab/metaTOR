
# Change Log

All notable changes to this project will be documented in this file.

## [1.1.4] - 2021-07-23

- Correct broken release 1.1.3

## [1.1.3] - 2021-07-23

- Avoid raising error if no pairs extracted in contact map module.
- Add network heatmap plot in the ouput.

## [1.1.2] - 2021-07-13

- Add the possibility to read pairix index files in the contact map module to retrieve faster pairs from a 2D region.

## [1.1.1] - 2021-06-30

Correct some minor issues:

- Correct pairs files:
  - Correct position index (previously 0-based but should be 1-based).
  - Change pairs entry to have an upper triangle shape.
- Correct recursive bin issue in the contig data file were unbinned contigs or binned contigs in bins smaller than the size threshold given were given an id of an existing bin.
- Correct issue in the contact map were contigs needed to be size sorted if a size threshold were given.
- Add some multiple network files in the output if multiple input samples were given.

## [1.1.0] - 2021-05-26

- Add some Figures in the output:
  - Completion/contamination distribution of the MAGs.
  - Distribution of the quality of the MAGs in the assembly.
  - Distribution of the GC and coverage inside each MAG.
- Add test for the new figure module.
- Add a module contactmap to build HiC matrix from some metaTOR objects such as contigs, bins, or other objects with contigs labels from the user. It could be use with metator contactmap [arguments]. The output could be use then with hicstuff view or cooler show to visualize matrix or with instagraal to try to scaffold MAGs.
- Modify the pairs files:
  - Change pairs file extension to .pairs.
  - Reorder the columns of .pairs file: readID-chr1-pos1-chr2-pos2-strand1-strand2 instead of readID-chr1-pos1-strand1-chr2-pos2-strand2.
  - Use tabulation to separate column instead of space to match the official specification.
- Correction of some minor issues:
  - Stop the duplication of the logging message.
  -Silence the HTSlib index error (from pysam.AlignementFile) which was meaningless.
  - Remove pyfastx index.
  - Put as optional the output of the clustering matrix which needs high memory usage for large assembly and optimize memory usage to build it.
  - Modify Assignment method to dataframe to avoid assignment to copy of the dataframe (SettingWithCopyWarning).

## [1.0.4] - 2021-05-14

- Change the format of the pairs file to match with the official specification.
- Add the HiC coverage (mean number of contacts in 1kb region) in the final bin summary.
- Add possibility to start with the pairs file.
- Add info about general alignment if multiples files are given.
- Organize temporary folder.

## [1.0.3] - 2021-05-04
  
- Add clustering matrix in the output of the partition, validation and pipeline modules.
- Add output files documentation.
- Add codeco.

## [1.0.2] - 2021-05-01

- Modification of the ouput files:
  - Add bin coverage based on the depth file given in the bin_summary file.
  - Add a final bin column in the contig data file after the validation to specify the name of the final bin of the contig. A "ND" value is given if the contig isn't binned.
  - Add the pairs file of the alignment to be able to produce HiC matrices from it using hicstuff for example.
  - Add a log file with all information of the run.
  - With the pipeline workflow remove intermediate contig data files.
- Modify call to input data to use only header names and not order of the column.
- Change the parameter of the bowtie2 alignement back to --very-sensitive-local.
- Add the total size of the different category of bins in the log.
- Correct minor issues in the doc.

## [1.0.1] - 2021-04-27

Correct bug in validation module:

- louvain_recursif: correct type error to identify contigs of a subnetwork which could triggers empty subnetwork error.

## [1.0.0] - 2021-04-26

New version of metaTOR in Python3

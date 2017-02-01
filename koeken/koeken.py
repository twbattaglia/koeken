#!/usr/bin/env python
"""
    koeken.py
    ~~~~~~~~~
    A Linear Discriminant Analysis effect size (LEfSe) wrapper.
    :copyright: (c) 2016 by Thomas W. Battaglia.
    :license: BSD, see LICENSE for more details.
"""
__author__ = 'Thomas W. Battaglia'
__copyright__ = 'Copyright 2016'
__license__ = 'BSD'
__version__ = '0.3.0'
__email__ = 'tb1280@nyu.edu'
__status__ = 'Development'
description = """
Koeken
=================================
This package is meant to provide better integration of the LEfSe into a
typical metagneomic analysis workflow. The statistical analysis can be
applied to QIIME, PICRUSt or humann2 datasets, without the need to manually
add metadata. It was developed to reduce the risk of incorrect metadata
and increased the reproducibility of analytical studies.
"""

try:
	import pandas as pd
except ImportError:
	raise ImportError('The module "pandas" could not be found.')
try:
	import qiime
except ImportError:
	raise ImportError('The module "qiime" could not be found.')
try:
	import biom
except ImportError:
	raise ImportError('The module "biom-format" could not be found.')
try:
	import rpy2
except ImportError:
	raise ImportError('The module "rpy2" could not be found.')
import util
import sys
import re
import argparse


# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(
        description = description,
		add_help = True,
        formatter_class = argparse.RawTextHelpFormatter,
		prog = "koeken")
    parser.add_argument(
        "-v", "--version",
        action = "version",
		version = "%(prog)s v"+ __version__)
    parser.add_argument(
        "-i", "--input",
		dest = "input_fp",
        metavar = "<file>",
		type = str,
        required = True,
        help = "File path to input feature table.\n")
    parser.add_argument(
        "-o", "--output",
		dest = "output_dir",
        metavar = "<folder>",
        type = str,
		required = True,
        help = "Folder path to store output files.\n")
    parser.add_argument(
        "-m", "--mapping",
        metavar = "<file>",
		type = str,
        required = True,
        help = "Path to sample metadata.")
    parser.add_argument(
        "-f", "--format",
        metavar="<string>",
		type = str,
		choices = ["qiime", "picrust", "humann2"],
        required = True,
        help = "Set format type to run analysis. Can be of type (qiime, picrust"
			   " or humann2)")

	# Metadata Params
    parser.add_argument(
        "--class",
        metavar = "<string>",
		dest = "classid",
        type = str,
		required = True,
        help = "Select which variable in metadata to use as class\n")
    parser.add_argument(
        "--subclass",
        metavar = "<string>",
        type = str,
        help = "Select which variable in metadata to use as subclass.")
    parser.add_argument(
        "--subject",
        metavar = "<string>",
        type = str,
		default = '#SampleID',
        help = "Select which variable in metadata to use as subject. "
			   "[DEFAULT: first column in mapping file]")

	# Splitting Params
    parser.add_argument(
		"--compare",
        metavar = "<string>",
        type = str,
		default = "",
		nargs = '+',
        help = "Select the factors within the class variable that you would \n"
			   "like to compare. Be sure to use a space to separate the \n"
			   "different factors. (e.g Treatment1 Treatment2 Control)")
    parser.add_argument(
        "--split",
        metavar = "<string>",
        type = str,
        help = "Select which variable in metadata to split the table by. \n"
			   "Typically this variable represents the timepoint in the data.")
    parser.add_argument(
        "--no-split",
		dest = "no_split",
		action = "store_true",
        help = "Select to only compare all timepoint together [DEFAULT: False]")

	# QIIME Params
    parser.add_argument(
        "-l", "--level",
        metavar = "<number>",
        type = int ,
		choices = [2, 3, 4, 5, 6, 7],
		default = 6,
        help = "Level at which to collapse taxa (Only for QIIME) [DEFAULT : 6]\n"
			   "Level 1 = Kingdom (e.g Bacteria)\n"
			   "Level 2 = Phylum (e.g Actinobacteria)\n"
			   "Level 3 = Class (e.g Actinobacteria)\n"
			   "Level 4 = Order (e.g Actinomycetales)\n"
			   "Level 5 = Family (e.g Streptomycetaceae)\n"
			   "Level 6 = Genus (e.g Streptomyces)\n"
	    	   "Level 7 = Species (e.g mirabilis)\n")
    parser.add_argument(
        "-n", "--unclassified",
        metavar = "<string>",
        type = str ,
		default = "unclassified",
        help = "The name used to replace taxa without a complete phylogeny \n"
			   "[DEFAULT: unclassified]")

	# LEfSe Options
    parser.add_argument(
        "--pvalue",
        metavar = "<number>",
        type = float,
		default = 0.05,
        help = "Maximum p-value for LEfSe analysis. [DEFAULT: 0.05]" )
    parser.add_argument(
        "--lda",
		metavar = "<number>",
        type = float,
		default = 2.0,
        help = "Minimum LDA score for LEfSe analysis. [DEFAULT: 2]")
    parser.add_argument(
        "--strictness",
		metavar = "<number>",
        type = int,
		choices = [0,1],
		default = 0,
        help = "Set the type of analysis for 2+ group comparisons. \n"
			   "0 = All-against-all (more strict) [DEFAULT]. \n"
			   "1 = One-against-all (less strict)")
    parser.add_argument(
        "--image-type",
		metavar = "<string>",
		dest = "image_type",
        type = str,
		choices = ["png", "pdf", "svg"],
		default = "pdf",
        help = "Format-type of output cladogram image. [DEFAULT: pdf] ")
    parser.add_argument(
        "--dpi",
		metavar = "<number>",
		type = int,
    	default = 300,
        help = "Resolution of output cladogram image. [DEFAULT: 300] ")
    args = parser.parse_args()
    return args


# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------
def main():
	print("\n=================================================================")
	print('\033[95m' + "koeken" + ' v' + __version__ + '\033[0m')
	print('Written by: ' + '\033[94m' + __author__ + ' (' + __email__ + ') ' +
        '\033[0m' + '\n')

	print('\033[92m' + 'LEfSe Publication:\n'
         'Metagenomic biomarker discovery and explanation.\n'
         'Nicola Segata, Jacques Izard, Levi Waldron, Dirk Gevers\n'
	     'Larisa Miropolsky, Wendy S Garrett, and Curtis Huttenhower')
	print('Genome Biology, 12:R60, 2011' + '\033[0m')
	print("=================================================================\n")

	# Get arguments
	args = parse_arguments()
	print(args)

	# Check data for errors
	util.check_map(args)

	# Create folders to store output
	util.create_dir(args.output_dir)
	util.create_dir(args.output_dir + "/formatted")
	util.create_dir(args.output_dir + "/results")
	util.create_dir(args.output_dir + "/split_tables")
	util.create_dir(args.output_dir + "/cladograms")

	# Run summarize taxa on BIOM file
	if args.format == "qiime" or args.format == "picrust":

		# Run summarize taxa command on BIOM file
		util.summarize_taxa(args)

		# Import data as a pandas dataframe
		sumtbl_df = pd.read_table(args.output_dir + "/summarize_table.txt")
		map_df = pd.read_table(args.mapping)

		# Set subject ID column as first row as default
		if args.subject != None:
			subjectID_pos = sumtbl_df.columns.get_loc(args.subject)
		else:
			subjectID_pos = 1

		# Find the cols and respective positions for input variables in the table.
		### Row 1 = Subject/First row in mapping file
		### Row 2 = Class
		### Row 3 = Subclass/None
		### Row 4-: = Bacteria/Features...
		classID_pos = sumtbl_df.columns.get_loc(args.classid)
		bacteria_pos = range(( len(map_df.columns) ), len(sumtbl_df.columns))

		# Cases for addition of subject class
		if args.subclass is None:
			to_keep = [subjectID_pos, classID_pos] + bacteria_pos
		else:
			subclassID_pos = sumtbl_df.columns.get_loc(args.subclass)
			to_keep = [subjectID_pos, classID_pos, subclassID_pos] + bacteria_pos

		# Subset the data if particular group comparisons are given
		if args.compare != "":
			sumtbl_df=sumtbl_df[sumtbl_df[args.classid].isin(args.compare)]

		# Remove greengenes taxa names to make it prettier
		# TODO fix the double || error
		if args.format=="qiime":

            # Replace blank names (e.g g__)
			sumtbl_df=sumtbl_df.rename(columns=lambda x: re.sub('__$', "__" + args.unclassified, x))

            # Remove greengenes ID (e.g k__|p__|c__ )
			sumtbl_df=sumtbl_df.rename(columns=lambda x: re.sub('.__', '', x))

            # Remove names
			#sumtbl_df=sumtbl_df.rename(columns=lambda x: re.sub('', args.unclassified, x))

		# Run the analysis with all timepoints merged
		all_split=args.output_dir + "/split_tables/all_levels.txt"
		all_format=args.output_dir + "/formatted/all_levels.txt"
		all_results=args.output_dir + "/results/all_levels.txt"
		all_clade=args.output_dir + "/cladograms/all_levels."+ args.image_type

        # Transpose table and write to disk
		sumtbl_df_tsp = sumtbl_df.iloc[:,to_keep].transpose()
		sumtbl_df_tsp.to_csv(all_split, sep='\t', header=False, index=True)

		# Format data step
		# TODO add conditonal for subclass
		util.format_lefse(input_fp = all_split,
						  output_fp = all_format,
						  name = "all samples",
						  subclass = None)

		# Running lefse step
		util.run_lefse(input_fp = all_format,
					   output_fp = all_results,
					   name = "all samples",
					   args = args)

		# Plot cladogram step
		util.plot_cladogram(input_fp = all_results,
							output_fp = all_clade,
							name = "all samples",
							args = args)

		# Aesthetically pleasing
		print('=========================================================\n')

		if not args.no_split:
			# Split table to iterate over multiple timepoints
			grouped_df = sumtbl_df.groupby(str(args.split))

			# Iterate over each splitted mapping file
			for name, group in grouped_df:
				# Set vars
				current_split=args.output_dir + "/split_tables/" + str(name) + "_split.txt"
				current_format=args.output_dir + "/formatted/" + str(name) + "_format.txt"
				current_results=args.output_dir + "/results/" + str(name) + ".txt"
				current_clado=args.output_dir + "/cladograms/" + str(name) +'.' + args.image_type

				# Subset table to current working timepoint
				table = group.iloc[:,to_keep].transpose()

				# Remove any rows with 0 sums
				table_filtered = table.loc[~(table==0).all(axis=1)]

				# Write Input tables to file
				table_filtered.to_csv(current_split,
									  sep = '\t',
									  header = False,
									  index = True)

				# Formatting step
				# TODO add conditonal for subclass
				util.format_lefse(input_fp = current_split,
								  output_fp = current_format,
								  name = name,
								  subclass = None)

				# Run lefse step
				util.run_lefse(input_fp = current_format,
							   output_fp = current_results,
							   name = name,
							   args = args)

				# Plot cladogram step
				util.plot_cladogram(input_fp = current_results,
									output_fp = current_clado,
									name = name,
									args = args)

				# Aesthetically pleasing
				print('=========================================================\n')
	else:
		print("Running humann2 workflow")

if __name__ == '__main__':
    main()

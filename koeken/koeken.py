#!/usr/bin/env python
"""
	pannenkoek.py
	~~~~~~~~~
	A simple command line application to run LEfSe multiple times
	on the same OTU table.
	:copyright: (c) 2015 by Thomas W. Battaglia.
	:license: BSD, see LICENSE for more details.
"""

__author__ = 'Thomas W. Battaglia'
__copyright__ = 'Copyright 2015'
__license__ = 'BSD'
__version__ = '0.2.1'
__email__ = 'tb1280@nyu.edu'
__status__ = 'Development'

import os
import os.path
import shlex
import sys
import glob
import argparse
import subprocess
import re
#import pretty_lefse as pl

try:
	import pandas as pd
except ImportError:
	raise ImportError('The module """pandas"" was not found. '
	'Please install with """pip install pandas""" and try again')

try:
	import qiime
except ImportError:
	raise ImportError('The module """QIIME"" was not found. '
	'Use either MacQIIME or """pip install qiime""" before running the command.')







def get_args():
	"""Gets arguments from command line inputs"""
	parser = argparse.ArgumentParser(description = 'Performs Linear Discriminant Analysis (LEfSe) on A Longitudinal Dataset.', add_help = True)
	parser.add_argument('-v', '--version', action = 'version', version = __version__)


	"""Arguments for summarize_taxa.py inputs"""
	parser.add_argument('-i', '--input', action = "store", dest = "input", help = 'Location of the OTU Table for main analysis. (Must be .biom format)', required = True, type=str)
	parser.add_argument('-o', '--output', action = "store", dest = "output", help = 'Location of the folder to place all resulting files. If folder does not exist, the program will create it.', required = True, type=str)
	parser.add_argument('-m', '--map', action = "store", dest = "map", help = 'Location of the mapping file associated with OTU Table.', required = True, type=str)
	parser.add_argument('-l', '--level', action = "store", dest = "level", default = 6, help = 'Level for which to use for summarize_taxa.py. [default = 6]', type=int, choices=[2,3,4,5,6,7])


	"""Arguments for LEfSe inputs"""
	parser.add_argument('-cl', '--class', action = "store", dest = "classid", help = 'Location of the OTU Table for main analysis. (Must be .biom format)', required = True, type = str)
	parser.add_argument('-sc', '--subclass', action = "store", dest = "subclassid", default = 'NA', help = 'Directory to place all the files.', type = str)
	parser.add_argument('-su', '--subject', action = "store", dest = "subjectid", default = '#SampleID', help = 'Only change if your Sample-ID column names differs. [default] = #SampleID.', type = str)


	"""Arguments for LEfSe discriminant analysis"""
	parser.add_argument('-p', '--pval', action = "store", dest = "p_cutoff", default = 0.05, help = 'Change alpha value for the Anova test (default 0.05)')
	parser.add_argument('-e','--effect', action = "store", dest = "lda_cutoff", default = 2, help = 'Change the cutoff for logarithmic                  LDA score (default 2.0).', type = float)
	parser.add_argument('-str', '--strict', action = "store", dest = "strictness", default = 0, choices=[0,1], help = 'Change the strictness of the comparisons. Can be changed to less strict(1). [default = 0](more-strict).', type = int)


	"""Arguments for organizing data"""
	parser.add_argument('-c', '--compare', action = "store", dest = "compare", default = "NA", nargs = '+', help = 'Which groups should be kept to be compared against one another. [default = all factors]', required = False, type = str)
	parser.add_argument('-sp', '--split', action = "store", dest = "split", default = "NA", help = 'The name of the timepoint variable in you mapping file. This variable will be used to split the table for each value in this variable.', required = True, type = str)


	"""Arguments for other types of analyses"""
	parser.add_argument('-pc', '--clade', action = "store_true", dest = "clade", help = 'Plot Lefse Cladogram for each output time point. Outputs are placed in a new folder created in the lefse results location. Default file type in (.png)', default = False)
	parser.add_argument('-g', '--graphlan', action = "store_true", dest = "graphlan", help = 'Convert LEfSe Results to Graphlan compatible tree (Experimental)', default = False)

    #parser.add_argument('-py', '--pretty', action = "store_true", dest = "pretty", help = 'Enable this option if you would like to generate a pretty lefse table from the results. This table will consist of a table of all LEfSe results across each longitudinal timepoint and the corresponding effect size. This command can currently only work with 2 groups. See the github page for more information.', default = True)

    #parser.add_argument('-con', '--control', action = "store", dest = "plcontrol", help = 'Enable this option if you would like to generate a pretty lefse table from the results. This table will consist of a table of all LEfSe results across each longitudinal timepoint and the corresponding effect size. This command can currently only work with 2 groups. See the github page for more information.')

	return parser.parse_args()




def main(input, output, map, level, classid, subclassid, subjectid, compare, split, p_cutoff, lda_cutoff, strictness, clade, graphlan):
	"""Run QIIME's summarize_taxa.py command to generate relative abundance tables with associated metadata.
	Then split the table by input parameter and """

	"""Check to see if output directories exist or not and create them."""
	if not os.path.exists(output):
		os.makedirs(output)
	else:
		print "Output folder already exists. Warning: Errors may be produced."
		print "Please delete or change output folder before running again!. "+ '\n'

	"""Output location from summarize_taxa.py step"""
	sumtaxa_dir = '{}/{}{}/'.format(output, "summarize_taxa_L", str(level))
	if not os.path.exists(sumtaxa_dir):
		os.makedirs(sumtaxa_dir)

	"""Output location for all LEfSe analyses"""
	lefse_dir = '{}/{}'.format(output, "lefse_output")
	if not os.path.exists(lefse_dir):
		os.makedirs(lefse_dir)

	"""Output location for LEfSe formatting step"""
	format_dir = '{}/{}/'.format(lefse_dir, "format_lefse")
	if not os.path.exists(format_dir):
		os.makedirs(format_dir)

	"""Output location for LEfSe analysis step"""
	run_dir = '{}/{}/'.format(lefse_dir, "run_lefse")
	if not os.path.exists(run_dir):
		os.makedirs(run_dir)

	"""Output location for LEfSe Cladogram step"""
	if clade == True:
		clado_dir = '{}/{}/'.format(run_dir, "cladograms")
		if not os.path.exists(clado_dir):
			os.makedirs(clado_dir)

	"""Output location for Pretty-LEfSe step"""
    #pretty_dir = '{}/{}'.format(lefse_dir, "pretty_lefse")
    #if not os.path.exists(pretty_dir):
    #os.makedirs(pretty_dir)


	"""Run summarize_taxa.py command"""
	print "Running QIIME's summarize_taxa.py... "+ '\n'
	summarize_cmd = "summarize_taxa.py -i %s -o %s -m %s -L %d -d '|'" % (input, sumtaxa_dir, map, level)
	subprocess.call(shlex.split(summarize_cmd))

	"""Get filename of generated summarize_taxa.py outputs"""
	sumtaxa_loc = glob.glob(sumtaxa_dir + '*_L'+ str(level) +'.txt')

	"""Create panda dataframes from summarize_taxa output file and mapping file"""
	sumtaxa_df = pd.read_table(sumtaxa_loc[0])
	map_df = pd.read_table(map)


	"""Find the cols and respective positions for input variables on the table."""
	""" Row 1 = Subject"""
	""" Row 2 = Class"""
	""" Row 3 = Subclass/None"""
	""" Row 4-: = Bacteria"""
	subjectID_pos = sumtaxa_df.columns.get_loc(subjectid)
	classID_pos = sumtaxa_df.columns.get_loc(classid)
	bacteria_pos = range((len(map_df.columns)),len(sumtaxa_df.columns)) # Find the number of cols in mapping file

	"""Cases for addition of subject class"""
	if subclassid == "NA":
		to_keep = [subjectID_pos,classID_pos] + bacteria_pos
	else:
		subclassID_pos = sumtaxa_df.columns.get_loc(subclassid)
		to_keep = [subjectID_pos, classID_pos, subclassID_pos] + bacteria_pos


	"""Subset the data if particular group comparisons are given"""
	if compare != "NA":
		sumtaxa_df = sumtaxa_df[sumtaxa_df[str(classid)].isin(compare)]


	""" Remove greengenes taxa names to makeit prettier """
	sumtaxa_df = sumtaxa_df.rename(columns = lambda x: re.sub('.__', '', x))

	"""For each timepoint, remove unwanted columns and perform LEfSe analysis"""
	grouped_df = sumtaxa_df.groupby(str(split))
	for name, group in grouped_df:

		"""Write Input tables to file"""
		table = group.iloc[:,to_keep].transpose()
		table_out = '{}{}{}'.format(sumtaxa_dir, name, '_input.txt')
		table.to_csv(table_out, sep = '\t', header = False, index = True)

		"""Run format_input.py from LEfSe package"""
		format_file_out = format_dir + os.path.basename(table_out).replace('_input.txt', '_format.txt')
		print 'Timepoint: ' + str(name)
		print 'Formatting Table...'
		print 'Formatting Input: ' + table_out
		print 'Formatting Output: ' + format_file_out
		if subclassid == "NA":
			subprocess.call(['format_input.py', table_out, format_file_out, '-u 1', '-c 2', '-o 1000000', '-f', 'r'])
		else:
			subprocess.call(['format_input.py', table_out, format_file_out, '-u 1', '-c 2', '-s 3', '-o 1000000', '-f', 'r'])


		"""Run run_lefse.py from LEfSe package"""
		run_file_out = run_dir + os.path.basename(format_file_out).replace('_format.txt', '.txt')
		print 'Running Analysis...'
		print 'Analysis Input: ' + format_file_out
		print 'Analysis Output: ' + run_file_out
		subprocess.call(['run_lefse.py', format_file_out, run_file_out, '-a', str(p_cutoff), '-l', str(lda_cutoff), '-y', str(strictness)])


		"""Check to see if cladogram option was chosen"""
		if clade == True:
			"""Run plot_cladogram.py from LEfSe package"""
			clade_file_out = clado_dir + os.path.basename(format_file_out).replace('_format.txt', '.png')
			print 'Plotting Cladogram...'
			#print 'Plot Input: ' + run_file_out
			#print 'Plot Output: ' + clade_file_out
			subprocess.call(['plot_cladogram.py', run_file_out, clade_file_out, '--format', 'png', '--dpi', '300'])

		print('\n')


		"""Check to see if graphlan option was chosen"""
		#if graphlan == True:
			#"""Convert files to graphlan compatibile files"""
			#clade_file_out = run_dir + os.path.basename(format_file_out).replace('_format.txt', '.png')
			#print 'Convert to Graphlan Format...'
			#print 'Graphlan Input: ' + run_file_out
			#print 'Graphlan Output Folder: ' + clade_file_out

			#subprocess.call(['export2graphlan.py', run_file_out, table_out, '--annotations', '2,3', 'png', '--external_annotations', '4,5,6', --fname_row 0])


		#export2graphlan.py -i hmp_aerobiosis_small.txt -o hmp_aerobiosis_small.res -t tree.txt -a annot.txt --title "HMP aerobiosis" --annotations 2,3 --external_annotations 4,5,6 --fname_row 0 --skip_rows 1,2 --ftop 200


		#"""Run Analysis on all combined timepoints"""
		#print 'Timepoint: All'
		#all_table = sumtaxa_df.iloc[:,to_keep]
		#print 'Formatting Table...'
		#print 'Formatting Input: ' + table_out
		#print 'Formatting Output: ' + format_file_out


    #print "Running Pretty LEfSe..."
    #frames = [ pl.process_table(file_loc = f, write_values = False, output = pretty_dir, control = pl_control) #for f in glob.glob(run_dir + "/*.txt") ]
    #result = pd.concat(frames, axis = 1)
    #result.to_csv(pretty_dir + "/pretty_table.txt", sep = '\t', header = True, index = True)
	print 'Analysis Completed.'




if __name__ == '__main__':

	"""Get input arguments and run command"""
	args = get_args()

	"""Credits"""
	print "Koeken" + ' v' + __version__ + ': ' + "Linear Discriminant Analysis (LEfSe) on A Longitudinal Microbial Dataset."
	print 'Written by ' + __author__ + ' (' + __email__ + ')' + '\n'
	print 'LEfSe Credits: "Metagenomic biomarker discovery and explanation"'
	print 'Nicola Segata, Jacques Izard, Levi Waldron, Dirk Gevers, Larisa Miropolsky, Wendy S Garrett, and Curtis Huttenhower'
	print 'Genome Biology, 12:R60, 2011'+ '\n'


	"""Error Check to see if class/subclass/split metadata columns exist in the mapping file"""
	map_chk = pd.read_table(args.map)
	if str(args.classid) not in map_chk.columns.values.tolist():
		raise ValueError('Warning. There is no class variable with that column name in your mapping file. Please verify the class ID chosen is actually a column name in your mapping file.')
	if str(args.split) not in map_chk.columns.values.tolist():
		raise ValueError('Warning. There is no split variable with that column name in your mapping file. Please verify the subclass ID chosen is actually a column name in your mapping file.')


	"""Run command"""
	main(args.input, args.output, args.map, args.level, args.classid, args.subclassid, args.subjectid, args.compare, args.split, args.p_cutoff, args.lda_cutoff, args.strictness, args.clade, args.graphlan)
    #args.pretty, args.plcontrol)

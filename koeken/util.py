#!/usr/bin/env python
'''
	util.py
	~~~~~~~~~
	Utility functions for koeken.
	:copyright: (c) 2016 by Thomas W. Battaglia.
	:license: BSD, see LICENSE for more details.
'''
import os
import subprocess
import pandas as pd
from qiime.summarize_taxa import make_summary, add_summary_mapping
from qiime.parse import parse_mapping_file
from qiime.format import (
    write_add_taxa_summary_mapping, format_add_taxa_summary_mapping)
from biom.table import Table
from biom import load_table


def create_dir(path):
	if not os.path.exists(path):
		os.makedirs(path)
		print('Created directory: ' + path)


def summarize_taxa(args):
	summarize_fp = args.output_dir + "/summarize_table.txt"
	otu_table = load_table(args.input_fp)
	mapping_file = open(args.mapping, 'U')
	mapping, header, comments = parse_mapping_file(mapping_file)

	if args.format == "picrust":

		# Collapse table by phylogeny level
		print("Summarizing PICRUSt functions...")
		summary, tax_order = add_summary_mapping(otu_table, mapping, 3,
		False, "KEGG_Pathways")

		# Write file to disk
		write_add_taxa_summary_mapping(summary, tax_order, mapping,
		header, summarize_fp , "|")

	elif args.format == "qiime":

		# Collapse table by phylogeny level
		otu_table = otu_table.norm(axis = 'sample', inplace = False)

		# Summarize taxa grouping
		print("Summarizing OTU table...")
		summary, tax_order = add_summary_mapping(otu_table, mapping,
		args.level, False, "taxonomy")

		# Write file to disk
		write_add_taxa_summary_mapping(summary, tax_order, mapping,
									   header, summarize_fp , "|")
	else:
		ValueError('Warning. humann2 is not supported in this function.')


def check_map(args):
	map_chk = pd.read_table(args.mapping)

	# Err if no class ID in mapping file columns
	if str(args.classid) not in map_chk.columns.values.tolist():
		raise ValueError('Warning. There is no class variable with that column '
		'name in your mapping file. Please verify the class ID chosen exists '
		'as a column name in your mapping file.')

	# Err if no variable to split by
	if not args.no_split and str(args.split) not in map_chk.columns.values.tolist():
		raise ValueError('Warning. There is no split variable with that '
		'column name in your mapping file. Please verify the splitting ID '
		'chosen exists as a column name in your mapping file.')

	# Err if comparison param did not use spaces.
	if(args.compare) != '':
		if (len(args.compare) == 1):
			raise ValueError('Warning. Comparison needs to be in the format '
			'of Group1 Group2. Only use spaces for separation.')


def format_lefse(input_fp, output_fp, name, subclass=None):
	print('Formatting data for: ' + str(name))
	if subclass is not None:
		subprocess.call(['lefse-format_input.py', input_fp, output_fp,
					     '-u 1', '-c 2', '-s 3', '-o 1000000', '-f', 'r'])
	else:
		subprocess.call(['lefse-format_input.py', input_fp, output_fp,
						 '-u 1', '-c 2', '-o 1000000', '-f', 'r'])

def run_lefse(input_fp, output_fp, name, args):
	print('Running LEfSe on: ' + str(name))
	subprocess.call(['run_lefse.py', input_fp, output_fp,
					 '-a', str(args.pvalue),
					 '-l', str(args.lda),
					 '-y', str(args.strictness)])

def plot_cladogram(input_fp, output_fp, name, args):
	print('Plotting cladogram for: ' + str(name))
	subprocess.call(['lefse-plot_cladogram.py', input_fp, output_fp,
					 '--format', args.image_type,
					 '--dpi', str(args.dpi),
					 '--title', str(name)])

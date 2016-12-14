#!/usr/bin/env python
'''
	util.py
	~~~~~~~~~
	Utility functions for koeken.
	:copyright: (c) 2016 by Thomas W. Battaglia.
	:license: BSD, see LICENSE for more details.
'''
import os
import pandas as pd
from qiime.summarize_taxa import make_summary, add_summary_mapping
from qiime.parse import parse_mapping_file
from qiime.format import (
    write_add_taxa_summary_mapping, format_add_taxa_summary_mapping)
from biom.table import Table
from biom import load_table


def summarize_taxa(args):
	summarize_fp = args.output_dir + "/summarize_taxa.txt"
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
		print("Summarizing OTU table")
		summary, tax_order = add_summary_mapping(otu_table, mapping,
		args.level, False, "taxonomy")

		# Write file to disk
		write_add_taxa_summary_mapping(summary, tax_order, mapping,
									   header, summarize_fp , "|")
	else:
		ValueError('Warning. humann2 is not supported in this function.')


def check_format(input):
    print('this is of type')

def create_dir(path):
	if not os.path.exists(path):
		os.makedirs(path)
		print('Created directory: ' + path)

def check_map(args):
	map_chk = pd.read_table(args.mapping)

	# Err if no class ID in mapping file columns
	if str(args.classid) not in map_chk.columns.values.tolist():
		raise ValueError('Warning. There is no class variable with that column '
		'name in your mapping file. Please verify the class ID chosen exists '
		'as a column name in your mapping file.')

	# Err if no variable to split by
	if str(args.split) not in map_chk.columns.values.tolist():
		raise ValueError('Warning. There is no split variable with that '
		'column name in your mapping file. Please verify the splitting ID '
		'chosen exists as a column name in your mapping file.')

	# Err if comparison param did not use spaces.
	if(args.compare) != '':
		if (len(args.compare) == 1):
			raise ValueError('Warning. Comparison needs to be in the format '
			'of Group1 Group2. Only use spaces for separation.')


def run_lefse(input, output, lda, pvalue):
	print('running lefse')


def format_lefse(input, output, subject = False):
	print('formatting lefse')


def plot_cladogram(input, output, format, dpi):
	print('plotting cladograms')

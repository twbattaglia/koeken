#!/usr/bin/env python
"""
	koeken.py
	~~~~~~~~~
	A Linear Discriminant Analysis effect size (LEfSe) wrapper.
	:copyright: (c) 2016 by Thomas W. Battaglia.
	:license: BSD, see LICENSE for more details.

"""
__author__ = 'Thomas W. Battaglia'
__copyright__ = 'Copyright 2015'
__license__ = 'BSD'
__version__ = '0.3.0'
__email__ = 'tb1280@nyu.edu'
__status__ = 'Development'
description = """
Koeken (command-line LEfSe)
===========================
Run LEfSe on multiple times.
"""

import sys
import argparse
import pandas
#import qiime



# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=description,
		add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
		prog="koeken")
    parser.add_argument(
        "-v", "--version",
        action="version",
		version="%(prog)s v"+ __version__)
    parser.add_argument(
        "-i", "--input",
        metavar="<path>",
		type=str,
        #required=True,
        help="Path to OTU/PICRUSt or humann2 feature table.")
    parser.add_argument(
        "-o", "--output",
        metavar="<path>",
        type=str,
		#required=True,
        help="Folder path to output files.")
    parser.add_argument(
        "-m", "--mapping",
        metavar="<path>",
		type=str,
        #required=True,
        help="Path to sample metadata.")
    parser.add_argument(
        "-f", "--format",
        metavar="<str>",
		type=str,
		choices=["qiime", "picrust", "humann2"],
        #required=True,
        help="Path to sample metadata. Must specifiy the type of input. Can be of type (qiime, picrust, humann2) ")

	# Metadata Params
    parser.add_argument(
        "--class",
        metavar="<string>",
        type=str,
        help="Class-definition")
    parser.add_argument(
        "--subclass",
        metavar="<string>",
        type=str,
        help="Subclass-definition")
    parser.add_argument(
        "--subject",
        metavar="<string>",
        type=str,
		default = '#SampleID',
        help="Subject-definition")

	# Splitting Params
    parser.add_argument(
		"--compare",
        metavar="<string>",
        type=str,
		nargs = '+',
        help="Compare-definition")
    parser.add_argument(
        "--split",
        metavar="<string>",
        type=str,
        help="Split-definition")
    parser.add_argument(
        "--no-split",
        action="store_true",
        help="No split-definition")

	# QIIME Params
    parser.add_argument(
        "-l", "--level",
        metavar="<int>",
        type=int,
		choices=[2, 3, 4, 5, 6, 7],
        help="Level-definition")

	# LEfSe Options
    parser.add_argument(
        "--pvalue",
        metavar="<float>",
        type=float,
		default=0.5,
        help="Pvalue-definition")
    parser.add_argument(
        "--lda",
		metavar="<float>",
        type=float,
		default=2.0,
        help="LDA-definition")
    parser.add_argument(
        "--strictness",
		metavar="<int>",
        type=int,
		choices=[0, 1],
		default=0,
        help="Strictness-definition")
    parser.add_argument(
        "--clade",
        action="store_true",
        help="Clade-definition")
    parser.add_argument(
        "--image-type",
		metavar="<str>",
        type=str,
		choices=["png", "pdf", "svg"],
		default="pdf",
        help="Image-type definition")
    parser.add_argument(
        "--dpi",
		metavar="<int>",
        type=int,
		default=300,
        help="DPI definition")
    args = parser.parse_args()
    return args


# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------
def main():
    args=parse_arguments()
	


if __name__ == "__main__":
    main()

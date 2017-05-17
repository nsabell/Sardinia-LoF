#!/usr/bin/env python
# Nathan Abell

# Last Usage: parse_VEP_vcfs.py <input_vcf> <output_name>

import sys
import HTSeq
import gzip
import argparse
import re

parser = argparse.ArgumentParser(description='Recover')
parser.add_argument('-i', dest='input_vcf', required=True, help='Path to gzipped, VEP-annotated vcf')
parser.add_argument('-o', dest='output_name', required=True, help='Path to gzipped output text file')
args = parser.parse_args()

with gzip.open(args.input_vcf, 'r') as vcf, gzip.open(args.output_name, 'w') as output:

	header = 'variant\tmaf\tclass'
	output.write(header + '\n')

	for line in vcf:

		if line.startswith("#"):
			if line.startswith("##INFO=<ID=CSQ"):
				vep_field_names = line.split('Format: ')[-1].strip('">\n').split('|')

		else:
			# Parsing code from loftee github
			fields = line.strip().split('\t')
			varID = "_".join(fields[0:2]) + "_" + "_".join(fields[3:5])
			info_field = dict([(x.split('=', 1)) for x in re.split(';(?=\w)', fields[7]) if x.find('=') > -1])

			if 'CSQ' in info_field:
				annotations = [dict(zip(vep_field_names, x.split('|'))) for x in info_field['CSQ'].split(',')]
			lof_annotations = [x for x in annotations if 'LoF' in x and x['LoF'] == 'HC']

			if lof_annotations:
				output.write(varID + "\t" + info_field['AF'] + "\tHC_LoF" + "\n")









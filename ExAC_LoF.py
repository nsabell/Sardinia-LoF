#!/usr/bin/env python
# Nathan Abell

# Last Usage: nohup ExAC_AF_byPop.py > ExAC_AF_withLoF_byPop.log 2>&1 &
# Retrieves allele frequencies for all variants broken out by population
# Input: the ExAC v3 VCF
# Output: a log file and a tab-delimited file, rows = variants, columns = populations

import sys
import HTSeq

exacVCF = HTSeq.VCF_Reader("/users/nsabell/data_shared/ExAC/release0.3/ExAC.r0.3.sites.vep.vcf.gz")
#exacVCF = HTSeq.VCF_Reader("/users/nsabell/sardinia/ExAC_TEST.r0.3.sites.vep.vcf.gz")

outFile = open('ExAC_AF_withLoF_byPop.txt', 'w')

annoFields = "Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|HGNC_ID|BIOTYPE|CANONICAL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|EXON|INTRON|DOMAINS|HGVSc|HGVSp|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF"
annoFields = annoFields.strip('">').split('|')

numVars = 0
for var in exacVCF:

	# Read in entry from VCF file
	varTuples = var.info.split(";")
	varInfo = {}

	# Format relevant fields
	for item in varTuples:
		if len(item.split("=")) == 2:
			newItem = item.split("=")
			varInfo[newItem[0]] = newItem[1]
	
	# For each alt allele, read allele counts and numbers into separate dicts
	for i in range(len(var.alt)):
		alleleCounts = {
			'AFR':varInfo['AC_AFR'].split(',')[i],
			'AMR':varInfo['AC_AMR'].split(',')[i],
			'EAS':varInfo['AC_EAS'].split(',')[i],
			'FIN':varInfo['AC_FIN'].split(',')[i],
			'NFE':varInfo['AC_NFE'].split(',')[i],
			'SAS':varInfo['AC_SAS'].split(',')[i],
			'OTH':varInfo['AC_OTH'].split(',')[i],
		}

		alleleNumber = {
			'AFR':varInfo['AN_AFR'],
			'AMR':varInfo['AN_AMR'],
			'EAS':varInfo['AN_EAS'],
			'FIN':varInfo['AN_FIN'],
			'NFE':varInfo['AN_NFE'],
			'SAS':varInfo['AN_SAS'],
			'OTH':varInfo['AN_OTH'],
		}

		# Format ID of variant
		ID_string = str(var.chrom) + "_" + str(var.pos.end - 1) + "_" + str(var.ref) + "_" + str(var.alt[i])

		# Compute allele frequencies
		alleleFrequencies = {}
		for pop in alleleCounts.keys():
			if alleleNumber[pop] == '0':
				alleleFrequencies[pop] = float(0)
			else: 
				alleleFrequencies[pop] = float(alleleCounts[pop]) / float(alleleNumber[pop])

		# Extract loss-of-function annotation, accounting for missing CSQ fields
		lofInfo = []

		if 'CSQ' in varInfo.keys():
			annotations = [dict(zip(annoFields, x.split('|'))) for x in varInfo['CSQ'].split(',')]

			for d in annotations:
				if d['Allele'] == var.alt[i]:
					lofInfo.append(d['LoF'])
		else:
			lofInfo.append('CSQ_missing')

		# Write header if appropriate
		if numVars == 0:
			outFile.write('ID\t' + '\t'.join(str(x) for x in alleleFrequencies.keys()) + '\t' + 'LoF_Info' + '\n')

		# Write tab-delimited frequencies and LoF annotations to output file
		outFile.write(ID_string + '\t' + '\t'.join(str(x) for x in alleleFrequencies.values()) + '\t' + ','.join(str(x) for x in lofInfo) + '\n')
	
	# Increment counter and write to log every 100k VCF lines
	numVars = numVars + 1
	if numVars % 100000 == 0:
		print numVars, " ExAC VCF lines (variant positions) processed."
		sys.stdout.flush()













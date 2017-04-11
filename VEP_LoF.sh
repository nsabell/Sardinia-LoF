#!/bin/bash

# Variables
VEP="/users/nsabell/software/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl"
HUMAN_ANCESTOR_FA="/users/nsabell/software/ensembl-tools-release-86/scripts/variant_effect_predictor/human_ancestor.fa.gz"
PERL_LIB="/users/nsabell/software/perl"

# Command formatting
OUTDIR="/users/nsabell/sardinia/vep_lof_2120"
for chrVCF in /users/nsabell/sardinia/vcfs/ss2120/*vcf.gz; do
	name=`basename $chrVCF | cut -d . -f 1-3`;
	cmd="export PERL5LIB=$PERL_LIB ; $VEP --offline -i $chrVCF -o $OUTDIR/$name.vepLoF.vcf --vcf --stats_text --cache --merged --plugin LoF,human_ancestor_fa:$HUMAN_ANCESTOR_FA > $OUTDIR/$name.vepLoF.log 2>&1 &";
	echo $cmd;
done > vep_lof_2120.cmds


#!/bin/bash

# Variables
VEP="/users/zappala/software/VEP/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl"
HUMAN_ANCESTOR_FA="/users/nsabell/software/ensembl-tools-release-86/scripts/variant_effect_predictor/human_ancestor.fa.gz"
PERL_LIB="/users/nsabell/software/perl:$PERL5LIB"

# Remove older VEP annotations from ExAC
gunzip -c /users/nsabell/data_shared/ExAC/release0.3/ExAC.r0.3.sites.vep.vcf.gz | sed 's/CSQ=.*//' | sed 's/##INFO=<ID=CSQ.*\n//' | gzip > /users/nsabell/sardinia/ExAC.r0.3.sites.noVEP.vcf.gz &

# Command formatting
OUTDIR="/users/nsabell/sardinia/vep_lof_2120"

for chrVCF in /users/nsabell/sardinia/vcfs/ss2120/*vcf.gz; do
	name=`basename $chrVCF | cut -d . -f 1-3`;
	cmd="export PERL5LIB=$PERL_LIB ; $VEP --offline -i $chrVCF -o $OUTDIR/$name.vepLoF.vcf --vcf --stats_text --cache --merged --maf_1kg --maf_exac --maf_esp --plugin LoF,human_ancestor_fa:$HUMAN_ANCESTOR_FA > $OUTDIR/$name.vepLoF.log 2>&1 &";
	echo $cmd;
done > vep_lof_2120.cmds

echo "export PERL5LIB=$PERL_LIB ; $VEP --offline -i /users/nsabell/sardinia/ExAC.r0.3.sites.noVEP.vcf.gz -o $OUTDIR/ExAC.r0.3.sites.vepLoF.vcf --vcf --stats_text --cache --merged --maf_1kg --maf_exac --maf_esp --plugin LoF,human_ancestor_fa:$HUMAN_ANCESTOR_FA > $OUTDIR/ExAC.r0.3.sites.vepLoF.log 2>&1 &" >> vep_lof_2120.cmds

# Execute the commands and zip everything
chmod 777 vep_lof_2120.cmds
./vep_lof_2120.cmds
for i in ./vep_lof_2120/*vcf ; do
	gzip $i;
done
gzip $OUTDIR/ExAC.r0.3.sites.vepLoF.vcf


#!/usr/bin/env bash

# SNP calling, based on a SNP matrix created by Fernando Rabanal and adapted by Elizaveta Grigoreva
# The SNP calling here is done in parallel for each chromosome.

# Input:
#    Directory of aligned BAM files, with corresponding .bai files
#    Reference genome
#    Matrix of SNPs in the parents
# Output:
#    'Targets file' giving SNP positions in the parents
#    five zipped VCF files.


# SLURM
#SBATCH --job-name=call_SNPs
#SBATCH --output=/groups/nordborg/projects/ddm1_resilencing/09_slurm_logs/%x-%a.out
#SBATCH --error=/groups/nordborg/projects/ddm1_resilencing/09_slurm_logs/%x-%a.err
#SBATCH --mem=10GB
#SBATCH --qos=medium
#SBATCH --time=1-20:00:00
#SBATCH --array=1


# === Input files ===

ml build-env/f2022 bcftools/1.17-gcc-12.2.0
ml build-env/f2022 samtools/1.18-gcc-12.3.0

# Directory containing aligned BAM files
workdir=/groups/nordborg/projects/ddm1_resilencing/03_processing
indir=$workdir/004_alignment/all_plates_96w
# Location of the reference genome to map to.
genome=/groups/nordborg/common_data/TAIR10/TAIR10_chromosome_files/TAIR10_chr_all.fas
# SNP matrix file
parental_SNP_matrix=/groups/nordborg/projects/ddm1_resilencing/11_resources/snp_matrix_fernando/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.heterozygous_acc_removed.missingmasked.morethan10pmissingremoved.bcftools.vcf.gz
# Which chromosome to work with
chr=Chr${SLURM_ARRAY_TASK_ID}

# === Output files ===

# Output directory
outdir=$workdir/006_snp_call/bcftools_mpileup
mkdir -p $outdir
# file with location of bamfiles (one per line)
bam_list=${indir}/bam_list.txt
# chromosome specific files
outfile=$outdir/F1_snp_matrix_${chr}.vcf.gz
targets_file=$outdir/targets_file_${chr}.tsv.gz

# === Script ===

# Get genotype likelihoods, and use them to call SNPs
# create targets file
echo "Creating targets file."
#bcftools query -r ${chr} -f'%CHROM\t%POS\t%REF,%ALT\n' $parental_SNP_matrix | bgzip -c > $targets_file
#tabix -s1 -b2 -e2 $targets_file

ulimit -n 10000
# Find bam files
find $indir -name "*Aligned.sortedByCoord.out.bam" > $bam_list
echo "Calling the SNPs."
bcftools mpileup --min-MQ 20 -a FORMAT/DP,FORMAT/AD --skip-indels -f $genome -r $chr -b $bam_list -Ou | \
bcftools call -m --constrain alleles --targets-file $targets_file --variants-only -Oz --output $outfile
tabix $outfile

if [ $? -eq 0 ]
then
    echo "Script completed with exit code 0."
fi

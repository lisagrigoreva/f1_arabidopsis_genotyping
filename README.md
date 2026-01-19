## Arabidopsis genotyping 
This repository contains code to check genotyping for A. thaliana using a set of expected parents from 1001G VCF. 

### Prepare input files
------
To reduce computational time and also accuracy we will work with genic regions only.For this you will need to extract only genic regions from both vcf matricies . Make sure you use updated SNP matrix without heterozygous calls and individuals
TAIR10.gff file should be converted to the TAIR10.bed

```
gff2bed < TAIR10_GFF3_genes.gff | grep 'gene' > TAIR10_GFF3_genes.bed
```
- Using obtained gene bed file prepare 2 matricies: 1001G genic regions and genotype expected genic regions
```
input_progeny=genotyped.samples.snps.variant_non_variant.vcf.gz
output_progeny=genotyped.samples.snps.variant_non_variant.genes.vcf.gz

input_parents=1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.heterozygous_acc_removed.hetmasked.vcf.gz
output_parents=1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.heterozygous_acc_removed.hetmasked.pos_of_interest_genes.vcf.gz

bcftools index $input_progeny 
bcftools view -R regions.bed $input_progeny -Oz -o $output_progeny
# Extract progeny SNPs from 1001G matrix
zcat $output_progeny | cut -f 1,2 >positions_progeny.txt
bcftools view -R positions_progeny.txt $input_parents -Oz -o $output_parents
```

- Convert both files (paternal and progeny ) to the hdf5 format using scikit-allel
```
python
import allel
allel.vcf_to_hdf5('input.vcf', 'output.hdf5', fields='*', overwrite=True)
```

Required input files:
- genic regions in hdf5 format of 1001G
- genic regions in hdf5 format of samples 
- expected parents

| ID in vcf file | Expected parent 1 from 1001G | Expected parent 2 from 1001G |
|----------------|------------------------------|------------------------------|
| 294314         | 6909                         | 7130                         |
| 294315         | 6909                         | 8348                         |
| 294317         | 6909                         | 15592                        |

### Prepare progeny files 
------
- Make files containing sample id, chromosome, position, parent1,parent2, genotype (ref/alt), allelic depth (ref depth/alt depth) and total depth
```
python prepare_progeny.py -i genotypes_progeny.hdf5 -p parents.txt -o output_dir
```
Expected output looks like this:

| Sample_ID | Chromosome | Position | Parent1 | Parent2 | Genotype | Allelic_depth | Total_depth |
|-----------|------------|----------|---------|---------|----------|---------------|-------------|
| 299410    | Chr1       | 7020     | 6909    | 9409    | 0/0      | 7/0           | 8           |
| 299410    | Chr1       | 7035     | 6909    | 9409    | 0/0      | 7/0           | 7           |

 
### Run genotyping and accession probability assignment  

```
python calc_probabilities.py genotypes_parents.hdf5 ./output_folder_previous_step/ 
```


Binomial likelihood test is applied for each parental combination and progeny file. **Expected** transition probabilities for each site are set based on the 1001G crosses. For this, only intersecting positions between progeny and parental 1001G file are considered. For every position expected $\theta$ is reported

- 0/0 results to $\theta=0.0001$
- 0/1 results to $\theta=0.5$
- 1/1 results to $\theta=0.9999$

 
**Observed** transition probabilities are calculated based on the number of reads attributed to the alternative allele (y in equation) in RNASeq data. In VCF file field AD for each site has a format (ref, alt). Thus, alt were taken. 
PMF is calculated for each site as:
```
binom.pmf(alt_counts, total_counts, thetha)
```

Where:

- alt_counts - number of alternative alleles from AD field
- total_counts - depth from DP field
- thetha - thetha for this site in every parental combination 

Finally, the product of probabilites been summed up across all sites in RNASeq data. Additionally, its been normilized to the number of SNPs, because for every progeny file ir differs:

$\frac{1}{u} \sum_i \ln \Pr(y_i \mid \theta)$ 
Where  $u$ is the number of SNPs in each progeny file.
basically:
```
binom.pmf(alt_counts, total_counts, thetha)/u
```
A small PMF value means that the observed alt_counts in the progeny is unlikely given the expected allele frequency ($\theta$)
Normilized likelihood ratios from every parental combination has been attributed for each progeny file.
Then, PMF values were transformed to the probabilities
$\frac{\exp{\lambda_i}}{\sum_i\exp(\lambda)}$

Where exp is PMF and sum is normilized by the number of observed SNPs in each parental-progeny file. 

We are getting to values to focus on:
1) Normilized likelihood ( the closer to 0 the most likely that observed parents are real)
2) Probability. (The higher the probability the most likely that observed parents are real)









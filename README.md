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


## Preprocessing 

### Prepare progeny files 
------
- Make files containing sample id, chromosome, position, parent1,parent2, genotype (ref/alt), allelic depth (ref depth/alt depth) and total depth
```
python 001prepare_progeny.py genotypes_progeny.hdf5 parents.txt output_dir
```
Expected output looks like this:

| Sample_ID | Chromosome | Position | Parent1 | Parent2 | Genotype | Allelic_depth | Total_depth |
|-----------|------------|----------|---------|---------|----------|---------------|-------------|
| 299410    | Chr1       | 7020     | 6909    | 9409    | 0/0      | 7/0           | 7           |
| 299410    | Chr1       | 7035     | 6909    | 9409    | 0/0      | 7/0           | 7           |

 
### Prepare parental files
- Make files containing sample id, chromosome, position and $\theta for each parental combination from the list   
```
python 002prepare_parents.py genotypes_parents.hdf5 parents.txt output_dir
```


## Likelihood assignment 

To assign the most probable parents, a binomial likelihood test is applied for each parental combination and progeny file. Expected allele frequencies ($\theta) for each site are set based on parental genotypes from the 1001G dataset. Only intersecting positions between progeny and parental files are considered for each parental combination.

**Expected** Allele combinations ($\theta)
For each position, the expected $\theta is determined from the parental genotype classification:

- 0/0 results to $\theta=0.0001$
- 0/1 results to $\theta=0.5$
- 1/1 results to $\theta=0.9999$

**Observed** data from progeny
For each site in the progeny RNA-seq data, we observe:
- alt_counts: Number of reads supporting the alternative allele (from AD field)
- total_counts: Total read depth (from DP field)

**PMF Calculation** 
For each site i, the probability of observing the progeny's allele counts given the parental genotype is calculated using the binomial PMF:

binom.pmf(alt_counts, total_counts, $\theta$)

Where:

- alt_counts - number of alternative alleles for the site from AD field
- total_counts - depth for the site (DP field)
- $\theta - Expected allele frequency at site i based on parental genotype

A small PMF value indicates that the observed alt_counts in the progeny is unlikely given the expected allele frequency ($\theta) from that parental pair.
To avoid numerical underflow, PMF is converted to the log-likelihood.

#### Log-Likelihood Calculation

To avoid numerical underflow when multiplying many small probabilities, we work with log-likelihoods.

**Total log-likelihood** for a parental pair:

$$\text{LL} = \sum_{i=1}^{u} \ln(\text{PMF})$$

**Normalized log-likelihood** (to account for different numbers of SNPs across progeny):

$$\text{LL}_{\text{norm}} = \frac{1}{u} \sum_{i=1}^{u} \ln(\text{PMF}_i)$$

Where:
- $u$ = Number of overlapping SNPs between progeny and parental pair
- $\ln(\text{PMF})$ = Natural logarithm of the PMF at site


#### Probability Calculation

To convert log-likelihoods into probabilities that sum to 1 across all tested parental pairs, we use normilized exponential function:

$$P(\text{parent pair } j \mid \text{data}) = \frac{\exp(\text{LL}_j)}{\sum_{k=1}^{m} \exp(\text{LL}_k)}$$

Where:
- $\text{LL}_j$ = Total log-likelihood for parent pair $j$
- $m$ = Total number of parent pairs tested

This gives the posterior probability that each parental pair is the true biological parent, given the observed progeny data.

Output Metrics
For each progeny, we report:

- Probability (%): The posterior probability that a given parental pair is correct
Higher is better (>95% indicates high confidence)
Values close to 100% indicate a clear, confident match
- Normalized Log-Likelihood: The average log-likelihood per SNP
Closer to 0 is better (indicates good fit)
Values around -0.2 to -0.3 indicate good matches
Values < -1.0 indicate poor matches
- Total Log-Likelihood: The cumulative log-likelihood across all SNPs
Used for probability calculation
Accounts for total evidence strength
- Number of SNPs: Number of overlapping positions used in the calculation
More SNPs generally provide more confident assignments








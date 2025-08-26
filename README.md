## Arabidopsis genotyping 
This repository contains code to check genotyping for A. thaliana using a set of expected parents from 1001G matrix. 

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
- genic regions in h5 format of 1001G
- genic regions in h5 format of samples 
- expected parents

| ID in vcf file | Expected parent 1 from 1001G | Expected parent 2 from 1001G |
|----------------|------------------------------|------------------------------|
| 294314         | 6909                         | 7130                         |
| 294315         | 6909                         | 8348                         |
| 294317         | 6909                         | 15592                        |


### Run genotyping
------
- Make files with depth and alleles for every site for expected parents 
```
python prepare_genotypes.py -i genotypes.hdf5 -p parents.txt -o output_dir
```
- Make files with expected parents pairs and add thetha for every site. If the site is homozygous (0/0 or 1/1), θ=0+0.0001; if the site is heterozygous, θ is set to 0.5
  
```
python parent_classifier.py genotypes.hdf5 parents.txt output_dir
```
- Create all pairwise combinations of parental-accessions pairs and add PMF to every intersection 
```
python create_intersections_with_pmf.py  genotypes_progeny_output_dir/ all_snps_classified_parents_output_dir/ intersections/
```





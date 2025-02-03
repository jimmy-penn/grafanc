# GrafAnc Software Documentation

GrafAnc is a software tool that can be used for genetic ancestry inference.

A C++ executable `grafanc` and other files are included in the package `GrafAnc.tar.gz` and visible as separate files after the user executes the command:

``` sh
tar zxvf GrafAnc.tar.gz
```

GrafAnc calculates genetic distances from each human individual (or sample) to some reference populations and estimates subject ancestry and ancestral proportions based on these distances.

Current version of GrafAnc takes a genotype data set as input and generates 18 scores to infer genetic ancestry for each individual, at both continental and subcontinental levels. The three GD scores, GD1, GD2, and GD3, are used in inferring genetic ancestry at continental level. GrafAnc assumes that each individual is an admixture of three continental ancestries: European (E), African (F), and Asian (A), and estimates ancestral proportions *P*<sub>e</sub>, P<sub>f</sub>, P<sub>a</sub> based on GD1 and GD2 scores using barycentric coordinates. Other scores are used to infer genetic ancestry at sub-continental levels.

### Input files

GrafAnc takes genotype data sets in either PLINK format (`.fam`, `.bim`, `.bed`) or VCF format (`.vcf` or `.vcf.gz`). Variants can be labeled using either RS IDs, and/or chromosome positions in either GRCh37 or GRCh38. Users don't need to specify which genome build to use. GrafAnc will figure out how to read the variant information.

GrafAnc uses 282,424 pre-selected biallelic SNPs (called ancestry SNPs in this document) to do ancestry inference. It is ideal that all the ancestry SNPs are included in the input file. If some of the ancestry SNPs are missing, GrafAnc still makes unbiased estimation of genetic ancestry, but with lower accuracy, depending on how many ancestry SNPs are missing. Usually at least 1,000 ancestry SNPs with non-missing genotypes are need to infer genetic ancestry at continental level, and at least 10,000 SNPs for subcontinental ancestry inference. GrafAnc will not generate any ancestry inference results if less than 100 ancestry SNPs are included in the input file.

GrafAnc infers genetic ancestry for each individual in the input file, one by one. It is acceptable to included duplicated samples or closely related individuals in the input file. It is also OK to split the file into multiple data sets, and combine the results into one file after running GrafAnc on all data sets. The results can be combined even if they were generated from genotype data sets collected using different sequencing or genotyping methods, and including different sets of variants.

GrafAnc tolerates genotype data missingness. SNPs and individuals don't need to be filtered out because of genotype data missingness, unless the user thinks that the missingness might reflect low data qualities for these individuals or SNPs

Linkage disequilibrium has already been considered when the ancestry SNPs were selected. So usually no LD pruning is needed, unless the user believes some of the ancestry SNPs need to be excluded for some stricter linkage equilibrium requirements.

### Running `grafanc` to infer genetic ancestry

`grafanc` is a command line executable that can be run under GNU/LINUX 64 bit systems. Brief usage is given when the program is executed without parameters:

``` sh
$ grafanc

Usage: grafanc <Binary PLINK set or VCF file> <output file>
        --maxmem  <size in MB>:  specify maximum memory in MiB to be used by GrafAnc. Default 8 MB
        --threads <number>:      specify maximum number of threads to use. Default 1
        --samples <number>:      specify maximum samples to be processed in each round
```

By default GrafAnc uses at most 8 MB memory. The user can specify maximum memory (in MB) to use using `--maxmem` option.

GrafAnc runs with one thread unless the user specify a different number of threads to user using `--threads` option.

The `--samples` option may be needed when the input data set is very large and there is not enough memory to process all the samples at once. To reduce memory usage, GrafAnc will split the data set into small ones, each of which contains a small set of samples. GrafAnc calculates how many samples need to be included in each data set. However, the user can specify how many samples should be included in each set using the `--samples` option when necessary.

The following command determines population structures and saves results to the output file:

``` sh
$ grafanc data/TGP.bed results/TGP_pops.txt 
```

The SNPs in the PLINK set are entered as RS IDs. If the SNP IDs are not provided, or provided but not in RS IDs, it is acceptable to `grafanc` if GRCh37 or GRCh38 chromosome positions are included in the genotype data set, e.g.,

``` sh
$ grafanc data/TGP_40_50K_gb37.bed results/TGP_40_50K_gb37_pops.txt
$ grafanc data/TGP_40_50K_gb38.bed results/TGP_40_50K_gb38_pops.txt
```

The input data set can also be a VCF file (zipped or not), e.g.,

``` sh
$ grafanc data/TGP_50_10K.vcf.gz results/TGP_50_10K_pops.txt 
```

If the VCF file includes many more SNPs than those being used by GrafAnc, e.g., containing whole genome sequencing data, `grafanc` can still read the data and do ancestry inference. However, it is recommend that the Perl script `ExtractAncSnpsFromVcfGz.pl` be used to extract the genotypes before `grafanc` is run (see instructions below for usage of the Perl script).

`grafanc` and the Perl scripts included in the package can be called from other directories, e.g.,

``` sh
$ cd results
$ ../grafanc ../data/TGP_50_10K_gb37_flipAlt5.vcf.gz TGP_50_10K_flip_pops.txt
```

GrafAnc C++ executable and Perl scripts need to find some information included in the `data` directory when being run. If the executable is moved away from the `data` directory, the user can set environment variable `GRAFPATH` to the directory where GrafAnc `data` directory and Perl packages (`.pm` files) are located and call `grafanc` and Perl scripts from any location.

### Running `ExtractAncSnpsFromVcfGz.pl` to extract genotypes of ancestry SNPs from one or more VCF files

When the VCF file contains many more SNPs than those used by GrafAnc, or genotypes are saved into multiple files, e.g., one file for one chromosome, the user can use software tools such as `bcftools` to extract genotypes of ancestry SNPs and save them into one file and then pass it to `grafanc`. The 282,424 ancestry SNPs can be found in file `AncSnpPopAFs.txt`, in which the first column is the chromosome, second and third columns are GRCh37 and GRCh38 positions, and the fourth one is the RS ID.

The user can also use the script `ExtractAncSnpsFromVcfGz.pl` (included in the package) to extract genotypes of ancestry SNPs from VCF files for `grafanc`.

`ExtractAncSnpsFromVcfGz.pl` takes two required parameters. The first parameter specifies the input VCF file, and second one is the output VCF file.

``` sh
$ ExtractAncSnpsFromVcfGz.pl
Usage: ExtractAncSnpsFromVcfGz.pl <vcf_or_vcf_gz_file> <output_vcf_file> [keyword]
```

For example, one can run the following commands to extract genotypes of ancestry SNPs and save the data to the output file, then pass the file to `grafanc` for ancestry inference:

``` sh
$ ExtractAncSnpsFromVcfGz.pl data/TGP_10_chr4.vcf.gz results/TGP_10_chr4_anc.vcf
$ grafanc results/TGP_10_chr4_anc.vcf results/TGP_10_chr4_pops.txt
```

If the genotypes of the same set of samples are saved in a set of VCF files, e.g., one file for one chromosome, `ExtractAncSnpsFromVcfGz.pl` can find all these files and extract genotypes of ancestry SNPs and save the results into one VCF file, so that it can be used by `grafanc`. If all file names differ from one another only by an embedded integer, one can run `ExtractAncSnpsFromVcfGz.pl` with the "keyword" before the integer as the optional third parameter to let the script extract genotypes from all these VCF files, e.g.,

``` sh
$ ExtractAncSnpsFromVcfGz.pl data/TGP_10_chr4.vcf.gz results/TGP_10_cmb.vcf chr
$ grafanc results/TGP_10_cmb.vcf results/TGP_10_cmb_pops.txt
```

The above command extracts genotypes of ancestry SNPs from the four files with name like "TGP_10_chr\<`number>`.vcf.gz" and save the results to the output file "TGP_10_cmb.vcf".

### Output files

The output file generated by `grafanc` is a table containing statistic scores and other information for all samples (in rows) in the input data set. The table contains the following columns:

| Column              | Description                                                                                |
|------------------------|------------------------------------------------|
| Sample              | Sample ID                                                                                  |
| #SNPs               | Number of ancestry SNPs with non-missing genotypes                                         |
| GD1, GD2, GD3       | Scores for inferring genetic ancestry at continental level                                 |
| EA1, EA2, EA3, EA4  | Scores designed to infer population structure for East Asians                              |
| AF1, AF2, AF3       | For Africans                                                                               |
| EU1, EU2, EU3       | For Europeans                                                                              |
| SA1, SA2            | For South Asians                                                                           |
| AF1, AF2, AF3       | For Africans                                                                               |
| IC1, IC2, IC3       | For inter-continental populations                                                          |
| Pe, Pf, Pa          | Ancestry components estimated by GrafAnc (percentages of European, African and East Asian) |
| RawPe, RawPf, RawPa | Raw Pe, Pf, Pa values, continuous and not bounded between 0 and 1                          |
| AncGroupID          | Ancestry group ID assigned by GrafAnc                                                      |

GrafAnc assigns ancestry groups at two levels, similar to the self-reported populations in 1000 Genomes Project (1KGP), and in UK Biobank (UKBB). Continental level ancestry information (or super-populations like in 1KGP) is coded in the first digit of `AncGroupID`:

| AncGroupID | 3-letter code | Description                  |
|------------|---------------|------------------------------|
| 100        | AFR           | African                      |
| 200        | MEN           | Middle East and North Africa |
| 300        | EUR           | European                     |
| 400        | SAS           | South Asian                  |
| 500        | EAS           | East Asian                   |
| 600        | AMR           | Admixed American             |
| 700        | OCN           | Oceania                      |
| 800        | MIX           | Multi-ancestry               |

The other two digits in `AncGroupID` code for sub-continental ancestry groups. See file `AncGroupIdCodeDesc.xlsx` (included in the package) for code descriptions. The ancestry group assignment was made based on the reference data from UKBB, 1KGP, and Human Genome Diversity Project (HGDP), especially the "Country of Birth" phenotype from UKBB. Please note that there are no clear boundaries between human populations. No ancestry group assigned by GrafAnc is expected to include only those individuals with ancestry backgrounds described in the above two tables. Actually in the tables we only provide our best estimation of which subpopulations are most likely to be included in each ancestry group assigned by GrafAnc.

When assigning ancestry groups, GrafAnc assumes there are enough (e.g., 10,000+ SNPs for continental, and 5,000+ SNPs for subcontinental populations) 50ancestry SNPs with non-missing genotypes in the input data set. When there are fewer SNPs with genotypes, the ancestry group assignments are less reliable. In order to estimate the variance of the GrafAnc scores when difference number of genotyped ancestry SNPs are available, we tested GrafAnc using simulated data and saved results to file `<tobeadded>.xlsx`.

### Visualization of GrafAnc results

GrafAnc scores can be visualized using 2D or 3D scatterplots. Although any combinations of GrafAnc scores can be plotted on a scatterplot, it is suggested that the user try the following combinations first:

| X-axis | Y-axis | Purpose                                                         |
|------------------|------------------|-------------------------------------|
| GD1    | GD2    | Global populations                                              |
| GD1    | IC1    | Global population, separating South Asians from Latin Americans |
| GD1    | GD3    | Latin Americans, Native Americans, also for data QC             |
| EU1    | EU2    | European subpopulations                                         |
| EU1    | EU3    | European subpopulations                                         |
| EU1    | IC2    | Middle East and North Africa individuals                        |
| EU1    | IC3    | Middle East and North Africa individuals                        |
| AF1    | AF2    | African subpopulations                                          |
| AF1    | AF3    | African subpopulations                                          |
| EA1    | EA2    | East Asian subpopulations                                       |
| EA1    | EA3    | East Asian and North Asian subpopulations                       |
| EA1    | EA4    | East Asian and Southeast Asian subpopulations                   |
| SA1    | SA2    | South Asian subpopulations                                      |

### UKBB GrafAnc results used as references

GrafAnc results obtained using different data sets can be combined and plotted on the same scatterplot. We combined the results obtained using UKBB, 1KGP and HGDP and calculated the mean values for each subcontinental population, after removing outliers. The results are saved to `PopMeanScores.txt (to be added)`. Users can combine these results to those obtained using their own data, no matter what variants are included.

## References

Jin Y, Schäffer AA, Feolo M, Holmes JB and Kattman BL (2019). [GRAF-pop: A Fast Distance-based Method to Infer Subject Ancestry from Multiple Genotype Datasets without Principal Components Analysis.](https://www.g3journal.org/content/9/8/2447.long) G3: Genes \| Genomes \| Genetics. August 1, 2019 vol. 9 no. 8 2447-2461.

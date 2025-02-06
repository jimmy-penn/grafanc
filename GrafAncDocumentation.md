# GrafAnc Software Documentation

GrafAnc calculates genetic distances from each human individual (or sample) to some reference populations and estimates subject ancestry and ancestral proportions based on these distances.

Current version of GrafAnc takes a genotype data set as input and generates 18 scores to infer genetic ancestry for each individual, at both continental and subcontinental levels. The three GD scores, GD1, GD2, and GD3, are used in inferring genetic ancestry at continental level. GrafAnc assumes that each individual is an admixture of three continental ancestries: European (E), African (F), and Asian (A), and estimates ancestral proportions *P*<sub>*e*</sub>, *P*<sub>*f*</sub>, *P*<sub>*a*</sub> based on GD1 and GD2 scores using barycentric coordinates. Other scores are used to infer genetic ancestry at sub-continental levels.

### Input files

GrafAnc takes genotype data sets in either PLINK format (`.fam`, `.bim`, `.bed`) or VCF format (`.vcf` or `.vcf.gz`). Variants can be labeled using either RS IDs, and/or chromosome positions in either GRCh37 or GRCh38. Users don't need to specify which genome build to use. GrafAnc will figure out how to read the variant information.

GrafAnc uses 282,424 pre-selected biallelic SNPs (called ancestry SNPs in this document) to do ancestry inference. It is ideal that all the ancestry SNPs are included in the input file. If some of the ancestry SNPs are missing, GrafAnc still makes unbiased estimation of genetic ancestry, but with lower accuracy, depending on how many ancestry SNPs are missing. Usually in order to infer population structure confidently, at least 10,000 ancestry SNPs with non-missing genotypes are needed for continental populations, and at least 50,000 SNPs for subcontinental ones. GrafAnc will not generate any ancestry inference results if less than 100 ancestry SNPs are included in the input file.

GrafAnc infers genetic ancestry for each individual in the input file, independent of other individuals. It is acceptable to included duplicated samples or closely related individuals in the input file. It is also OK to split the file into multiple data sets and combine the results into one file after running GrafAnc on all data sets. The results can be combined even if they were generated from genotype data sets including different variants.

GrafAnc tolerates genotype data missingness. SNPs and individuals don't need to be filtered out because of genotype data missingness, unless the user thinks that the missingness might reflect low data qualities for these individuals or SNPs.

Linkage disequilibrium has already been considered when the ancestry SNPs were selected. So usually no LD pruning is needed, unless the user believes some of the ancestry SNPs need to be excluded for some stricter linkage equilibrium requirements.

### Running `grafanc` to infer genetic ancestry

`grafanc` is a command line executable that can be run under GNU/LINUX 64 bit systems. Brief usage is given when the program is executed without parameters:

``` sh
$ grafanc

Usage: grafanc <Binary PLINK set or VCF file> <output file>
        --maxmem  <size in MB>:  specify maximum memory in MB to be used by GrafAnc. Default 8 MB
        --threads <number>:      specify maximum number of threads to use. Default 1
        --samples <number>:      specify maximum samples to be processed in each round
```

By default GrafAnc uses at most 8 MB memory. The user can specify maximum memory (in MB) using `--maxmem` option.

GrafAnc runs with one thread unless the user specifies a different number of threads using `--threads` option.

The `--samples` option may be needed when the input data set is very large and there is not enough memory to process all the samples at once. To reduce memory usage, GrafAnc will split the data set into small ones, each of which contains a small set of samples. GrafAnc calculates how many samples need to be included in each data set. However, the user can specify how many samples should be included in each set using the `--samples` option when necessary.

The following command determines population structures and saves results to the output file:

``` sh
$ cd cpp
$ mkdir results
$ grafanc ../medium_testing/input/TGP_40_50K_rs.bed data/TGP.bed results/TGP_40_50K_pops.txt 
```

The SNPs in the PLINK set are entered as RS IDs. If the SNP IDs are not provided, or provided but not in RS IDs, it is acceptable to `grafanc` if GRCh37 or GRCh38 chromosome positions are included, e.g.,

``` sh
$ grafanc ../medium_testing/input/TGP_40_50K_gb37.bed results/TGP_40_50K_gb37_pops.txt
$ grafanc ../medium_testing/input/TGP_40_50K_gb38.bed results/TGP_40_50K_gb38_pops.txt
```

The input data set can also be a VCF file (zipped or not), e.g.,

``` sh
$ grafanc ../medium_testing/input/TGP_50_10K.vcf.gz results/TGP_50_10K_pops.txt 
```

If the VCF file includes many more SNPs than those being used by GrafAnc, e.g., containing whole genome sequencing data, `grafanc` can still read the data and do ancestry inference. However, it is recommend that the Perl script `ExtractAncSnpsFromVcfGz.pl` be used to extract the genotypes before `grafanc` is run (see instructions below for usage of the Perl script).

`grafanc` and the Perl scripts can be called from other directories, e.g.,

``` sh
$ cd results
$ ../grafanc ../../medium_testing/input/TGP_50_10K_gb37_flipAlt5.vcf.gz TGP_50_10K_flip_pops.txt
```

GrafAnc C++ executable and Perl scripts need to find some information included in the `data` directory when being run. If the executable is moved away from the `data` directory, the user can set environment variable `GRAFPATH` to the directory where GrafAnc `data` directory and Perl packages (`.pm` files) are located and call `grafanc` and Perl scripts from any directory.

### Running `ExtractAncSnpsFromVcfGz.pl` to extract genotypes of ancestry SNPs from one or more VCF files

When the VCF file contains many more SNPs than those used by GrafAnc, or genotypes are saved into multiple files, e.g., one file for one chromosome, the user can use software tools such as `bcftools` to extract genotypes of ancestry SNPs and save them into one file and then pass it to `grafanc`. The 282,424 ancestry SNPs can be found in file `cpp/data/AncSnpPopAFs.txt`, in which the first column is the chromosome, second and third columns are GRCh37 and GRCh38 positions, and the fourth one is the RS ID.

The user can also use the script `perl/ExtractAncSnpsFromVcfGz.pl` to extract genotypes of ancestry SNPs from VCF files for `grafanc`.

`ExtractAncSnpsFromVcfGz.pl` takes two required parameters. The first parameter specifies the input VCF file, and second one is the output VCF file.

``` sh
$ ExtractAncSnpsFromVcfGz.pl
Usage: ExtractAncSnpsFromVcfGz.pl <vcf_or_vcf_gz_file> <output_vcf_file> [keyword]
```

For example, one can run the following commands to extract genotypes of ancestry SNPs and save the data to the output file, then pass the file to `grafanc` for ancestry inference:

``` sh
$ ExtractAncSnpsFromVcfGz.pl ../medium_testing/input/TGP_10_chr4.vcf.gz results/TGP_10_chr4_anc.vcf
$ grafanc results/TGP_10_chr4_anc.vcf results/TGP_10_chr4_pops.txt
```

If the genotypes of samples are saved in a set of VCF files, e.g., one file for one chromosome, `ExtractAncSnpsFromVcfGz.pl` can find all these files and extract genotypes of ancestry SNPs and save the results into one VCF file, so that it can be used by `grafanc`. If all file names differ from one another only by an embedded integer, one can run `ExtractAncSnpsFromVcfGz.pl` with the "keyword" before the integer as the optional third parameter to let the script extract genotypes from all these VCF files, e.g.,

``` sh
$ ExtractAncSnpsFromVcfGz.pl ../medium_testing/input/TGP_10_chr4.vcf.gz results/TGP_10_cmb.vcf chr
$ grafanc results/TGP_10_cmb.vcf results/TGP_10_cmb_pops.txt
```

The above command extracts genotypes of ancestry SNPs from the four files with name like "TGP_10_chr\<`number>`.vcf.gz" and save the results to the output file `TGP_10_cmb.vcf`.

### Output files

The output file generated by `grafanc` is a table containing statistic scores and other information for all samples (in rows) from the input data set. The table contains the following columns:

| Column              | Description                                                                                |
|----------------------|--------------------------------------------------|
| Sample              | Sample ID                                                                                  |
| #SNPs               | Number of ancestry SNPs with non-missing genotypes                                         |
| GD1, GD2, GD3       | Scores for inferring genetic ancestry at continental level                                 |
| EA1, EA2, EA3, EA4  | Scores designed to infer population structure for East Asians                              |
| AF1, AF2, AF3       | For Africans                                                                               |
| EU1, EU2, EU3       | For Europeans                                                                              |
| SA1, SA2            | For South Asians                                                                           |
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

The other two digits in `AncGroupID` code for sub-continental ancestry groups. The following table shows the ancestry backgrounds of individuals assigned to each ancestry group by GrafAnc. `Main regions` are the regions of ancestry origins or ethnic groups of most of the individuals in each ancestry group. The last column includes the subcontinental population groups in UKBB, 1KGP and HGDP that are assigned to each ancestry group based on the *average* GrafAnc scores. In this column, HGDP and 1KGP populations are prefixed with `HGDP-` and `TGP-`, respectively. Those without `-` are UKBB `Country of Birth` groups.

| AncGroupID | Main regions        | Countries or ethnic groups in UKBB, 1KGP, HGDP                                                                                                        |
|------------------|------------------|------------------------------------|
| 101        | Nigeria             | Nigeria, HGDP-Yoruba, TGP-ESN, TGP-YRI                                                                                                                |
| 102        | Western Africa      | Gambia, Ghana, Guinea, Ivory Coast, Liberia, Senegal, Sierra Leone, Togo, HGDP-Mandenka, TGP-GWD, TGP-MSL                                             |
| 103        | Central Africa      | Angola, Burundi, Cameroon, Congo, Malawi, Rwanda, Tanzania, Uganda, Zambia                                                                            |
| 104        | Kenya               | Kenya, HGDP-Bantu Kenya, HGDP-Biaka, TGP-LWK                                                                                                          |
| 105        | Southern Africa     | Botswana, South Africa, Zimbabwe, HGDP-Bantu South Africa, HGDP-Mbuti                                                                                 |
| 106        | Other Africa        | Mozambique                                                                                                                                            |
| 201        | Northern Africa     | Algeria, Libya, Morocco, Tunisia, HGDP-Mozabite                                                                                                       |
| 202        | Northeastern Africa | Eritrea, Ethiopia, Somalia, Sudan                                                                                                                     |
| 203        | Middle East 1       | Cyprus, Egypt, Israel, Jordan, Kuwait, Lebanon, Palestine, Syria, Yemen, HGDP-Bedouin, HGDP-Druze, HGDP-Palestinian                                   |
| 204        | Middle East 2       | Armenia, Azerbaijan, Iran, Iraq, Turkey                                                                                                               |
| 301        | Finland             | Finland, TGP-FIN                                                                                                                                      |
| 302        | Northern Europe     | Denmark, Iceland, Norway, Sweden                                                                                                                      |
| 303        | Western Europe      | Austria, Belgium, British, France, Germany, Irish, Luxembourg, Netherlands, Switzerland, United Kingdom, HGDP-French, HGDP-Orcadian, TGP-CEU, TGP-GBR |
| 304        | Southern Europe     | Italy, Portugal, Spain, HGDP-Basque, HGDP-Bergamo Italian, HGDP-Sardinian, HGDP-Tuscan, TGP-IBS, TGP-TSI                                              |
| 305        | Northeastern Europe | Estonia, Latvia, Lithuania, Poland, Russia, Ukraine, HGDP-Russian                                                                                     |
| 306        | Southeastern Europe | Bulgaria, Croatia, Czech, Hungary, Macedonia, Romania, Serbia, Slovakia, Slovenia                                                                     |
| 307        | Balkans             | Albania, Greece, Kosovo                                                                                                                               |
| 308        | Other Europe        | HGDP-Adygei                                                                                                                                           |
| 401        | Asian India         | Mauritius, TGP-ITU, TGP-PJL                                                                                                                           |
| 402        | Gujarati in India   | TGP-GIH                                                                                                                                               |
| 403        | Northern South Asia | Afghanistan, India, Kashmir, Pakistan, HGDP-Balochi, HGDP-Brahui, HGDP-Burusho, HGDP-Kalash, HGDP-Makrani, HGDP-Pathan, HGDP-Sindhi                   |
| 404        | Sri Lanka           | Sri Lanka, TGP-STU                                                                                                                                    |
| 405        | Bangladesh          | Bangladesh, TGP-BEB                                                                                                                                   |
| 501        | Japan Ryukyu        |                                                                                                                                                       |
| 502        | Japan Main Islands  | Japan, HGDP-Japanese, TGP-JPT                                                                                                                         |
| 503        | Korea               | North Korea, South Korea                                                                                                                              |
| 504        | Northern Asia       | Kyrgyzstan, Mongolia, Nepal, HGDP-Hazara, HGDP-Uygur, HGDP-Yakut                                                                                      |
| 505        | Northern China 1    | China, HGDP-Northern Han, TGP-CHB                                                                                                                     |
| 506        | Northern China 2    | HGDP-Daur, HGDP-Hezhen, HGDP-Mongolian, HGDP-Oroqen, HGDP-Xibo                                                                                        |
| 507        | Southern China 1    | Brunei, Hong Kong, Indonesia, Malaysia, Singapore, Taiwan, HGDP-Han, HGDP-Miao, HGDP-She, HGDP-Tujia, TGP-CHS                                         |
| 508        | Southern China 2    | HGDP-Naxi, HGDP-Tu, HGDP-Yi                                                                                                                           |
| 509        | Southeast Asia      | Cambodia, Philippines, Vietnam , HGDP-Dai, HGDP-Lahu, TGP-CDX, TGP-KHV                                                                                |
| 510        | Thailand            | Laos, Thailand, HGDP-Cambodian                                                                                                                        |
| 511        | Other East Asia     |                                                                                                                                                       |
| 601        | African American    | Antigua, Aruba, Barbados, Caribbean, Guianas, West Indies, TGP-ACB, TGP-ASW                                                                           |
| 602        | Latin American 1    | Brazil, TGP-PUR                                                                                                                                       |
| 603        | Latin American 2    | Argentina, Bolivia, Chile, Colombia, Ecuador, Mexico, Peru, Venezuela, TGP-CLM, TGP-MXL                                                               |
| 604        | Native American     | HGDP-Colombian, HGDP-Karitiana, HGDP-Maya, HGDP-Pima, HGDP-Surui, TGP-PEL                                                                             |
| 700        | Oceania             | HGDP-Bougainville, HGDP-Papuan Highlands                                                                                                              |
| 800        | Multiple Races      |                                                                                                                                                       |

The ancestry group assignment was made based on the reference data from UKBB, 1KGP, and Human Genome Diversity Project (HGDP), especially the "Country of Birth" phenotype from UKBB. Please note that there are no clear boundaries between human populations. No ancestry group assigned by GrafAnc is expected to include only those individuals with ancestry backgrounds described in the above two tables. What we provide in the tables are our best estimation of which subpopulations are most likely to be included in each ancestry group assigned by GrafAnc.

### Estimating reliabilities of GrafAnc results

When assigning ancestry groups, GrafAnc assumes there are enough (see above) ancestry SNPs with non-missing genotypes in the input data set. When there are fewer SNPs with genotypes, the ancestry group assignments are less reliable. In order to estimate the variances of the GrafAnc scores when difference number of genotyped ancestry SNPs are available, using 1KGP data we tested GrafAncSouth by randomly setting some genotypes as missing. As expected, the standard deviations of a GrafAnc scores are directly proportional to the inverse of square root of numbers of SNPs with non-missing genotypes. We have done linear regression analysis and obtained the parameters for all GrafAnc scores:

| Score | *a*     | *b*   |
|-------|---------|-------|
| GD1   | -0.0001 | 0.83  |
| GD2   | -0.0001 | 0.95  |
| GD3   | -0.0002 | 1.31  |
| EA1   | -0.0023 | 17.94 |
| EA2   | -0.0021 | 20.68 |
| EA3   | -0.0018 | 12.74 |
| EA4   | -0.0038 | 15.30 |
| AF1   | -0.0035 | 30.08 |
| AF2   | -0.0025 | 46.89 |
| AF3   | -0.0020 | 18.11 |
| EU1   | -0.0028 | 17.74 |
| EU2   | -0.0014 | 14.56 |
| EU3   | -0.0016 | 22.42 |
| SA1   | -0.0030 | 21.85 |
| SA2   | -0.0009 | 12.01 |
| IC1   | -0.0008 | 6.83  |
| IC2   | -0.0008 | 10.95 |
| IC3   | -0.0017 | 15.99 |

For each score calculated from a sample, let *n* be the number of SNPs with non-missing ancestry SNPs, the standard deviation $\sigma$ can be estimated using the following equation:

$\sigma$ = *a* + *b* / $\sqrt{n}$

where parameters *a* and *b* can be found from the above table.

### Visualization of GrafAnc results

GrafAnc scores can be visualized using 2D or 3D scatterplots. Although any combinations of GrafAnc scores can be plotted on a scatterplot, it is suggested that the user try the following combinations first:

| X-axis | Y-axis | Purpose                                                         |
|------------------|------------------|------------------------------------|
| GD1    | GD2    | Global populations                                              |
| GD1    | IC1    | Global population, separating South Asians from Latin Americans |
| GD1    | GD3    | Latin Americans, Native Americans, also for data QC             |
| EU1    | EU2    | European subpopulations                                         |
| EU1    | EU3    | European subpopulations                                         |
| EU1    | IC2    | Middle East and North Africa separated from other populations   |
| EU1    | IC3    | Middle East and North Africa separated from other populations   |
| AF1    | AF2    | African subpopulations                                          |
| AF1    | AF3    | African subpopulations                                          |
| EA1    | EA2    | East Asian subpopulations                                       |
| EA1    | EA3    | East Asian and North Asian subpopulations                       |
| EA1    | EA4    | East Asian and Southeast Asian subpopulations                   |
| SA1    | SA2    | South Asian subpopulations                                      |

Ancestry proportions *P*<sub>*e*</sub>, *P*<sub>*f*</sub>, *P*<sub>*a*</sub> are calculated based on the GD1 and GD2 values using Barycentric cooridinates, referenced to the triangle formed by the centriods of European, African and East Asian (see the following table).

| Population     | GD1    | GD2    |
|----------------|--------|--------|
| European (E)   | 1.4758 | 1.4370 |
| African (F)    | 1.0800 | 1.1000 |
| East Asian (A) | 1.7045 | 1.1000 |

In many cases, it is useful to plot the reference triangle on the GD1-GD2 scatterplots.

### Using UKBB, 1KGP, HGDP GrafAnc results as references

GrafAnc results obtained using different data sets can be combined and plotted on the same scatterplot. We combined the results obtained using UKBB, 1KGP and HGDP and calculated the mean values for each subcontinental population, after removing outliers. The results are saved to [PopMeanScores.txt](https://github.com/jimmy-penn/grafanc/blob/master/PopMeanScores.txt). Users can combine these results to those obtained using their own data, no matter what variants are included. Note that some of the scores in the table are left empty, which is because these populations were used as training sets for calculating these scores.

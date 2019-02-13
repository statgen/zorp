tail -n +2 has-fields-.txt | sort -k1,1n -k2,2n | bgzip > sample.gz

tabix sample.gz -s1 -b2 -e2



## Association file datasources
EPACTS, SAIGE, RAREMETAL, rvtests, many many many others




- [Plink](http://zzz.bwh.harvard.edu/plink/anal.shtml#glm) (1.0.7), linear association examples
     CHR       Chromosome
     SNP       SNP identifier
     BP        Physical position (base-pair)
     A1        Tested allele (minor allele by default) 
     TEST      Code for the test (see below)
     NMISS     Number of non-missing individuals included in analysis
     BETA/OR   Regression coefficient (--linear) or odds ratio (--logistic)
     STAT      Coefficient t-statistic 
     P         Asymptotic p-value for t-statistic

- EPACTS

- RAREMETAL (TODO: what we really want is raremetworker for the non-group gwas results)
    Header lines begins with '#' or 'CHROM'. After header part, we listed the meaning of each column as following:

    CHROM Chromosome
    POS Position
    REF Reference Allele
    ALT Alternative Allele
    N_INFORMATIVE Count of individuals with genotype and phenotype
    FOUNDER_AF(RAREMETALWORKER only) Allele frequency among founders
    ALL_AF(RAREMETALWORKER only) Allele frequency across entire sample
    AF(rvtests only) Allele frequency (for related samples, this is adjusted allele frequency)
    INFORMATIVE_ALT_AC Copies of the rare allele among samples with genotype and phenotype
    CALL_RATE Fraction of called genotypes
    HWE_PVALUE Exact Hardy-Weinberg equilibrium p-value
    N_REF Count of reference homozygotes
    N_HET Count of heterozygotes
    N_ALT Count of alternative allele homozygotes
    U_STAT Score statistic numerator
    SQRT_V_STAT Score statistic denominator
    ALT_EFF_SIZE Estimated effect size
    PVALUE P-value 


    https://genome.sph.umich.edu/wiki/RAREMETAL_Documentation#Gene-level_Tests_Meta-Analysis_Output - two kinds
    Long tables
 ##Method=Burden
 ##STUDY_NUM=2
 ##TotalSampleSize=14308
 #GROUPNAME      NUM_VAR VARs    MAFs    SINGLEVAR_EFFECTs       SINGLEVAR_PVALUEs       AVG_AF  MIN_AF  MAX_AF  EFFECT_SIZE     PVALUE
 NOC2L   7       1:880502:C:T;1:881918:G:A;1:887799:C:T;1:888659:T:C;1:889238:G:A;1:891591:C:T;1:892380:G:A        0.000166722,0.0242172,0.0109203,0.0355845,0.0333729,0.00700233,0.00200067       -0.183575,-0.00228307,-0.0598337,0.0220595,0.0229464,-0.0302768,-0.0200417      0.790161,0.953446,0.515806,0.503548,0.499251,0.791773,0.926625  0.0161807       0.000166722     0.0355845       0.00667875      0.662531
 KLHL17  2       1:897285:A:G;1:898869:C:T       0.0148408,0.00108369    -0.0502034,-0.0256403   0.528269,0.934606       0.00796222      0.00108369      0.0148408       -0.0484494      0.528878

    Short tables
 ##Method=Burden
 ##STUDY_NUM=2
 ##TotalSampleSize=14308
 #GROUPNAME      NUM_VAR VARs    AVG_AF  MIN_AF  MAX_AF  EFFECT_SIZE     PVALUE
 NOC2L   7       1:880502:C:T;1:881918:G:A;1:887799:C:T;1:888659:T:C;1:889238:G:A;1:891591:C:T;1:892380:G:A      0.0161807       0.000166722     0.0355845       0.00667875      0.662531
 KLHL17  2       1:897285:A:G;1:898869:C:T       0.00796222      0.00108369      0.0148408       -0.0484494      0.528878


- Rvtests "Explanation of outputs"- https://github.com/zhanxw/rvtests#models (no example??)
    https://genome.sph.umich.edu/wiki/Summary_Statistics_Files_Specification_for_RAREMETAL_and_rvtests    
    N_INFORMATIVE: Number of samples that are analyzed for association.
    AF: allele frequency. For related individuals, we use BLUE estimator. For case-control study, we list overall frequency (adjusted by relatedness if possible), case frequency and control frequency separated by a colon.
    INFORMATIVE_ALT_AC: The number of alternative alleles in the analyzed samples.
    CALL_RATE: The fraction of non-missing alleles. For case-control study, we calculate call rate for all samples, case samples and control samples separated by a colon.
    HWE_PVALUE: Hardy-Weinberg equilibrium. For related individuals, this statistic can be inflated. For case-control study, we calculate HWE pvalues for all samples, case samples and controls samples separated by a colon.
    N_REF/N_HET/N_ALT: Number of samples carrying homozygous reference/heterozygous/homozygous alternative alleles. For case-control study, we calculate these three statistics for all samples, case samples and controls samples separated by a colon.
    U_STAT, SQRT_V_STAT: U and V statistics are score statistics and their covariance matrix. Details can be found in Dajiang Liu (2014) Nature Genetics.
    ALT_EFFSIZE: for continuous outcome, this is the estimated effect size; for binary outcome, this is the estimated log odds-ratio. We apply a new correction method when binary trait associations for related individuals are analyzed in standard linear mixed models. The log odds ratio is approximately correct for related individual analysis as well.


- Epacts   `out/test.single.b.score.epacts.gz` 
https://genome.sph.umich.edu/wiki/EPACTS#Output_Files
#CHROM	BEGIN	END	MARKER_ID	NS	AC	CALLRATE	MAF	PVALUE	SCORE	N.CASE	N.CTRL	AF.CASE	AF.CTRL
20	68303	68303	20:68303_A/G_Upstream:DEFB125	266	1	1	0.0018797	NA	NA	NA	NA	NA	NA
20	68319	68319	20:68319_C/A_Upstream:DEFB125	266	1.4467e-36	1	0	NA	NA	NA	NA	NA	NA
20	68396	68396	20:68396_C/T_Nonsynonymous:DEFB125	266	1	1	0.0018797	NA	NA	NA	NA	NA	NA
20	76635	76635	20:76635_A/T_Intron:DEFB125	266	1.534e-37	1	0	NA	NA	NA	NA	NA	NA
20	76689	76689	20:76689_T/C_Synonymous:DEFB125	266	0	1	0	NA	NA	NA	NA	NA	NA
20	76690	76690	20:76690_T/C_Nonsynonymous:DEFB125	266	1	1	0.0018797	NA	NA	NA	NA	NA	NA
20	76700	76700	20:76700_G/A_Nonsynonymous:DEFB125	266	0	1	0	NA	NA	NA	NA	NA	NA
20	76726	76726	20:76726_C/G_Nonsynonymous:DEFB125	266	0	1	0	NA	NA	NA	NA	NA	NA
20	76771	76771	20:76771_C/T_Nonsynonymous:DEFB125	266	3	1	0.0056391	0.68484	0.40587	145	121	0.013793	0.0082645

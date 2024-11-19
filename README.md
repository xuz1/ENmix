# ENmix: Quality control and analysis tools for Illumina DNA methylation BeadChip

# Introduction

The `r Biocpkg("ENmix")` package provides a set of quality control,
preprocessing/correction and data analysis tools for Illumina Methylation Beadchips.
It includes functions to read in raw idat data,
background correction, dye bias correction, probe-type bias adjustment,
 along with a number of additional tools.
These functions can be used to remove unwanted experimental noise and thus to
improve accuracy and reproducibility of methylation measures.
`r Biocpkg("ENmix")` functions
are flexible and transparent. Users have option to choose a single pipeline
command to finish all data pre-processing steps (including quality control,
background correction,
dye-bias adjustment, between-array normalization and probe-type bias correction) or
to use individual functions sequentially to perform data pre-processing in a more
customized manner. In addition the `r Biocpkg("ENmix")` package has selectable
complementary functions for efficient data visualization (such as QC plots, data
distribution plot, manhattan plot and Q-Q plot), quality control (identifing and filtering
low quality data points, samples, probes, and outliers, along with
imputation of missing values), identification of probes with multimodal
distributions due to SNPs or other factors, exploration of data variance
 structure using principal component regression analysis plot, preparation
of experimental factors related surrogate control variables
to be adjusted in downstream
statistical analysis, an efficient algorithm oxBS-MLE to estimate
5-methylcytosine and 5-hydroxymethylcytosine level; estimation of celltype
proporitons; methlation age calculation and differentially methylated
region (DMR) analysis.

Most `r Biocpkg("ENmix")` package can also support the data structure used by
several other related R packages, such as `r Biocpkg("minfi")`,
`r Biocpkg("wateRmelon")` and `r Biocpkg("ChAMP")`,
providing straightforward integration of
ENmix-corrected datasets for subsequent data analysis.

`r Biocpkg("ENmix")` readidat function does not depend on array annotation R packages.
It can directly read in Illuminal manifest file, which makes it easier to work with
newer array, such as MethylationEPICv2.0 and  mouse Beadchip.

The software is designed to support large scale data analysis, and provides
 multi-processor parallel computing options for most functions.


# List of functions

<b>Data acquisition</b>

<ul>
<li>`readidat()`:  Read idat files into R</li>
<li>`readmanifest()`:  Read array manifest file into R</li>
</ul>

<b>Quality control</b>

<ul>
<li>`QCinfo()`:      Extract and visualize QC information</li>
<li>`plotCtrl()`:    Generate internal control plots</li>
<li>`getCGinfo()`:   Extract CpG probe annotation information</li>
<li>`calcdetP()`:    Compute detection P values</li>
<li>`qcfilter()`:  Remove low quality values, samples or CpGs; remove outlier samples and perform imputation</li>
<li>`nmode()`:   Identify "gap" probes, i.e. those with multimodal distribution from underlying caused by underlying SNPs</li>
<li>`dupicc()`:     Calculate Introclass correlation coefficient (ICC) using data for duplicates</li>
<li>`freqpoly()`:  Frequency polygon plot for single variable</li>
<li>`multifreqpoly()`:  Frequency polygon plot for multiple variables</li>
</ul>

<b>Preprocessing</b>

<ul>
<li>`mpreprocess()`:      Preprocessing pipeline</li>
<li>`preprocessENmix()`:  ENmix background correction and dye bias correction</li>
<li>`relic()`:            RELIC dye bias correction</li>
<li>`norm.quantile()`:    Quantile normalization</li>
<li>`rcp()`:              RCP probe design type bias correction</li>
</ul>

<b>Differential methylated region (DMR) analysis</b>

<ul>
<li>`ipdmr()`:    ipDMR differentially methylated region analysis</li>
<li>`combp()`:    Combp differentially methylated region analysis</li>
</ul>

<b>Other functions</b>

<ul>
<li>`oxBS.MLE()`:  MLE estimates of 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC)</li>
<li>`estimateCellProp()`:    Estimate white blood cell type proportions</li>
<li>`methyAge()`:     Calculate methylation age</li>
<li>`methscore()`:    calculate various methylation predictors, including DNA methylation age, exposures and plasma protein levels.</li>
<li>`predSex()`:      Estimate sample sex</li>
<li>`ctrlsva()`:      Derive surrogate variables to control for experimental confounding using non-negative internal control probes</li>
<li>`pcrplot()`:      Principal component regression plot </li>
<li>`mhtplot()`:      P value manhattan plot </li>
<li>`p.qqplot()`:    P value Q-Q plot  </li>
<li>`B2M()`:          Convert Beta value to M value </li>
<li>`M2B()`:          Convert M value to Beta value  </li>
</ul>

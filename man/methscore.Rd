\name{methscore}
\alias{methscore}
\title{
DNA Methylation predictors
}
\description{
To calculate various methylation predictors, including DNA methylation age, plasma protein levels, exposures etc.
}
\usage{
methscore(datMeth,datPheno=NULL,fastImputation=FALSE,normalize=TRUE,GrimAgeComponent=NULL,UserRef=NULL,ForceUserRef=FALSE)
}
\arguments{
  \item{datMeth}{Methylation beta value matrix with CpG names(row names) and sample ids(column names).}
  \item{datPheno}{Phenotype data, must included columns SampleID, Age (year) and Female (0:male,1:female)}
  \item{fastImputation}{If "TRUE" mean methylation values will used
      for imputation, if "FALSE", KNN nearest neighbor method will be used.}
  \item{normalize}{TRUE or FALSE, if TRUE, user data will be normalized to a reference data using a modified RCP method}
  \item{GrimAgeComponent}{A data frame of grimage component from methylation age online calculator 
	  (https://dnamage.clockfoundation.org/user/login). It must include the following variables:
	 ("SampleID","DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS", "DNAmPAI1", "DNAmTIMP1"). 
         If the file is not provided, bAge will be calculated using PC Grimage components}
 \item{UserRef}{User provided methylation reference data will be used for some predictors (see details), must include variables cg and meth_mean}
 \item{ForceUserRef}{If TRUE, UserRef will be forced to use in normalization and imputation for all estimates}
}
\details{
	The original publications (see references) provided reference methylation mean values for all methylation predictors except for the following predictors: 
	PEDBE,EpiToc,EpiToc2,Zhang10CpG,Horvath2,MiAge,DNAmTL,PEDBE,GACPC,GARPC,GARRPC,Bohlin,Knight. Reference methylatin derived from Sister Study data were used 
	for these predictors by default. If UserRef was provided, it will be used to replace the default reference data. If set ForceUserRef=TRUE, UserRef will be 
	used for all predictors in normalization and CpG imputation. 
}
\value{
  A data frame with rows for samples and columns for methylation predictors. 
  Output file "summary_CpG_count.csv" has information about CpG counts and reference for each predictor.
}

\references{
Belsky et al. \emph{DunedinPACE, a DNA methylation biomarker of the pace of aging}. Elife. 2022 Jan 14;11:e73420. doi: 10.7554/eLife.73420

Bernabeu et al. \emph{Refining epigenetic prediction of chronological and biological age}. Genome Med. 2023 Feb 28;15(1):12. doi: 10.1186/s13073-023-01161-y

Bohlin et al. \emph{Prediction of gestational age based on genome-wide differentially methylated regions}.Genome Biol. 2016 Oct 7;17(1):207. doi: 10.1186/s13059-016-1063-4.

Hannum et al. \emph{Genome-wide methylation profiles reveal quantitative views of human aging rates}.Mol Cell. 2013 Jan 24;49(2):359-367. doi: 10.1016/j.molcel.2012.10.016. Epub 2012 Nov 21.

Higgins-Chen et al. \emph{A computational solution for bolstering reliability of epigenetic clocks: Implications for clinical trials and longitudinal tracking}.	Nat Aging. 2022 Jul; 2(7): 644–661.

Horvath. \emph{DNA methylation age of human tissues and cell types}. Genome Biol. 2013;14(10):R115. doi: 10.1186/gb-2013-14-10-r115.

Horvath et al. \emph{Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies}.Aging (Albany NY). 2018 Jul 26;10(7):1758-1775. doi: 10.18632/aging.101508.

Knight et al. \emph{An epigenetic clock for gestational age at birth based on blood methylation data}. Genome Biol. 2016 Oct 7;17(1):206. doi: 10.1186/s13059-016-1068-z.

Lee et al. \emph{Placental epigenetic clocks: estimating gestational age using placental DNA methylation levels}.Aging (Albany NY). 2019 Jun 30; 11(12): 4238–4253.

Levine et al. \emph{An epigenetic biomarker of aging for lifespan and healthspan}. Aging (Albany NY). 2018 Apr; 10(4): 573–591.

Li et al. \emph{Derivation and validation of an epigenetic frailty risk score in population-based cohorts of older adults}. Nat Commun. 2022 Sep 7;13(1):5269. doi: 10.1038/s41467-022-32893-x.

Lu et al. \emph{DNA methylation-based estimator of telomere length}.Aging (Albany NY). 2019 Aug 31; 11(16): 5895–5923. 

McCartney et al. \emph{Epigenetic prediction of complex traits and death}.Genome Biol. 2018 Sep 27;19(1):136. doi: 10.1186/s13059-018-1514-1.

McEwen et al. \emph{The PedBE clock accurately estimates DNA methylation age in pediatric buccal cells}.Proc Natl Acad Sci U S A. 2020 Sep 22; 117(38): 23329–23335.

Teschendorff et al. \emph{A comparison of epigenetic mitotic-like clocks for cancer risk prediction}.Genome Med. 2020 Jun 24;12(1):56. doi: 10.1186/s13073-020-00752-3.

Yang et al. \emph{Correlation of an epigenetic mitotic clock with cancer risk}. Genome Biol. 2016 Oct 3;17(1):205. doi: 10.1186/s13059-016-1064-3.

Youn et al. \emph{The MiAge Calculator: a DNA methylation-based mitotic age calculator of human tissue types}.Epigenetics. 2018; 13(2): 192–206.

Zhang et al. \emph{Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing}.Genome Med. 2019 Aug 23;11(1):54. doi: 10.1186/s13073-019-0667-1.

Zhang et al. \emph{DNA methylation signatures in peripheral blood strongly predict all-cause mortality}.Nat Commun. 2017 Mar 17;8:14617. doi: 10.1038/ncomms14617.
}

\author{
Zongli Xu
}

\examples{
\donttest{
require(minfiData)
path <- file.path(find.package("minfiData"),"extdata")
#based on rgDataset
rgSet <- readidat(path = path,recursive = TRUE)
meth=getmeth(rgSet)
beta=getB(meth)
pheno=data.frame(SampleID=colnames(beta),Age=c(23,45,52,36,58,43),Female=c(1,0,1,1,0,1))
mage=methscore(datMeth=beta,datPheno=pheno)
}
}

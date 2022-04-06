\name{hccframe}
\alias{hccframe}
\docType{data}
\title{
Liver Cancer Methylation Data
}
\description{
These data are a subset of subjects and CpG sites reported in the original paper where liver samples were assayed using the Illumina GoldenGate Methylation BeadArray Cancer Panel I. Technical replicate samples
were removed to ensure all samples were independent. The matched cirrhotic samples from subjects with hepatocellular carcinoma (HCC, labeled Tumor) were also excluded. Therefore methylation levels in liver tissue are provided for independent subjects whose liver was Normal (N=20),
cirrhotic but not having HCC (N=16, Cirrhosis non-HCC), and HCC (N=20, Tumor). 
}
\usage{data(hccframe)}
\format{
  A data frame with 56 observations on the following 46 variables.
  \describe{
    \item{\code{group}}{an ordered factor with levels \code{Normal} < \code{Cirrhosis non-HCC} < \code{Tumor}}
    \item{\code{CDKN2B_seq_50_S294_F}}{a numeric vector representing a  CpG site proportion methylation for CDKN2B}
    \item{\code{DDIT3_P1313_R}}{a numeric vector representing a  CpG site proportion methylation for DDIT3}
    \item{\code{ERN1_P809_R}}{a numeric vector representing a  CpG site proportion methylation for ERN1}
    \item{\code{GML_E144_F}}{a numeric vector representing a  CpG site proportion methylation for GML}
    \item{\code{HDAC9_P137_R}}{a numeric vector representing a  CpG site  proportion methylation for HDAC9}
    \item{\code{HLA.DPA1_P205_R}}{a numeric vector representing a  CpG site  proportion methylation for HLA.DPA1}
    \item{\code{HOXB2_P488_R}}{a numeric vector representing a  CpG site  proportion methylation for HOXB2}
    \item{\code{IL16_P226_F}}{a numeric vector representing a  CpG site  proportion methylation for IL16}
    \item{\code{IL16_P93_R}}{a numeric vector representing a  CpG site  proportion methylation for IL16}
    \item{\code{IL8_P83_F}}{a numeric vector representing a  CpG site  proportion methylation for IL8}
    \item{\code{MPO_E302_R}}{a numeric vector representing a  CpG site  proportion methylation for MPO}
    \item{\code{MPO_P883_R}}{a numeric vector representing a  CpG site  proportion methylation for MPO}
    \item{\code{PADI4_P1158_R}}{a numeric vector representing a  CpG site  proportion methylation for PADI4}
    \item{\code{SOX17_P287_R}}{a numeric vector representing a  CpG site  proportion methylation for SOX17}
    \item{\code{TJP2_P518_F}}{a numeric vector representing a  CpG site  proportion methylation for TJP2}
    \item{\code{WRN_E57_F}}{a numeric vector representing a  CpG site  proportion methylation for WRN}
    \item{\code{CRIP1_P874_R}}{a numeric vector representing a  CpG site  proportion methylation for CRIP1}
    \item{\code{SLC22A3_P634_F}}{a numeric vector representing a  CpG site  proportion methylation for SLC22A3}
    \item{\code{CCNA1_P216_F}}{a numeric vector representing a  CpG site  proportion methylation for CCNA1}
    \item{\code{SEPT9_P374_F}}{a numeric vector representing a  CpG site  proportion methylation for SEPT9}
    \item{\code{ITGA2_E120_F}}{a numeric vector representing a  CpG site  proportion methylation for ITGA2}
    \item{\code{ITGA6_P718_R}}{a numeric vector representing a  CpG site  proportion methylation for ITGA6}
    \item{\code{HGF_P1293_R}}{a numeric vector representing a  CpG site  proportion methylation for HGF}
    \item{\code{DLG3_E340_F}}{a numeric vector representing a  CpG site  proportion methylation for DLG3}
    \item{\code{APP_E8_F}}{a numeric vector representing a  CpG site  proportion methylation for APP}
    \item{\code{SFTPB_P689_R}}{a numeric vector representing a  CpG site  proportion methylation for SFTPB}
    \item{\code{PENK_P447_R}}{a numeric vector representing a  CpG site  proportion methylation for PENK}
    \item{\code{COMT_E401_F}}{a numeric vector representing a  CpG site  proportion methylation for COMT}
    \item{\code{NOTCH1_E452_R}}{a numeric vector representing a  CpG site  proportion methylation for NOTCH1}
    \item{\code{EPHA8_P456_R}}{a numeric vector representing a  CpG site  proportion methylation for EPHA8}
    \item{\code{WT1_P853_F}}{a numeric vector representing a  CpG site  proportion methylation for WT1}
    \item{\code{KLK10_P268_R}}{a numeric vector representing a  CpG site  proportion methylation for KLK10}
    \item{\code{PCDH1_P264_F}}{a numeric vector representing a  CpG site  proportion methylation for PCDH1}
    \item{\code{TDGF1_P428_R}}{a numeric vector representing a  CpG site  proportion methylation for TDGF1}
    \item{\code{EFNB3_P442_R}}{a numeric vector representing a  CpG site  proportion methylation for EFNB3}
    \item{\code{MMP19_P306_F}}{a numeric vector representing a  CpG site  proportion methylation for MMP19}
    \item{\code{FGFR2_P460_R}}{a numeric vector representing a  CpG site  proportion methylation for FGFR2}
    \item{\code{RAF1_P330_F}}{a numeric vector representing a  CpG site  proportion methylation for RAF1}
    \item{\code{BMPR2_E435_F}}{a numeric vector representing a  CpG site  proportion methylation for BMPR2}
    \item{\code{GRB10_P496_R}}{a numeric vector representing a  CpG site  proportion methylation for GRB10}
    \item{\code{CTSH_P238_F}}{a numeric vector representing a  CpG site  proportion methylation for CTSH}
    \item{\code{SLC6A8_seq_28_S227_F}}{a numeric vector representing a  CpG site proportion methylation for SLC6A8}
    \item{\code{PLXDC1_P236_F}}{a numeric vector representing a  CpG site proportion methylation for PLXDC1}
    \item{\code{TFE3_P421_F}}{a numeric vector representing a  CpG site proportion methylation for TFE3}
    \item{\code{TSG101_P139_R}}{a numeric vector representing a  CpG site proportion methylation for TSG101}
  }
}
\source{
The full dataset is available as GSE18081 from Gene Expression Omnibus at 
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18081}
}
\references{
Archer KJ, Mas VR, Maluf DG, Fisher RA. High-throughput assessment of CpG site methylation for distinguishing between HCV-cirrhosis and HCV-associated hepatocellular carcinoma. Molecular Genetics and Genomics, 283(4): 341-349, 2010.
}
\seealso{
See Also as \code{ \link{ordinalgmifs} }
}
\examples{
data(hccframe)
}
\keyword{datasets}

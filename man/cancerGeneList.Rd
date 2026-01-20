\name{cancerGeneList}
\alias{cancerGeneList}
\docType{data}
\title{
Comprehensive List of Cancer Driver Genes
}
\description{
A dataset containing a comprehensive list of cancer driver genes, including Oncogenes and Tumor Suppressor Genes (TSGs), annotated from various resources such as OncoKB, COSMIC, and MSK-IMPACT.
}
\usage{data("cancerGeneList")}
\format{
  A data frame with observations on the following variables:
  \describe{
    \item{\code{Hugo Symbol}}{HGNC gene symbol}
    \item{\code{Entrez Gene ID}}{Entrez gene ID}
    \item{\code{GRCh37 Isoform}}{Ensembl transcript ID for GRCh37}
    \item{\code{GRCh37 RefSeq}}{RefSeq transcript ID for GRCh37}
    \item{\code{GRCh38 Isoform}}{Ensembl transcript ID for GRCh38}
    \item{\code{GRCh38 RefSeq}}{RefSeq transcript ID for GRCh38}
    \item{\code{Gene Type}}{Classification of the gene (e.g., ONCOGENE, TSG, ONCOGENE_AND_TSG)}
    \item{\code{# of occurrence within resources (Column J-P)}}{Number of resources listing this gene}
    \item{\code{OncoKB Annotated}}{Whether the gene is annotated in OncoKB}
    \item{\code{MSK-IMPACT}}{Whether the gene is in MSK-IMPACT panel}
    \item{\code{MSK-HEME}}{Whether the gene is in MSK-HEME panel}
    \item{\code{FOUNDATION ONE}}{Whether the gene is in Foundation One panel}
    \item{\code{FOUNDATION ONE HEME}}{Whether the gene is in Foundation One Heme panel}
    \item{\code{Vogelstein}}{Whether the gene is listed in Vogelstein et al.}
    \item{\code{COSMIC CGC (v99)}}{Whether the gene is in the COSMIC Cancer Gene Census v99}
    \item{\code{Gene Aliases}}{Alternative gene names/aliases}
  }
}
\source{
This dataset was downloaded from the OncoKB Cancer Gene List 
(\url{https://www.oncokb.org/cancer-genes}).

OncoKB is a precision oncology knowledge base developed and maintained by 
Memorial Sloan Kettering Cancer Center (MSKCC). The knowledge base contains 
information about the effects of specific cancer gene alterations on 
treatment response and prognosis.

The cancer gene list integrates data from multiple sources including:
\itemize{
  \item OncoKB (\url{https://www.oncokb.org/})
  \item COSMIC Cancer Gene Census (\url{https://cancer.sanger.ac.uk/census})
  \item MSK-IMPACT and MSK-HEME panels
  \item Foundation Medicine panels
  \item Vogelstein et al. cancer gene classification
}
}
\references{
Chakravarty D, Gao J, Phillips SM, et al. OncoKB: A Precision Oncology 
Knowledge Base. JCO Precision Oncology 2017;1:1-16. 
doi:10.1200/PO.17.00011

Sondka Z, Bamford S, Cole CG, et al. The COSMIC Cancer Gene Census: 
describing genetic dysfunction across all human cancers. 
Nature Reviews Cancer 2018;18:696-705. doi:10.1038/s41568-018-0060-1

Vogelstein B, Papadopoulos N, Velculescu VE, et al. Cancer genome 
landscapes. Science 2013;339:1546-1558. doi:10.1126/science.1235122
}
\keyword{datasets}

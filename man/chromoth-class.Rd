\name{chromoth-class}
\docType{class}

% Class:
\alias{class:chromoth}
\alias{chromoth-class}
\alias{chromoth}
\alias{getSVs}
\alias{getSegs}



\title{chromoth objects}

\description{
	The chromoth object is an object returned by the function \link{detect_chromothripsis}.
	}

\details{
	This object stores the results returned by the function \link{detect_chromothripsis}. 
	This object may be coerced to a data.frame using \code{as(object,"data.frame")}. 
		The obtained data.frame object contains 5 columns, "chrom","clusterSize","SVpvalue","CNVpvalue", and "score". The "chrom" column is the chromosome name. The "clusterSize" column is the maximum cluster size that the SVs can be clustered on a chromosome. The "SVpvalue" column is the pvalue given by performing statistical tests based on structural variations. The "CNVpvalue" is the pvalue given by performing statistical tests based on copy number variaitons. The "score" equals to \code{-log10(SVpvalue)-log10(CNVpvalue)}.
}

\section{Accessors}{
	In the following, x is a chromoth object.
	  \describe{
	    \item{getSVs(x)}{
	        Get the structrual variations; Note that since \link{detect_chromothripsis} first try to filter possible false structrual variations, the structrual variations returned by this function might be different from the structrual variations provided to \link{detect_chromothripsis}.
    		}
	  \item{getSegs(x)}{
		Get the copy number segmentation.
		}
	}
}

\seealso{
	\link{plot_chromothripsis}, \link{detect_chromothripsis}
	}

\author{Ruibin Xi}

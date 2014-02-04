\name{get_elbow_limma}
\alias{get_elbow_limma}
\title{Calculates fold values from an MArrayLM object.}
\usage{
  get_elbow_limma(marraylm, columns = NULL)
}
\arguments{
  \item{marraylm}{is the MArrayLM object to analyze.}

  \item{columns}{is the list of sample columns to obtain
  the elbow fold cut-off values for.  This can be specified
  as a vector or a single value.}
}
\value{
  a matrix specified as follows \itemize{ \item{columns ---
  (1) \dQuote{up_limit}, the upper ELBOW fold-change
  cut-off value; (2) \dQuote{low_limit}, the lower ELBOW
  fold-change cut-off value} \item{rows --- one row per
  sample, specified by the parameter \dQuote{columns.}} }
}
\description{
  Performs the ELBOW fold change test on an MArrayLM
  object.  This is a wrapper class to help integrate the
  ELBOW method into Bioconductor. followed tutorial from:
  http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
}
\details{
  Some of the code in this method is based on how toptable
  accesses the MArrayLM object to read fold values.
  Therefore, any MArrayLM object which works with the
  toptable method should also work with this method.
}
\examples{
########
		# LOAD DATA INTO LIMMA
		########
		library("limma")

		# load a filtered expression set into R
		# NOTE: see the vignette for instructions on preparing
		#       a filtered dataset with your own data.
		data(GSE20986_eset_exprs, package="ELBOW")
		data(GSE20986_design, package="ELBOW")

		# fit the linear model to the filtered expression set
		fit <- lmFit(GSE20986_eset_exprs, GSE20986_design)

		# set up a contrast matrix to compare tissues v cell line
		contrast.matrix <- makeContrasts(huvec_choroid = huvec - choroid, huvec_retina = huvec - retina, huvec_iris <- huvec - iris, levels=GSE20986_design)

		# Now the contrast matrix is combined with the per-probeset linear model fit.
		huvec_fits <- contrasts.fit(fit, contrast.matrix)
		huvec_ebFit <- eBayes(huvec_fits)

		########
		# GET THE ELBOW LIMIT (this function)
		########
		get_elbow_limma(huvec_ebFit)
}


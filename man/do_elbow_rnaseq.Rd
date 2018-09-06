\name{do_elbow_rnaseq}
\alias{do_elbow_rnaseq}
\title{Calculates fold values from an MArrayLM object.}
\usage{
  do_elbow_rnaseq(rnaSeq)
}
\arguments{
  \item{rnaSeq}{is the CountDataSet object to analyze.}
}
\value{
  a matrix specified as follows \itemize{ \item{columns ---
  (1) \dQuote{up_limit}, the upper ELBOW fold-change
  cut-off value; (2) \dQuote{low_limit}, the lower ELBOW
  fold-change cut-off value} \item{rows --- one row per
  sample, specified by the parameter \dQuote{columns.}} }
}
\description{
  Performs the ELBOW fold change test on an CountDataSet
  from DESeq object.  This is a wrapper class to help
  integrate the ELBOW method into Bioconductor. followed
  tutorial from:
  http://cgrlucb.wikispaces.com/Spring+2012+DESeq+Tutorial
}
\examples{
# install the DESeq libraries
		#if (!requireNamespace("BiocManager", quietly=TRUE))
    		#install.packages("BiocManager")
		#BiocManager::install("DESeq")

		## download the table
		library("DESeq")

		# the following bam file dataset was obtained from:
		# http://cgrlucb.wikispaces.com/file/view/yeast_sample_data.txt
		# it has been downloaded into this package for speed convenience.
		filename <- system.file("extdata", "yeast_sample_data.txt", package = "ELBOW")

		count_table <- read.table(filename, header=TRUE, sep="\\t", row.names=1)
		expt_design <- data.frame(row.names = colnames(count_table), condition = c("WE","WE","M","M","M"))
		conditions = expt_design$condition
		data <- newCountDataSet(count_table, conditions)
		data <- estimateSizeFactors(data)
		data <- as(data, "CountDataSet")
		## data <- estimateVarianceFunctions(data)
		data <- estimateDispersions(data)

		# this next step is essential, but it takes a long time...
		# so, just like a good cooking show we will skip this step
		# and load a finished version.
		#results <- nbinomTest(data, "M", "WE")

		# The below two code lines load a copy of the above dataset
		# which has already been processed by:
		#     results <- nbinomTest(data, "M", "WE")
		# For your own real data, you must use:
		#     results <- nbinomTest(data, "M", "WE")'
		# Instead of the two lines below:
		data(yeast_nbinomTest_results, package="ELBOW")
		results <- yeast_nbinomTest_results

		# obtain the elbow limit for the dataset
		# the final step in the analysis pipeline
		do_elbow_rnaseq(results)
}


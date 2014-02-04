\name{plot_dataset}
\alias{plot_dataset}
\title{plots data to see if it matches an elbow distribution.}
\usage{
  plot_dataset(my_data, column = NULL, upper_limit = NULL,
    lower_limit = NULL)
}
\arguments{
  \item{my_data}{is the object, or array, to plot.}

  \item{column}{the column in the object/array to plot, if
  applicable.}

  \item{upper_limit}{the upper ELBOW limit for the data
  (optional)}

  \item{lower_limit}{the lower ELBOW limit for the data
  (optional)}
}
\description{
  A debug function for plotting data.  This function plots
  an array or R object, and optionally draws the ELBOW
  limits on the plot. This function is mainly useful for
  the development and maintenance of the ELBOW method R
  package.
}
\examples{
# install the DESeq libraries
		#source("http://www.bioconductor.org/biocLite.R")
		#biocLite("DESeq")

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

		# plot the results w/o elbow (useful if do_elbow doesn't work)
		#plot_dataset(results, "log2FoldChange")

		# plot the results w/elbow
		limits <- do_elbow_rnaseq(results)
		plot_dataset(results, "log2FoldChange", limits$up_limit, limits$low_limit)
}


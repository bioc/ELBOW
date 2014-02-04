\name{yeast_nbinomTest_results}
\docType{data}
\alias{yeast_nbinomTest_results}
\title{Pre-processed Yeast sample RNA-seq data}
\description{
	Pre-processed Yeast sample RNA-seq data, obtained from the website:
	\url{http://cgrlucb.wikispaces.com/Spring+2012+DESeq+Tutorial}
	
	For more information, please contact the author of the webpage.
	
	The content from the tutorial's webpage is licensed as follows:
	Contributions to http://cgrlucb.wikispaces.com/ are licensed under a Creative Commons Attribution Share-Alike 3.0 License. Creative Commons Attribution Share-Alike 3.0 License.
}
\usage{yeast_nbinomTest_results}
\format{
	An 8 column table defined as such:
	\itemize{
		\item{id --- the Gene ID reprented by the row}
		\item{baseMean --- mean estimated from both conditions (with normalization by size factors)}
		\item{baseMeanA --- mean estimated from the first condition specified as an argument to nbinomTest (in this case, \dQuote{M})}
		\item{baseMeanB --- mean estimated from the second condition (in this case, \dQuote{WE})}
		\item{foldChange --- ratio of the estimated means}
		\item{log2FoldChange --- log2 of the ratio}
		\item{pval --- uncorrected p-value from the negative binomial test}
		\item{padj --- p-value adjusted for multiple testing using Benjamini-Hochberg to estimate the false discovery rate}
		\item{resVarA --- dispersion estimated from the first condition}
		\item{resVarB --- dispersion estimated from the second condition}
	}
}
\source{http://cgrlucb.wikispaces.com/Spring+2012+DESeq+Tutorial}

\references{Lee, H. \dQuote{Spring 2012 DESeq Tutorial.} \url{http://cgrlucb.wikispaces.com/Spring+2012+DESeq+Tutorial}}
\keyword{datasets}

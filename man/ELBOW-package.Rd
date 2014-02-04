\name{ELBOW}
\alias{pkg_elbowlib-package}
\alias{pkg_elbowlib}
\alias{elbowlib}
\docType{package}
\title{
	Evaluating foLd change By the lOgit Way
}
\description{
	The "elbow" method an improved fold change test for determining cut off for biologically significant changes in expression levels in transcriptomics.
}
\details{
	\tabular{ll}{
		Package: \tab ELBOW\cr
		Type: \tab Package\cr
		Version: \tab 1.0\cr
		Date: \tab 2013-08-08\cr
		License: \tab Creative Commons 3.0 Attribution
			+ ShareAlike\cr
			\tab (see: \url{http://creativecommons.org/licenses/by-sa/3.0/})\cr
	}
	Elbow an improved fold change test that uses cluster analysis and pattern recognition to set cut off limits that are derived directly from intrareplicate variance without assuming a normal distribution for as few as 2 biological replicates. Elbow also provides the same consistency as fold testing in cross platform analysis. Elbow has lower false positive and false negative rates than standard fold testing when both are evaluated using T testing and Statistical Analysis of Microarray using 12 replicates (six replicates each for initial and final conditions). Elbow provides a null value based on initial condition replicates and gives error bounds for results to allow better evaluation of significance.

	\bold{Abstract Reference:} \cr
Conference Proceeding: Zhang, X., Bjorklund, N. K., Rydzak, T., Sparling, R., Alvare, G., Fristensky, B. (April 2013) \dQuote{The Elbow Method for deciding significant fold change cutoffs of differentially expressed genes.} \emph{Recomb 2013 17th International Conference on Research in Computational Biology}.

	\bold{Paper Reference:} \cr
Zhang, X., Bjorklund, N. K., Alvare, G., Rydzak, T., Sparling, R., Fristensky, B. (2013) \dQuote{Elbow, an improved fold test method for transcriptomics.} Departments of Plant Science and Microbiology, University of Manitoba, Winnipeg, Canada, R3T 2N2

The corresponding author: Brian Fristensky \email{frist@cc.umanitoba.ca}
}
\author{
	Xiangli Zhang, Natalie Bjorklund, Graham Alvare, Tom Ryzdak,
		Richard Sparling, Brian Fristensky

	Maintainers: Graham Alvare \email{alvare@cc.umanitoba.ca},
		Xiangli Zhang \email{zhangju@cc.umanitoba.ca}
}

% ~ put references to the literature/web site here ~
% ~ use the following guide for citation formatting: http://developer.r-project.org/Rds.html ~
\references{
	Claeskens, G. and Hjort N. L. (2008) \emph{Model Selection and Model Averaging.} Cambridge, England: Cambridge University Press.
	
	Cui X. and Churchill G. A. (2003) \dQuote{Statistical tests for differential expression in cDNA microarray experiments.} \emph{Genome Biol}, Vol. \bold{4}, p. 210.
	
	Dalman, M. R., Deeter, A., Nimishakavi, G. and Duan, Z. H. (2012) \dQuote{Fold change and p-value cutoffs significantly alter microarray interpretations.} \emph{BMC Bioinformatics}, Vol. \bold{13} (Suppl. 2), p. S11
	
	Faraway, J. J. (2006) \emph{Extending the Linear Model with R: Generalized Linear, Mixed Effects and Nonparametric Regression Models.} Boca Raton, Florida: Chapman and Hall/CRC.
	
	Guo, L., Lobenhofer, E. K., et al. (2006) \dQuote{Rat toxicogenomic study reveals analytical consistency across microarray platforms.} \emph{Nature Biotechnol}, Vol. \bold{24}, pp. 1162--1169.
	
	Klebanov, L., Qiu, X., Welle, S., Yakovlev, A. (2007) \dQuote{Statistical methods and microarray data.} \emph{Nature Biotechnol}, Vol. \bold{25}, p.1.
	
	MAQC Consortium (2006) \dQuote{The MicroArray Quality Control (MAQC) project shows inter- and intraplatform reproducibility of gene expression measurements.} \emph{Nature Biotechnol}, September, Vol. \bold{24}, No. 9, pp. 1151--1161.
	
	Minty, J. J., Lesnefsky, A. A., et al. (2011) \dQuote{Evolution combined with genomic study elucidates genetic bases of isobutanol tolerance in Escherichia coli.} \emph{Microb Cell Fact.} Vol. \bold{10}, p. 18.

	\emph{NIST/SEMATECH e-Handbook of Statistical Methods}, April 2012, \href{http://www.itl.nist.gov/div898/handbook/}{http://www.itl.nist.gov/div898/handbook/}.
	
	Oshlack, A., Robinson, M. D., Young, M. D., (2010) \dQuote{From RNA-seq reads to differential expression results.} \emph{Genome Biol}, Vol. \bold{11}, p. 220.

	Shi, L., Tong, W., et al. (2005) \dQuote{Cross-platform comparability of microarray technology: intra-platform consistency and appropriate data analysis procedures are essential.} \emph{BMC Bioinformatics}, Vol. \bold{6} (Suppl. 2), p. S12.
	
	Sjogren, A., Kristiansson, E., Rudemo, M., Nerman, O. (2007) \dQuote{Weighted analysis of general microarray experiments.} \emph{BMC Bioinformatics}, Vol. \bold{8}, p. 387.
	
	Tan, P. K., Downey, T. J., et al. (2003) \dQuote{Evaluation of gene expression measurements from commercial microarray platforms.} \emph{Nucleic Acids Res}, Vol. \bold{31}, No.19, pp. 5676--5684.
	
	Thorndike, R. L. (December 1953) \dQuote{Who Belong in the Family?} Psychometrika, Vol. \bold{18}, No. 4, pp. 267--276.
	
	Tilstone, C. (2003) \dQuote{Vital statistics.} \emph{Nature}, Vol. \bold{424}, p. 611
	
	Tusher, V. G., Tibshirani, R., Chu, G. (2001) \dQuote{Significance analysis of microarrays applied to the ionizing radiation response.} \emph{Proc Natl Acad Sci}, Aug. 28, Vol. \bold{98} No. 9: pp. 5116--5121.
}

% ~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
% ~~ the R documentation directory ~~
\keyword{ package }
\keyword{ ELBOW }

% ~~ Optional links to other man pages, e.g. ~~
% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
\seealso{
		See \link[ELBOW:analyze_elbow]{analyze_elbow} for doing 
			a full ELBOW analysis and plot. \cr
		See \link[ELBOW:do_elbow]{do_elbow} if you want to extract
			only the ELBOW cut-off values.
}

% ~~ simple examples of the most important functions ~~
\examples{
	# read in the EcoliMutMA sample data from the package
	data(EcoliMutMA, package="ELBOW")
	csv_data <- EcoliMutMA
	# - OR - Read in a CSV file (uncomment - remove the #'s
	#        - from the line below and replace 'filename' with
	#        the CSV file's filename)
	# csv_data <- read.csv(filename)
	
	# set the number of initial and final condition replicates both to three
	init_count  <- 3
	final_count <- 3
	
	# Parse the probes, intial conditions and final conditions
	# out of the CSV file.  Please see: extract_working_sets
	# for more information.
	#
	# init_count should be the number of columns associated with
	#       the initial conditions of the experiment.
	# final_count should be the number of columns associated with
	#       the final conditions of the experiment.
	working_sets <- extract_working_sets(csv_data, init_count, final_count)
	
	probes <- working_sets[[1]]
	initial_conditions <- working_sets[[2]]
	final_conditions <- working_sets[[3]]

	# Uncomment to output the plot to a PNG file (optional)
	# png(file="output_plot.png")

	# Analyze the elbow curve.
	sig <- analyze_elbow(probes, initial_conditions, final_conditions)

	# write the significant probes to 'signprobes.csv'
	write.table(sig,file="signprobes.csv",sep=",",row.names=FALSE)
}

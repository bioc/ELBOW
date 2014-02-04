\name{analyze_elbow}
\alias{analyze_elbow}
\title{extracts all elbow statistics and plots and elbow curve.}
\usage{
  analyze_elbow(probes, initial_conditions,
    final_conditions, gtitle = "")
}
\arguments{
  \item{probes}{the data set of probes (as a data.frame).}

  \item{initial_conditions}{a 2D data.fram containing all
  of the replicate gene expression values corresponding to
  the experiment's initial conditions.}

  \item{final_conditions}{a 2D data.fram containing all of
  the replicate gene expression values corresponding to the
  experiment's final conditions.}

  \item{gtitle}{the title to display for the graph.}
}
\value{
  a list of all significant probes
}
\description{
  Analyzes: \itemize{ \item{the Elbow cut-offs} \item{the
  Elbow curve variance} \item{the upper and lower error
  Elbow curves} \item{\eqn{\log\chi^2} p-value for the
  Elbow curve model} }
}
\details{
  Then plots the data to a curve - AND - prints the
  statistics to both the terminal and the plot canvas.
}
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


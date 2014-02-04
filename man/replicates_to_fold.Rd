\name{replicates_to_fold}
\alias{replicates_to_fold}
\title{calculates fold values from a table of experimental data}
\usage{
  replicates_to_fold(probes, initial_conditions,
    final_conditions)
}
\arguments{
  \item{probes}{the data set of probes.}

  \item{initial_conditions}{a data set of replicates
  corresponding to initial conditions.}

  \item{final_conditions}{a data set of replicates
  corresponding to final conditions.}
}
\value{
  a data.frame comprised of probe names and fold-change
  values
}
\description{
  Converts a data.frame of probes with initial and final
  experimental conditions to a table with probes and fold
  change values.
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

		my_data <- replicates_to_fold(probes, initial_conditions, final_conditions)
}


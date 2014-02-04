\name{elbow_variance}
\alias{elbow_variance}
\title{calculates the variance for the upper and lower limits of an Elbow curve}
\usage{
  elbow_variance(probes, initial_conditions,
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
  a list containing the following keys: \itemize{
  \item{max_upper --- the maximum upper elbow limit
  \emph{(most positive)}} \item{min_upper --- the minimum
  upper elbow limit} \item{max_upper --- the maximum lower
  elbow limit} \item{min_upper --- the minimum lower elbow
  limit \emph{(most negative)}} }
}
\description{
  The following function determines the variation of cut
  off limits for significance between individual subsets.
  The first step is to generate the subsets, then the
  function determines the highest and lowest values for the
  upper limit and the highest and lowest values for the
  lower limit. These values are used to determine and upper
  and lower error value reflecting fold dataset variance.
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

		# compute the elbow for the dataset
		limits <- do_elbow(data.frame(my_data$fold))

		plus_minus <- elbow_variance(probes, initial_conditions, final_conditions)
}


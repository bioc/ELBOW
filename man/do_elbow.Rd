\name{do_elbow}
\alias{do_elbow}
\title{calculates the upper and lower elbow cut-off fold-values}
\usage{
  do_elbow(fold_values)
}
\arguments{
  \item{fold_values}{a data.frame containing all of the
  fold values to calculate the elbow for.}
}
\value{
  a list containing the following: \itemize{ \item{up_limit
  --- the upper (most positive) limit of the elbow curve}
  \item{low_limit --- the lower (most negative) limit of
  the elbow curve} }
}
\description{
  calculates the upper and lower elbow cut-off fold-values
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
}


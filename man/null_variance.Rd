\name{null_variance}
\alias{null_variance}
\title{determines the variance for a null Elbow curve}
\usage{
  null_variance(my_data, upper_limit, lower_limit,
    initial_conditions)
}
\arguments{
  \item{my_data}{A table to analyze the null variance of.
  The columns in the table should be as follows: \itemize{
  \item{probes --- one column containing the names of the
  probes} \item{initial conditions --- one or more columns
  containing the initial conditions of the experiment}
  \item{final conditions --- one or more columns containing
  the final conditions of the experiment} }}

  \item{upper_limit}{the actual upper limit cut-off for the
  Elbow curve}

  \item{lower_limit}{the actual lower limit cut-off for the
  Elbow curve}

  \item{initial_conditions}{a data set of replicates
  corresponding to initial conditions.}
}
\value{
  A list containing the following keys: \itemize{
  \item{prmax --- the error limit based on the maximum
  value for each probe} \item{prmed --- the null (median)
  value for each probe} \item{prmin --- the error limit
  based on the minimum value for each probe} }
}
\description{
  determine an upper and lower error limit using all
  possible subsets of differences between initial
  conditions and then use the median value for each probe
  to create the null set. Maximum and minimum values are
  used to set upper and lower error bounds.
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

		max_upper_variance <- plus_minus$max_upper
		min_upper_variance <- plus_minus$min_upper
		max_lower_variance <- plus_minus$max_lower
		min_lower_variance <- plus_minus$min_lower

		# rounded number for nice appearance graph
		upper_limit <- round(limits[[1]],digits=2)

		# rounded number nice appearance for graph
		lower_limit <- round(limits[[2]],digits=2)

		p_limits <- null_variance(my_data, upper_limit, lower_limit, initial_conditions)

		prowmin <- p_limits[[1]]
		prowmax <- p_limits[[2]]
		prowmedian <- p_limits[[3]]
}


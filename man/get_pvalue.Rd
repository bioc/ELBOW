\name{get_pvalue}
\alias{get_pvalue}
\title{The calculated \eqn{\log\chi^2} p-value for the fit
			of the elbow curve to the model.}
\usage{
  get_pvalue(my_data, upper_limit, lower_limit)
}
\arguments{
  \item{my_data}{A table (data.frame) to analyze the elbow
  variance of. The columns in the table should be as
  follows: \itemize{ \item{probes --- one column containing
  the names of the probes} \item{fold --- the fold values
  for the table} }}

  \item{upper_limit}{the upper limit/cut-off for the
  elbow.}

  \item{lower_limit}{the lower limit/cut-off for the
  elbow.}
}
\value{
  The calculated \eqn{\log\chi^2} p-value for the elbow
  curve.
}
\description{
  The following function determines a p value from the
  comparison of the slope of a regression line (i.e.
  logit).  The function calculates a fold value from a set
  of initial conditions (replicate set 1) substracted from
  a set of final conditions (replicate set 2).  The glm
  function is used to assess the fit of the dataset to the
  model, assuming family = binomial (link = "probit"),
  maxit = 1000. Determining the \eqn{\log\chi^2} p-value
  for the model, null equals no significance for the set of
  significant probes to the model, in which case the Elbow
  method cannot be used to assess the dataset.  If
  significant probes are important to the model, then the
  \eqn{\log\chi^2} p-value will be below 0.05, in which
  case the Elbow method can be used to assess the dataset.
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

		pvalue1 <- get_pvalue(my_data, upper_limit, lower_limit)
}


\name{plot_elbow}
\alias{plot_elbow}
\title{plots Elbow curve data.}
\usage{
  plot_elbow(my_data, upper_limit, lower_limit, pvalue1,
    prowmin, prowmax, prowmedian, max_upper_variance,
    min_upper_variance, max_lower_variance,
    min_lower_variance, gtitle = "")
}
\arguments{
  \item{my_data}{A table (data.frame) to plot.  The columns
  in the table should be as follows: \itemize{ \item{probes
  --- one column containing the names of the probes}
  \item{fold --- the fold values for the table} }}

  \item{upper_limit}{the upper limit/cut-off for the
  elbow.}

  \item{lower_limit}{the lower limit/cut-off for the
  elbow.}

  \item{pvalue1}{the \eqn{\log\chi^2} p-value for the elbow
  curve.}

  \item{prowmin}{the error limit based on the maximum value
  for each probe.}

  \item{prowmax}{the error limit based on the minimum value
  for each probe.}

  \item{prowmedian}{the null (median) value for each
  probe.}

  \item{max_upper_variance}{the maximum upper elbow limit
  \emph{(most positive)}.}

  \item{min_upper_variance}{the minimum upper elbow limit.}

  \item{max_lower_variance}{the maximum lower elbow limit.}

  \item{min_lower_variance}{the minimum lower elbow limit
  \emph{(most negative)}.}

  \item{gtitle}{the title to display for the graph.}
}
\description{
  Plots an elbow curve and its associated data: \itemize{
  \item{the upper and lower elbow limits for the curve}
  \item{the upper, lower, and median initial condition
  Elbow plots} \item{the \eqn{\log\chi^2} p-value for the
  Elbow curve} \item{the variance for the upper and lower
  Elbow cut-off values} }
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
		plot_elbow(my_data, upper_limit, lower_limit, pvalue1, prowmin, prowmax, prowmedian, max_upper_variance, min_upper_variance, max_lower_variance, min_lower_variance, "Title of the ELBOW Plot")
}


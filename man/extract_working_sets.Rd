\name{extract_working_sets}
\alias{extract_working_sets}
\title{extract the initial and final conditions from a table.}
\usage{
  extract_working_sets(first_my_data, replicate_count1,
    replicate_count2)
}
\arguments{
  \item{first_my_data}{A table to translate.  The columns
  in the table should be as follows: \itemize{ \item{probes
  --- one column containing the names of the probes}
  \item{initial conditions --- one or more columns
  containing the initial conditions of the experiment}
  \item{final conditions --- one or more columns containing
  the final conditions of the experiment} }}

  \item{replicate_count1}{the number of columns
  (replicates) corresponding to initial conditions.}

  \item{replicate_count2}{the number of columns
  (replicates) corresponding to final conditions.}
}
\value{
  a list containing the following: \itemize{ \item{probes
  --- the set of probes.} \item{initial --- a data set of
  replicates corresponding to initial conditions.}
  \item{final --- a data set of replicates corresponding to
  final conditions.} }
}
\description{
  Extract the initial and final conditions from a table.
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
}


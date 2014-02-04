#####
#TODO - FUTURE
#####
# How-to 'roxygenize' (from: http://stackoverflow.com/questions/12460347/processing-example-files-with-roxygen2-backslashes-are-duplicated-dontrun-bec):
#	require("roxygen2")
#	roxygenize(package.dir="/bigfiles/Elbow_Package/ELBOW", overwrite = TRUE, unlink.target = FALSE, roclets = c("collate", "namespace", "rd"))
#
# OR
#	Rscript --default-packages=roxygen2,stats,methods -e 'roxygenize(package.dir="~/Desktop/Elbow_Aug2013/pkg_elbowlib", overwrite = TRUE, unlink.target = FALSE, roclets = c("collate", "namespace", "rd")'
#
# How-to test roxygen documentation (via PDF):
# 	Generic:  rm /tmp/foo.pdf; R CMD Rd2pdf --pdf --title='Test of foo' -o /tmp/foo.pdf *.Rd
#	Specific: rm /tmp/foo.pdf; R CMD Rd2pdf --pdf --title='ELBOW' -o /tmp/foo.pdf ELBOW-package.Rd replicates_to_fold.Rd do_elbow.Rd elbow_variance.Rd null_variance.Rd get_pvalue.Rd plot_elbow.Rd extract_working_sets.Rd analyze_elbow.Rd EcoliMutMA.Rd EcoliWT.Rd
#
# Reference docs (manuals for roxygen, etc.)
# 	http://roxygen.org/useR/
# 	https://ebi-forecast.igb.illinois.edu/redmine/projects/pecan-1-2-5/wiki/Roxygen2
#
# Citation information
#	http://stat.ethz.ch/R-manual/R-devel/library/utils/html/citation.html
#	http://cran.r-project.org/manuals.html
#	http://images.webofknowledge.com/WOK46/help/WOS/N_abrvjt.html
#	http://developer.r-project.org/Rds.html
#
# How to download files:
#   download.file("http://cgrlucb.wikispaces.com/file/view/yeast_sample_data.txt", "yeast_sample_data.txt", "auto")

##########
#' Extract the initial and final conditions from a table.
#'
#' @export
#' @title extract the initial and final conditions from a table.
#' @param first_my_data A table to translate.
#'		The columns in the table should be as follows:
#'		\itemize{
#'			\item{probes --- one column containing the names of the probes}
#'			\item{initial conditions --- one or more columns containing the
#'				  initial conditions of the experiment}
#'			\item{final conditions --- one or more columns containing the
#'				  final conditions of the experiment}
#'		}
#' @param replicate_count1 the number of columns (replicates)
#'		corresponding to initial conditions.
#' @param replicate_count2 the number of columns (replicates)
#'		corresponding to final conditions.
#' @return a list containing the following:
#'		\itemize{
#'			\item{probes --- the set of probes.}
#'			\item{initial --- a data set of replicates
#'						corresponding to initial conditions.}
#'			\item{final --- a data set of replicates
#'						corresponding to final conditions.}
#'		}
#' @examples
#'		# read in the EcoliMutMA sample data from the package
#'		data(EcoliMutMA, package="ELBOW")
#'		csv_data <- EcoliMutMA
#'		# - OR - Read in a CSV file (uncomment - remove the #'s
#'		#        - from the line below and replace 'filename' with
#'		#        the CSV file's filename)
#'		# csv_data <- read.csv(filename)
#'		
#'		# set the number of initial and final condition replicates both to three
#'		init_count  <- 3
#'		final_count <- 3
#'	
#'		# Parse the probes, intial conditions and final conditions
#'		# out of the CSV file.  Please see: extract_working_sets
#'		# for more information.
#'		#
#'		# init_count should be the number of columns associated with
#'		#       the initial conditions of the experiment.
#'		# final_count should be the number of columns associated with
#'		#       the final conditions of the experiment.
#'		working_sets <- extract_working_sets(csv_data, init_count, final_count)
##########
extract_working_sets <- function (first_my_data, replicate_count1, replicate_count2) {
	# first column of dataset is probes
	probe_list <- subset(first_my_data, select = c(1))
	
	column_set1 <- replicate_count1 +1
			#uses x value to designate last column of first condition results 
			#(add one to skip over probe column)

	column_set2 <- replicate_count2 +2
			#uses x value to designate first column of second condition results

	column_set3 <- column_set1 + replicate_count2
			#uses y value to designate last column of second condition results

	initial_conditions <- first_my_data[,2:column_set1]
			#create subset of first condition columns

	final_conditions <- first_my_data[,column_set2:column_set3]
			#create subset of second condition columns
    	   
    # output = probes, initial_conditions, final_conditions
    return(list(probes = probe_list, initial = initial_conditions, final = final_conditions))
}

##########
#' Converts a data.frame of probes with initial and final
#' experimental conditions to a table with probes and
#' fold change values.
#'
#' @export
#' @title calculates fold values from a table of experimental data
#' @param probes the data set of probes.
#' @param initial_conditions a data set of replicates
#'		corresponding to initial conditions.
#' @param final_conditions a data set of replicates
#'		corresponding to final conditions.
#' @return a data.frame comprised of probe names and fold-change values
#' @examples
#'		# read in the EcoliMutMA sample data from the package
#'		data(EcoliMutMA, package="ELBOW")
#'		csv_data <- EcoliMutMA
#'		# - OR - Read in a CSV file (uncomment - remove the #'s
#'		#        - from the line below and replace 'filename' with
#'		#        the CSV file's filename)
#'		# csv_data <- read.csv(filename)
#'		
#'		# set the number of initial and final condition replicates both to three
#'		init_count  <- 3
#'		final_count <- 3
#'	
#'		# Parse the probes, intial conditions and final conditions
#'		# out of the CSV file.  Please see: extract_working_sets
#'		# for more information.
#'		#
#'		# init_count should be the number of columns associated with
#'		#       the initial conditions of the experiment.
#'		# final_count should be the number of columns associated with
#'		#       the final conditions of the experiment.
#'		working_sets <- extract_working_sets(csv_data, init_count, final_count)
#'		
#'		probes <- working_sets[[1]]
#'		initial_conditions <- working_sets[[2]]
#'		final_conditions <- working_sets[[3]]
#'		
#'		# Uncomment to output the plot to a PNG file (optional)
#'		# png(file="output_plot.png")
#'
#'		my_data <- replicates_to_fold(probes, initial_conditions, final_conditions)
##########
replicates_to_fold <- function(probes, initial_conditions, final_conditions) {
	#average of first condition columns
	averageA <- rowSums(initial_conditions)/ncol(initial_conditions)

	#average of second condition columns
	averageB <- rowSums(final_conditions)/ncol(final_conditions)
	
	print("rowsums")

	#fold values over all sets
	fold <- averageA-averageB

	print("fold")
	print(head(fold))

	#now have a dataset of fold values between conditions and matching probes
	first_my_data  <- data.frame(cbind(probes,fold))
	second_my_data <- data.frame(cbind(probes,fold))

	print("bound data")
	print(head(second_my_data))
	print("firsta_data")
	print(head(first_my_data))

	my_data <- data.frame(second_my_data[with(first_my_data,order(fold)),]) #sorted from lowest value to highest value

	print("sorted")
    
	headers <- list("probe_names","fold") #creates a set of header names that is always consistent

	print("headers")
    
	#step one of transfering column names
	col_names <-  colnames(my_data, do.NULL = TRUE, prefix = "col") 

	# headers are now column names for output files
	# if you get a message saying "fold" not found restart R
	colnames(my_data) <- headers 

	return(my_data)
}

##########
#' calculates the upper and lower elbow cut-off fold-values
#'
#' @export
#' @title calculates the upper and lower elbow cut-off fold-values
#' @param fold_values a data.frame containing all of the fold
#'		values to calculate the elbow for.
#' @return a list containing the following:
#'		\itemize{
#'			\item{up_limit --- the upper (most positive) limit of the elbow curve}
#'			\item{low_limit --- the lower (most negative) limit of the elbow curve}
#'		}
#' @examples
#'		# read in the EcoliMutMA sample data from the package
#'		data(EcoliMutMA, package="ELBOW")
#'		csv_data <- EcoliMutMA
#'		# - OR - Read in a CSV file (uncomment - remove the #'s
#'		#        - from the line below and replace 'filename' with
#'		#        the CSV file's filename)
#'		# csv_data <- read.csv(filename)
#'		
#'		# set the number of initial and final condition replicates both to three
#'		init_count  <- 3
#'		final_count <- 3
#'	
#'		# Parse the probes, intial conditions and final conditions
#'		# out of the CSV file.  Please see: extract_working_sets
#'		# for more information.
#'		#
#'		# init_count should be the number of columns associated with
#'		#       the initial conditions of the experiment.
#'		# final_count should be the number of columns associated with
#'		#       the final conditions of the experiment.
#'		working_sets <- extract_working_sets(csv_data, init_count, final_count)
#'		
#'		probes <- working_sets[[1]]
#'		initial_conditions <- working_sets[[2]]
#'		final_conditions <- working_sets[[3]]
#'		
#'		# Uncomment to output the plot to a PNG file (optional)
#'		# png(file="output_plot.png")
#'
#'		my_data <- replicates_to_fold(probes, initial_conditions, final_conditions)
#'	
#'		# compute the elbow for the dataset
#'		limits <- do_elbow(data.frame(my_data$fold))
##########
do_elbow <- function (fold_values) {
	
	# used to divide the data into two creating an upper and lower set of values
	halfway_point <- round(nrow(fold_values)/2,digits=0)

	#counts the rows in the dataset
	subx <- nrow(fold_values)

	# gives all upper y ordered from lowest to highest values starting at the midpoint
	yupper <- fold_values[halfway_point:subx,]

	# count the number in the upper set 
	# (if the dataset is symmetric about zero this will be all positive values)
	xupper1 <- length(yupper)

	# creating an x vector for graphing and predicting
	xupper <- 1:xupper1

	# takes the upper half of the fold data and makes a smooth curve of it.
	splupper <- smooth.spline(xupper, yupper)

	# gives a smoothed line of the y data
	predupper1 <- predict(splupper)

	##########
	# We are now getting the derivatives at each point on
	# the smoothed curve made from the fold data
	#   diff(predupper$y = dy/dx because dx = 1)
	#     result is a long flat line with noise at the end
	# The noise at the end corresponds to the elbow of the
	# curve for cut off values
	##########
	predupper <- data.matrix(data.frame(predupper1[2]))
	predupperde <- diff(predupper) 

	# lose the last fold value during diff. To keep row
	# numbers consistent we add one value
	predupperadd1 <- max(predupperde)

	# new extra value is added
	predupperadd <- c(predupperde,predupperadd1)

	# creates table of fold values and diff values of smooth splined curve
	combupperdata <- data.frame(cbind(yupper,predupperadd))

	##########
	# We take the mean of the entire diff smooth spline predicted curve. 
	# The mean value is very close to value of the flat portion of the
	# curve since most of them are on the flat 1.75*mean marks the beginning
	# of the noise at the end of the predposde result
	##########
	upper_limit1 <- subset(combupperdata, predupperadd >= 1.75*mean(combupperdata$predupperadd))

	# pulls out the value of the fold change = to the point where noice
	# begins, i.e. the elbow
	upper_limit <- upper_limit1[1,1] 

	# gives all lower y ordered from lowest to highest values
	# (i.e. starting at next to mid point)
	ylower1 <- fold_values[1:halfway_point,]
	
	# flipping data for the smooth spline function for making a curve
	ylower <- sort(ylower1, decreasing = TRUE)

	# creating an x value for number of rows
	xlower1 <-length(ylower)

	# creating an x vector
	xlower <- 1:xlower1

	# takes the fold data and makes a smooth curve of it.
	spllower <- smooth.spline(xlower, ylower)

	# gives a smoothed line of the y data
	predlower1 <- predict(spllower)

	##########
	# We are now getting the derivatives at each point on
	# the smoothed curve made from the fold data
	#   diff(predlower$y = dy/dx because dx = 1)
	#     result is a long flat line with noise at the end
	# The noise at the end corresponds to the elbow of the
	# curve for cut off values
	##########
	predlower <- data.matrix(data.frame(predlower1[2]))
	predlowerde <- diff(predlower) 

	# We lose the last fold value during diff. To keep row
	# numbers consistent we add one value
	predloweradd1 <- min(predlowerde)

	# New extra value is added
	predloweradd <- c(predlowerde,predloweradd1)

	# Table of fold values and diff values of smooth splined curve
	comblower <- cbind(ylower,predloweradd) 

	# need to dataframe for future steps
	comblowerdata <- data.frame(comblower)

	##########
	# We take the mean of the entire diff smooth spline predicted curve. 
	# The mean value is very close to value of the flat portion of the
	# curve since most of them are on the flat 1.75*mean marks the beginning
	# of the noise at the end of the predposde result = curve
	##########
	lower_limit1 <- subset(comblowerdata, predloweradd <= 1.75*mean(comblowerdata$predloweradd))

	# pulls out the value of the fold change = to
	# the point where noice begins, i.e. the elbow
	lower_limit <- lower_limit1[1,1]
	return(list(up_limit = upper_limit, low_limit = lower_limit))
}


##########
#' The following function determines the variation of cut off
#' limits for significance between individual subsets. The first
#' step is to generate the subsets, then the function determines
#' the highest and lowest values for the upper limit and the
#' highest and lowest values for the lower limit. These values
#' are used to determine and upper and lower error value reflecting
#' fold dataset variance.
#'
#' @export
#' @title calculates the variance for the upper and lower limits of an Elbow curve
#' @param probes the data set of probes.
#' @param initial_conditions a data set of replicates
#'		corresponding to initial conditions.
#' @param final_conditions a data set of replicates
#'		corresponding to final conditions.
#' @return a list containing the following keys:
#'		\itemize{
#'			\item{max_upper --- the maximum upper elbow limit \emph{(most positive)}}
#'			\item{min_upper --- the minimum upper elbow limit}
#'			\item{max_upper --- the maximum lower elbow limit}
#'			\item{min_upper --- the minimum lower elbow limit \emph{(most negative)}}
#'		}
#' @examples
#'		# read in the EcoliMutMA sample data from the package
#'		data(EcoliMutMA, package="ELBOW")
#'		csv_data <- EcoliMutMA
#'		# - OR - Read in a CSV file (uncomment - remove the #'s
#'		#        - from the line below and replace 'filename' with
#'		#        the CSV file's filename)
#'		# csv_data <- read.csv(filename)
#'		
#'		# set the number of initial and final condition replicates both to three
#'		init_count  <- 3
#'		final_count <- 3
#'	
#'		# Parse the probes, intial conditions and final conditions
#'		# out of the CSV file.  Please see: extract_working_sets
#'		# for more information.
#'		#
#'		# init_count should be the number of columns associated with
#'		#       the initial conditions of the experiment.
#'		# final_count should be the number of columns associated with
#'		#       the final conditions of the experiment.
#'		working_sets <- extract_working_sets(csv_data, init_count, final_count)
#'		
#'		probes <- working_sets[[1]]
#'		initial_conditions <- working_sets[[2]]
#'		final_conditions <- working_sets[[3]]
#'		
#'		# Uncomment to output the plot to a PNG file (optional)
#'		# png(file="output_plot.png")
#'
#'		my_data <- replicates_to_fold(probes, initial_conditions, final_conditions)
#'	
#'		# compute the elbow for the dataset
#'		limits <- do_elbow(data.frame(my_data$fold))
#'	
#'		plus_minus <- elbow_variance(probes, initial_conditions, final_conditions)
##########
elbow_variance <- function(probes, initial_conditions, final_conditions) {	
	# designating terms of final vector
	out <- vector("list", length = ncol(initial_conditions) * ncol(final_conditions)) 
	    
	counter <- 1
	for (i in 1:ncol(initial_conditions)) {
		for (j in 1:ncol(final_conditions)) {
			out[[counter]] <- cbind(initial_conditions[,i, drop = FALSE], final_conditions[,j, drop = FALSE])
			counter <- counter + 1
		}
	}
	
	# create a dataset where the elbow limits have been
	# determined over all replicate pairs.
	list_dataframes <- lapply(out, function(x) {
		my_data <- replicates_to_fold(probes, x[1], x[2])
		
		# pulls out the fold data
		# (making it a dataframe avoids all kinds of R bugs)
		fold_values <- data.frame(my_data[,"fold"])
		headers <- list("fold") #creates a set of header names that is always consistent
		    
		#step one of transfering column names
		col_names <-  colnames(fold_values, do.NULL = TRUE, prefix = "col") 
		
		# headers are now column names for output files
		# if you get a message saying "fold" not found restart R
		colnames(fold_values) <- headers
		
		return(do_elbow(fold_values))
	} ) 
	
	# when you give up - do.call will do this - it is in the package "Base"
	#print(list_dataframes)
	row_vector <- as.data.frame(matrix(unlist(list_dataframes), nrow = length(list_dataframes), byrow = TRUE))

	# determine the pair with the maximum number of probes
	upper_plus_value <- max(row_vector[1]) 
	lower_plus_value <- min(row_vector[1]) 
	
	# determine the pair with the minimum number of probes
	upper_minus_value <- max(row_vector[2])
	lower_minus_value <- min(row_vector[2])
	
	# max_upper = the maximum upper elbow limit (most positive)
	# min_upper = the minimum upper elbow limit
	# max_lower = the maximum lower elbow limit
	# min_lower = the minimum lower elbow limit (most negative)
	return(list(max_upper = upper_plus_value, min_upper = lower_plus_value, max_lower = upper_minus_value, min_lower = lower_minus_value))
}


##########
#' determine an upper and lower error limit using all
#' possible subsets of differences between initial
#' conditions and then use the median value for each
#' probe to create the null set. Maximum and minimum
#' values are used to set upper and lower error bounds.
#'
#' @export
#' @title determines the variance for a null Elbow curve
#' @param my_data A table to analyze the null variance of.
#'		The columns in the table should be as follows:
#'		\itemize{
#'			\item{probes --- one column containing the names of the probes}
#'			\item{initial conditions --- one or more columns containing the
#'				  initial conditions of the experiment}
#'			\item{final conditions --- one or more columns containing the
#'				  final conditions of the experiment}
#'		}
#' @param upper_limit the actual upper limit cut-off for the Elbow curve
#' @param lower_limit the actual lower limit cut-off for the Elbow curve
#' @param initial_conditions a data set of replicates
#'		corresponding to initial conditions.
#' @return A list containing the following keys:
#'		\itemize{
#'			\item{prmax --- the error limit based on the maximum value for each probe}
#'			\item{prmed --- the null (median) value for each probe}
#'			\item{prmin --- the error limit based on the minimum value for each probe}
#'		}
#' @examples
#'		# read in the EcoliMutMA sample data from the package
#'		data(EcoliMutMA, package="ELBOW")
#'		csv_data <- EcoliMutMA
#'		# - OR - Read in a CSV file (uncomment - remove the #'s
#'		#        - from the line below and replace 'filename' with
#'		#        the CSV file's filename)
#'		# csv_data <- read.csv(filename)
#'		
#'		# set the number of initial and final condition replicates both to three
#'		init_count  <- 3
#'		final_count <- 3
#'	
#'		# Parse the probes, intial conditions and final conditions
#'		# out of the CSV file.  Please see: extract_working_sets
#'		# for more information.
#'		#
#'		# init_count should be the number of columns associated with
#'		#       the initial conditions of the experiment.
#'		# final_count should be the number of columns associated with
#'		#       the final conditions of the experiment.
#'		working_sets <- extract_working_sets(csv_data, init_count, final_count)
#'		
#'		probes <- working_sets[[1]]
#'		initial_conditions <- working_sets[[2]]
#'		final_conditions <- working_sets[[3]]
#'		
#'		# Uncomment to output the plot to a PNG file (optional)
#'		# png(file="output_plot.png")
#'
#'		my_data <- replicates_to_fold(probes, initial_conditions, final_conditions)
#'	
#'		# compute the elbow for the dataset
#'		limits <- do_elbow(data.frame(my_data$fold))
#'	
#'		plus_minus <- elbow_variance(probes, initial_conditions, final_conditions)
#'	
#'		max_upper_variance <- plus_minus$max_upper
#'		min_upper_variance <- plus_minus$min_upper
#'		max_lower_variance <- plus_minus$max_lower
#'		min_lower_variance <- plus_minus$min_lower
#'	
#'		# rounded number for nice appearance graph
#'		upper_limit <- round(limits[[1]],digits=2)
#'	
#'		# rounded number nice appearance for graph
#'		lower_limit <- round(limits[[2]],digits=2)
#'	
#'		p_limits <- null_variance(my_data, upper_limit, lower_limit, initial_conditions)
#'	
#'		prowmin <- p_limits[[1]]
#'		prowmax <- p_limits[[2]]
#'		prowmedian <- p_limits[[3]]
##########
null_variance <- function (my_data, upper_limit, lower_limit, initial_conditions) {
	initial_condition_count <- ncol(initial_conditions)
	
	# significant probes and folds based on elbow limits
	sig <- subset(my_data, my_data$fold >= upper_limit | my_data$fold <= lower_limit) 
	
	# counts the rows in the dataset
	subx <- nrow(my_data)
	   
	# creates a vector of numbers that match the number of rows in y
	probe_order <- 1:subx
	
	#original_cond5 <- combn(initial_conditions,2) #make all possible pairs of columns excluding identical self/self pairs
	
	pcolumn_set1 <- round((initial_condition_count/2),digits=0)
	
	# get length of one half of initial_conditions, rounded to avoid 0.5 values
	    
	# uses pcolumn_set1 value to designate first column of second half initial_conditions
	pcolumn_set2 <- pcolumn_set1 +1
	    
	# uses y value to designate last column of second condition results
	pcolumn_set3 <- initial_condition_count - pcolumn_set1
	
	# create subset of first condition columns    
	initial_conditionsa <- data.frame(initial_conditions[,1:pcolumn_set1])
	
	# create subset of second condition columns
	initial_conditionsb <- data.frame(initial_conditions[pcolumn_set2:initial_condition_count])
	
	# designating terms of final vector
	pout <- vector("list", length = ncol(initial_conditionsa) * ncol(initial_conditionsb)) 
	    
	counter <- 1
	for (i in 1:ncol(initial_conditionsa)) {
	    for (j in 1:ncol(initial_conditionsb)) {
			pout[[counter]] <- cbind(initial_conditionsa[,i, drop = FALSE], initial_conditionsb[,j, drop = FALSE])
			counter <- counter + 1
		}
	}
	
	#function determines the fold value for all pairs
	pvalue_variance <- function(pset){
	    
	    pfold <- pset[,1] - pset[,2]
	    
	    pfold_sort <- sort(pfold, decreasing = FALSE)
	    
	}
	
	# creates a list of fold values for all possible
	# combination of pairs from original condition replicates  
	plist <- lapply(pout, pvalue_variance) 
	
	prows <- do.call(cbind,plist)
	
	# gets a median
	prowmedian <- apply(prows,1,median)
	
	# Takes the max of pos values and min of neg values 
	prowmax <- apply(prows,1, max)
	
	#Takes the min of pos values and max of neg values 
	prowmin <- apply(prows,1, min)
	
	return(list(prmin = prowmin, prmax = prowmax, prmed = prowmedian))
}


##########
#' The following function determines a p value from
#' the comparison of the slope of a regression line
#' (i.e. logit).  The function calculates a fold value
#' from a set of initial conditions (replicate set 1)
#' substracted from a set of final conditions
#' (replicate set 2).  The glm function is used to
#' assess the fit of the dataset to the model, assuming
#' family = binomial (link = "probit"), maxit = 1000.
#' Determining the \eqn{\log\chi^2} p-value for the model,
#' null equals no significance for the set of significant
#' probes to the model, in which case the Elbow method
#' cannot be used to assess the dataset.  If significant
#' probes are important to the model, then the \eqn{\log\chi^2}
#' p-value will be below 0.05, in which case the Elbow method
#' can be used to assess the dataset.
#'
#' @export
#' @title The calculated \eqn{\log\chi^2} p-value for the fit
#'			of the elbow curve to the model.
#' @param my_data A table (data.frame) to analyze the
#'		elbow variance of. The columns in the table
#'		should be as follows:
#'		\itemize{
#'			\item{probes --- one column containing the names of the probes}
#'			\item{fold --- the fold values for the table}
#'		}
#' @param upper_limit the upper limit/cut-off for the elbow.
#' @param lower_limit the lower limit/cut-off for the elbow.
#' @return The calculated \eqn{\log\chi^2} p-value for the elbow curve.
#' @examples
#'		# read in the EcoliMutMA sample data from the package
#'		data(EcoliMutMA, package="ELBOW")
#'		csv_data <- EcoliMutMA
#'		# - OR - Read in a CSV file (uncomment - remove the #'s
#'		#        - from the line below and replace 'filename' with
#'		#        the CSV file's filename)
#'		# csv_data <- read.csv(filename)
#'		
#'		# set the number of initial and final condition replicates both to three
#'		init_count  <- 3
#'		final_count <- 3
#'	
#'		# Parse the probes, intial conditions and final conditions
#'		# out of the CSV file.  Please see: extract_working_sets
#'		# for more information.
#'		#
#'		# init_count should be the number of columns associated with
#'		#       the initial conditions of the experiment.
#'		# final_count should be the number of columns associated with
#'		#       the final conditions of the experiment.
#'		working_sets <- extract_working_sets(csv_data, init_count, final_count)
#'		
#'		probes <- working_sets[[1]]
#'		initial_conditions <- working_sets[[2]]
#'		final_conditions <- working_sets[[3]]
#'		
#'		# Uncomment to output the plot to a PNG file (optional)
#'		# png(file="output_plot.png")
#'
#'		my_data <- replicates_to_fold(probes, initial_conditions, final_conditions)
#'	
#'		# compute the elbow for the dataset
#'		limits <- do_elbow(data.frame(my_data$fold))
#'	
#'		plus_minus <- elbow_variance(probes, initial_conditions, final_conditions)
#'	
#'		max_upper_variance <- plus_minus$max_upper
#'		min_upper_variance <- plus_minus$min_upper
#'		max_lower_variance <- plus_minus$max_lower
#'		min_lower_variance <- plus_minus$min_lower
#'	
#'		# rounded number for nice appearance graph
#'		upper_limit <- round(limits[[1]],digits=2)
#'	
#'		# rounded number nice appearance for graph
#'		lower_limit <- round(limits[[2]],digits=2)
#'	
#'		p_limits <- null_variance(my_data, upper_limit, lower_limit, initial_conditions)
#'	
#'		prowmin <- p_limits[[1]]
#'		prowmax <- p_limits[[2]]
#'		prowmedian <- p_limits[[3]]
#'	
#'		pvalue1 <- get_pvalue(my_data, upper_limit, lower_limit)
#'
##########
get_pvalue <- function(my_data, upper_limit, lower_limit) {
	subx <- nrow(my_data)
	
	probe_order <- 1:subx
	
	# creates a true false answer for sig versus not sig,
	# i.e. sig = True, not sig = false
	set1a <- my_data$fold >= upper_limit|my_data$fold <= lower_limit 

	#data.frame step to help R work
	set1df <- data.frame(set1a)
	
	# converts true to 1 and false to zero for logistic fitting
	set1df[set1df == TRUE] <- 1
	
	# analysis dataset
	set1 <- cbind(set1df,probe_order,my_data$fold)
	
	fit1 <- glm(set1a~my_data$fold + probe_order,family=binomial(link = "probit"), maxit = 1000, set1) 
    
	pvalue1 <- with(fit1, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
	return(pvalue1)
}


##########
#' Plots an elbow curve and its associated data:
#'	\itemize{
#'		\item{the upper and lower elbow limits for the curve}
#'		\item{the upper, lower, and median initial condition Elbow plots}
#'		\item{the \eqn{\log\chi^2} p-value for the Elbow curve}
#'		\item{the variance for the upper and lower Elbow cut-off values}
#'	}
#'
#' @export
#' @title plots Elbow curve data.
#' @param my_data A table (data.frame) to plot.
#'		The columns in the table should be as follows:
#'		\itemize{
#'			\item{probes --- one column containing the names of the probes}
#'			\item{fold --- the fold values for the table}
#'		}
#' @param upper_limit the upper limit/cut-off for the elbow.
#' @param lower_limit the lower limit/cut-off for the elbow.
#' @param pvalue1 the \eqn{\log\chi^2} p-value for the elbow curve.
#' @param prowmin the error limit based on the maximum value for each probe.
#' @param prowmax the error limit based on the minimum value for each probe.
#' @param prowmedian the null (median) value for each probe.
#' @param max_upper_variance the maximum upper elbow limit \emph{(most positive)}.
#' @param min_upper_variance the minimum upper elbow limit.
#' @param max_lower_variance the maximum lower elbow limit.
#' @param min_lower_variance the minimum lower elbow limit \emph{(most negative)}.
#' @param gtitle the title to display for the graph.
#' @examples
#'		# read in the EcoliMutMA sample data from the package
#'		data(EcoliMutMA, package="ELBOW")
#'		csv_data <- EcoliMutMA
#'		# - OR - Read in a CSV file (uncomment - remove the #'s
#'		#        - from the line below and replace 'filename' with
#'		#        the CSV file's filename)
#'		# csv_data <- read.csv(filename)
#'		
#'		# set the number of initial and final condition replicates both to three
#'		init_count  <- 3
#'		final_count <- 3
#'	
#'		# Parse the probes, intial conditions and final conditions
#'		# out of the CSV file.  Please see: extract_working_sets
#'		# for more information.
#'		#
#'		# init_count should be the number of columns associated with
#'		#       the initial conditions of the experiment.
#'		# final_count should be the number of columns associated with
#'		#       the final conditions of the experiment.
#'		working_sets <- extract_working_sets(csv_data, init_count, final_count)
#'		
#'		probes <- working_sets[[1]]
#'		initial_conditions <- working_sets[[2]]
#'		final_conditions <- working_sets[[3]]
#'		
#'		# Uncomment to output the plot to a PNG file (optional)
#'		# png(file="output_plot.png")
#'
#'		my_data <- replicates_to_fold(probes, initial_conditions, final_conditions)
#'	
#'		# compute the elbow for the dataset
#'		limits <- do_elbow(data.frame(my_data$fold))
#'	
#'		plus_minus <- elbow_variance(probes, initial_conditions, final_conditions)
#'	
#'		max_upper_variance <- plus_minus$max_upper
#'		min_upper_variance <- plus_minus$min_upper
#'		max_lower_variance <- plus_minus$max_lower
#'		min_lower_variance <- plus_minus$min_lower
#'	
#'		# rounded number for nice appearance graph
#'		upper_limit <- round(limits[[1]],digits=2)
#'	
#'		# rounded number nice appearance for graph
#'		lower_limit <- round(limits[[2]],digits=2)
#'	
#'		p_limits <- null_variance(my_data, upper_limit, lower_limit, initial_conditions)
#'	
#'		prowmin <- p_limits[[1]]
#'		prowmax <- p_limits[[2]]
#'		prowmedian <- p_limits[[3]]
#'	
#'		pvalue1 <- get_pvalue(my_data, upper_limit, lower_limit)
#'		plot_elbow(my_data, upper_limit, lower_limit, pvalue1, prowmin, prowmax, prowmedian, max_upper_variance, min_upper_variance, max_lower_variance, min_lower_variance, "Title of the ELBOW Plot")
#'
##########
plot_elbow <- function(my_data, upper_limit, lower_limit, pvalue1, prowmin, prowmax, prowmedian, max_upper_variance, min_upper_variance, max_lower_variance, min_lower_variance, gtitle = "") {
	subx <- nrow(my_data)
	
	fold_values <- my_data$fold
	
	probe_order <- 1:subx
	
	# make p-value 3 decimal digits
	pvalue <- signif(pvalue1, digits = 3)
	
	# rounded number for nice appearance graph
	upper_limit <- round(upper_limit,digits=2)
	lower_limit <- round(lower_limit,digits=2)
	min_upper_variance <- round(min_upper_variance,digits=2)
	max_upper_variance <- round(max_upper_variance,digits=2)
	max_lower_variance <- round(max_lower_variance,digits=2)
	min_lower_variance <- round(min_lower_variance,digits=2)
	
	# creation of error bounds text message
	error_message <- paste(upper_limit,min_upper_variance,max_upper_variance,lower_limit,max_lower_variance,min_lower_variance) 
	
	# creation of error bounds text message
	error_message1 <- paste("upper elbow limit =",upper_limit,"(replicate variance error",min_upper_variance," to ",max_upper_variance,")")
	
	# creation of error bounds text message
	error_message2 <- paste("lower elbow limit =",lower_limit,"( replicate variance error",max_lower_variance," to ",min_lower_variance,")")
	
	color_message <- paste("red=elbow, dash bl=median initial conds (gr=upper,pur=lower)")

	pvalue_message <- if(pvalue1 < 1) {paste("log chi squared p =", pvalue)} else if(pvalue1 > 1) {paste("log chi squared p =", pvalue)}else if (pvalue1 == 1){paste("p value cannot be computed due to initial condition variance - consider testing density/distribution of individual replicates for outliers")}	
	
	message4 <- paste("Report with a list of significant probes has been generated and labeled 'Significant_Probes.csv.'")
	
	# plots the original uncorrected data
	plot(probe_order,fold_values, type='l', lwd=2, col='red', main = "Elbow Plot Results") 
	
	# finds the midpoint on the x axis for placing the text
	text_x <- subx/2
	
	# This is for creating an upper input line by the data 
	text_yupper <- upper_limit +0.5
	
	# This is for creating an lower input line by the data
	text_ylower <- lower_limit +0.5
	
	# places x coordinate for text message
	error_messagex <- subx/12
	
	error_message1y <- min(fold_values) - 3*min(fold_values)/10
	
	# places y coordinate for text message
	error_message2y <- min(fold_values) - 2*min(fold_values)/10
	
	# places y coordinate for text message
	message3y <- min(fold_values) - min(fold_values)/10
	
	# places y coordinate for p_vaue message message
	filey <- max(fold_values)
	
	text_p_valuey <- min(fold_values)
	
	# places lines at the corrected upper and lower values
	abline(h = lower_limit, untf = TRUE)
	abline(h = upper_limit, untf = TRUE)
	
	# adds a second line representing the average value of
	# the variance from original condition, replicate on
	# order to be able to compare that line to the elbow line
	lines(probe_order,prowmax,type='l', lty = 3, lwd=2, col = "green")
	lines(probe_order,prowmin,type='l', lty = 3, lwd=2, col = "purple")
	lines(probe_order,prowmedian,type='l', lty = 'dashed', lwd=2)
	
	#places text
	text(error_messagex,error_message1y, error_message1, pos = 4)
	text(error_messagex,error_message2y, error_message2, pos = 4)
	text(error_messagex,message3y, color_message, pos = 4)
	text(error_messagex,filey,gtitle, pos = 4)
	text(error_messagex,text_p_valuey, pvalue_message, pos = 4)
	
	print(error_message1)
	print(error_message2)
	print(pvalue_message)
}

##########
#' Analyzes:
#'	\itemize{
#'		\item{the Elbow cut-offs}
#'		\item{the Elbow curve variance}
#'		\item{the upper and lower error Elbow curves}
#'		\item{\eqn{\log\chi^2} p-value for the Elbow curve model}
#'	}
#'
#' Then plots the data to a curve - AND - prints the statistics
#' to both the terminal and the plot canvas.
#'
#' @export
#' @title extracts all elbow statistics and plots and elbow curve.
#' @param probes the data set of probes (as a data.frame).
#' @param initial_conditions a 2D data.fram containing
#'		all of the replicate gene expression values
#'		corresponding to the experiment's initial conditions.
#' @param final_conditions a 2D data.fram containing
#'		all of the replicate gene expression values
#'		corresponding to the experiment's final conditions.
#' @param gtitle the title to display for the graph.
#' @return a list of all significant probes
#' @examples
#'		# read in the EcoliMutMA sample data from the package
#'		data(EcoliMutMA, package="ELBOW")
#'		csv_data <- EcoliMutMA
#'		# - OR - Read in a CSV file (uncomment - remove the #'s
#'		#        - from the line below and replace 'filename' with
#'		#        the CSV file's filename)
#'		# csv_data <- read.csv(filename)
#'		
#'		# set the number of initial and final condition replicates both to three
#'		init_count  <- 3
#'		final_count <- 3
#'	
#'		# Parse the probes, intial conditions and final conditions
#'		# out of the CSV file.  Please see: extract_working_sets
#'		# for more information.
#'		#
#'		# init_count should be the number of columns associated with
#'		#       the initial conditions of the experiment.
#'		# final_count should be the number of columns associated with
#'		#       the final conditions of the experiment.
#'		working_sets <- extract_working_sets(csv_data, init_count, final_count)
#'		
#'		probes <- working_sets[[1]]
#'		initial_conditions <- working_sets[[2]]
#'		final_conditions <- working_sets[[3]]
#'		
#'		# Uncomment to output the plot to a PNG file (optional)
#'		# png(file="output_plot.png")
#'		
#'		# Analyze the elbow curve.
#'		sig <- analyze_elbow(probes, initial_conditions, final_conditions)
#'		
#'		# write the significant probes to 'signprobes.csv'
#'		write.table(sig,file="signprobes.csv",sep=",",row.names=FALSE)
##########
analyze_elbow <- function(probes, initial_conditions, final_conditions, gtitle = "") {
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
	
	# create a table with the pos and neg probes at the significant limits
	sig <- subset(my_data, my_data$fold >= upper_limit | my_data$fold <= lower_limit) 
	
	p_limits <- null_variance(my_data, upper_limit, lower_limit, initial_conditions)
	
	prowmin <- p_limits[[1]]
	prowmax <- p_limits[[2]]
	prowmedian <- p_limits[[3]]
	
	#count1 <- subset(regression_dftf2,regression_dftf2$regression_tf2 ==F)
	
	#prowminpos <- subset(prowmin, 1, prowmin > 0 ) #separates the min of pos values
	
	#prowminneg <- subset(prowmin, prowmin < 0) #Separates the maximum values of negative
	
	pvalue1 <- get_pvalue(my_data, upper_limit, lower_limit)
	
	plot_elbow(my_data, upper_limit, lower_limit, pvalue1, prowmin, prowmax, prowmedian, max_upper_variance, min_upper_variance, max_lower_variance, min_lower_variance, gtitle)
	
	return(sig)
}

##########
#' Performs the ELBOW fold change test on an MArrayLM
#' object.  This is a wrapper class to help integrate the
#' ELBOW method into Bioconductor.
#' followed tutorial from: http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
#'
#' Some of the code in this method is based on how toptable accesses
#' the MArrayLM object to read fold values.  Therefore, any MArrayLM
#' object which works with the toptable method should also work with
#' this method.
#'
#' @export
#'
#' @title Calculates fold values from an MArrayLM object.
#'
#' @param marraylm is the MArrayLM object to analyze.
#'
#' @param columns is the list of sample columns to obtain the elbow
#'                fold cut-off values for.  This can be specified as
#'                a vector or a single value.
#'
#' @return a matrix specified as follows
#'		\itemize{
#'			\item{columns --- (1) \dQuote{up_limit}, the upper
#'                                ELBOW fold-change cut-off value;
#'			                  (2) \dQuote{low_limit}, the lower
#'                                ELBOW fold-change cut-off value}
#'			\item{rows --- one row per sample, specified by the
#'                         parameter \dQuote{columns.}}
#'		}
#' @examples
#'		########
#'		# LOAD DATA INTO LIMMA
#'		########
#'		library("limma")
#'		
#'		# load a filtered expression set into R
#'		# NOTE: see the vignette for instructions on preparing
#'		#       a filtered dataset with your own data.
#'		data(GSE20986_eset_exprs, package="ELBOW")
#'		data(GSE20986_design, package="ELBOW")
#'
#'		# fit the linear model to the filtered expression set
#'		fit <- lmFit(GSE20986_eset_exprs, GSE20986_design)
#'		
#'		# set up a contrast matrix to compare tissues v cell line
#'		contrast.matrix <- makeContrasts(huvec_choroid = huvec - choroid, huvec_retina = huvec - retina, huvec_iris <- huvec - iris, levels=GSE20986_design)
#'		
#'		# Now the contrast matrix is combined with the per-probeset linear model fit.
#'		huvec_fits <- contrasts.fit(fit, contrast.matrix)
#'		huvec_ebFit <- eBayes(huvec_fits)
#'		
#'		########
#'		# GET THE ELBOW LIMIT (this function)
#'		########
#'		get_elbow_limma(huvec_ebFit)
#'		
##########
get_elbow_limma <- function(marraylm, columns=NULL) {
	dataset <- marraylm$coefficients
	cnames  <- colnames(dataset)
	results <- NULL
	
	if (is.null(columns)) {
		columns <- 1:length(cnames)
	}
	
	for (cnum in unique(columns)) {
		coeff <- dataset[,cnum]
		matrixcoeff <- as.matrix(coeff)
		fold_values <- sort(matrixcoeff)
		elbow_data <- do_elbow(as.matrix(fold_values))
		if (is.null(results)) {
			results <- matrix(c(elbow_data$up_limit, elbow_data$low_limit), ncol=2)
		} else {
			results <- rbind(results, c(elbow_data$up_limit, elbow_data$low_limit))
		}
	}
	
	if (!is.null(results)) {
		rownames(results) <- cnames[columns]
		colnames(results) <- c("up_limit", "low_limit")
		return(results)
	}
}


##########
#' Performs the ELBOW fold change test on an CountDataSet
#' from DESeq object.  This is a wrapper class to help
#' integrate the ELBOW method into Bioconductor.
#' followed tutorial from: http://cgrlucb.wikispaces.com/Spring+2012+DESeq+Tutorial
#'
#' @export
#'
#' @title Calculates fold values from an MArrayLM object.
#'
#' @param rnaSeq is the CountDataSet object to analyze.
#'
#' @return a matrix specified as follows
#'		\itemize{
#'			\item{columns --- (1) \dQuote{up_limit}, the upper
#'                                ELBOW fold-change cut-off value;
#'			                  (2) \dQuote{low_limit}, the lower
#'                                ELBOW fold-change cut-off value}
#'			\item{rows --- one row per sample, specified by the
#'                         parameter \dQuote{columns.}}
#'		}
#'
#' @examples
#'		# install the DESeq libraries
#'		#source("http://www.bioconductor.org/biocLite.R")
#'		#biocLite("DESeq")
#'
#'		## download the table
#'		library("DESeq")
#'
#'		# the following bam file dataset was obtained from:
#'		# http://cgrlucb.wikispaces.com/file/view/yeast_sample_data.txt
#'		# it has been downloaded into this package for speed convenience.
#'		filename <- system.file("extdata", "yeast_sample_data.txt", package = "ELBOW")
#'	
#'		count_table <- read.table(filename, header=TRUE, sep="\t", row.names=1)
#'		expt_design <- data.frame(row.names = colnames(count_table), condition = c("WE","WE","M","M","M"))
#'		conditions = expt_design$condition
#'		data <- newCountDataSet(count_table, conditions)
#'		data <- estimateSizeFactors(data)
#'		data <- as(data, "CountDataSet")
#'		## data <- estimateVarianceFunctions(data)
#'		data <- estimateDispersions(data)
#'
#'		# this next step is essential, but it takes a long time...
#'		# so, just like a good cooking show we will skip this step
#'		# and load a finished version.
#'		#results <- nbinomTest(data, "M", "WE")
#'	
#'		# The below two code lines load a copy of the above dataset
#'		# which has already been processed by:
#'		#     results <- nbinomTest(data, "M", "WE")
#'		# For your own real data, you must use:
#'		#     results <- nbinomTest(data, "M", "WE")'
#'		# Instead of the two lines below:
#'		data(yeast_nbinomTest_results, package="ELBOW")
#'		results <- yeast_nbinomTest_results	
#'
#'		# obtain the elbow limit for the dataset
#'		# the final step in the analysis pipeline
#'		do_elbow_rnaseq(results)
#'		
##########
do_elbow_rnaseq <- function (rnaSeq) {
	rsvals <- rnaSeq[,"log2FoldChange"]
	is.na(rsvals) <- sapply(rsvals, is.infinite)
	rsvals <- na.omit(rsvals)
	return(do_elbow(data.frame(sort(rsvals))))
}


##########
#' A debug function for plotting data.  This function plots an array
#' or R object, and optionally draws the ELBOW limits on the plot.
#' This function is mainly useful for the development and maintenance
#' of the ELBOW method R package.
#'
#' @export
#'
#' @title plots data to see if it matches an elbow distribution.
#'
#' @param my_data is the object, or array, to plot.
#'
#' @param column the column in the object/array to plot, if applicable.
#'
#' @param upper_limit the upper ELBOW limit for the data (optional)
#'
#' @param lower_limit the lower ELBOW limit for the data (optional)
#'
#' @examples
#'		# install the DESeq libraries
#'		#source("http://www.bioconductor.org/biocLite.R")
#'		#biocLite("DESeq")
#'
#'		## download the table
#'		library("DESeq")
#'
#'		# the following bam file dataset was obtained from:
#'		# http://cgrlucb.wikispaces.com/file/view/yeast_sample_data.txt
#'		# it has been downloaded into this package for speed convenience.
#'		filename <- system.file("extdata", "yeast_sample_data.txt", package = "ELBOW")
#'	
#'		count_table <- read.table(filename, header=TRUE, sep="\t", row.names=1)
#'		expt_design <- data.frame(row.names = colnames(count_table), condition = c("WE","WE","M","M","M"))
#'		conditions = expt_design$condition
#'		data <- newCountDataSet(count_table, conditions)
#'		data <- estimateSizeFactors(data)
#'		data <- as(data, "CountDataSet")
#'		## data <- estimateVarianceFunctions(data)
#'		data <- estimateDispersions(data)
#'
#'
#'		# this next step is essential, but it takes a long time...
#'		# so, just like a good cooking show we will skip this step
#'		# and load a finished version.
#'		#results <- nbinomTest(data, "M", "WE")
#'	
#'		# The below two code lines load a copy of the above dataset
#'		# which has already been processed by:
#'		#     results <- nbinomTest(data, "M", "WE")
#'		# For your own real data, you must use:
#'		#     results <- nbinomTest(data, "M", "WE")'
#'		# Instead of the two lines below:
#'		data(yeast_nbinomTest_results, package="ELBOW")
#'		results <- yeast_nbinomTest_results	
#'		
#'		# plot the results w/o elbow (useful if do_elbow doesn't work)
#'		#plot_dataset(results, "log2FoldChange")
#'
#'		# plot the results w/elbow
#'		limits <- do_elbow_rnaseq(results)
#'		plot_dataset(results, "log2FoldChange", limits$up_limit, limits$low_limit)
#'
##########
plot_dataset <- function (my_data, column = NULL, upper_limit = NULL, lower_limit = NULL) {
	if (!is.null(column)) {
		fold_values <- my_data[,column]
	} else {
		fold_values <- my_data
	}
	fold_values <- sort(fold_values)
	plot(fold_values, type='l', lwd=2, col='red', main = "Elbow Plot Results") 
	if (!is.null(upper_limit)) {
		abline(h = lower_limit, untf = TRUE)
		abline(h = upper_limit, untf = TRUE)
	}
}

source("R/insilico_core.r")
source("R/diag.r")
#* Implement InSilicoVA methods. This function implements InSilicoVA model. The InSilicoVA model is fitted with MCMC implemented in Java. For more detail, see the paper on \url{http://arxiv.org/abs/1411.3042}.
#* @param data The original data to be used. It is suggested to use similar input as InterVA4, with the first column being death IDs and 245 symptoms.  The only difference in input is InsilicoVA takes three levels: ``present'', ``absent'', and ``missing (no data)''. Similar to InterVA software, ``present'' symptoms takes value ``Y''; ``absent'' symptoms take take value ``NA'' or ``''. For missing symptoms, e.g., questions not asked or answered in the original interview, corrupted data, etc., the input should be coded by ``.'' to distinguish from ``absent'' category. The order of the columns does not matter as long as the column names are correct. It can also include more  unused columns than the standard InterVA4 input. But the first column should be  the death ID. For example input data format, see \code{RandomVA1} and  \code{RandomVA2}.
#* @param data.type Type of questionnaire. ``WHO2012'' corresponds to the standard input of InterVA4, and  ``WHO2016'' corresponds to the standard input of InterVA5.
#* @param sci A data frame that contains the symptom-cause-information (aka Probbase) that InterVA uses to assign a cause of death.
#* @param isNumeric Indicator if the input is already in numeric form. If the input is coded numerically such that 1 for ``present'', 0 for ``absent'', and -1 for ``missing'', this indicator could be set to True to avoid conversion to standard InterVA format.
#* @param updateCondProb Logical indicator. If FALSE, then fit InSilicoVA model without re-estimating conditional probabilities.
#* @param keepProbbase.level Logical indicator when \code{updateCondProb} is FALSE. If TRUE, then only estimate the InterVA's conditional probability interpretation table; if FALSE, estimate the whole conditional probability matrix. Default to TRUE.
#* @param CondProb Customized conditional probability matrix to use.It should be strict the same configuration as InterVA-4 software. That is, it should be a matrix of 245 rows of symptoms and 60 columns of causes, arranged in the same order as in InterVA-4 specification. The elements in the matrix should be the conditional probability of corresponding symptom given the corresponding cause, represented in alphabetic form indicating levels. For example input, see \code{\link{condprob}}
#* @param CondProbNum Customized conditional probability matrix to use if specified fully by numerical values between 0 and 1. If it is specified, re-estimation of conditional probabilities will not be performed, i.e., \code{updateCondProb} will be set to FALSE.
#* @param datacheck Logical indicator for whether to check the data satisfying InterVA rules. Default set to be TRUE. If \code{warning.write} is set to true, the inconsistent input will be logged in file warning_insilico.txt and errorlog_insilico.txt. It's strongly suggested to be set to TRUE.
#* @param datacheck.missing Logical indicator for whether to perform data check before deleting complete missing symptoms. Default to TRUE.
#* @param warning.write Logical indicator for whether to save the changes made to data input by \code{datacheck}. If set to TRUE, the changes will be logged in file warning_insilico.txt and errorlog_insilico.txt in current working directory.
#* @param directory The directory to store the output from. It should be an valid existing directory or a folder to be created.
#* @param external.sep Logical indicator for whether to separate out external causes first. Default set to be TRUE. If set to TRUE, the algorithm will estimate external causes, e.g., traffic accident, accidental fall, suicide, etc., by checking the corresponding indicator only without considering other medical symptoms. It is strongly suggested to set to be TRUE.
#* @param Nsim Number of iterations to run. Default to be 4000.
#* @param thin Proportion of thinning for storing parameters. For example, if thin = k, the output parameters will only be saved every k iterations. Default to be 10
#* @param burnin Number of iterations as burn-in period. Parameters sampled in burn-in period will not be saved.
#* @param auto.length Logical indicator of whether to automatically increase chain length if convergence not reached.
#* @param conv.csmf Minimum CSMF value to check for convergence if auto.length is set to TRUE. For example, under the default value 0.02, all causes with mean CSMF at least 0.02 will be checked for convergence.
#* @param jump.scale The scale of Metropolis proposal in the Normal model. Default to be 0.1.
#* @param levels.prior Vector of prior expectation of conditional probability levels. They do not have to be scaled. The algorithm internally calibrate the scale to the working scale through \code{levels.strength}. If NULL the algorithm will use InterVA table as prior.
#* @param levels.strength Scaling factor for the strength of prior beliefs in the conditional probability levels. Larger value constrain the posterior estimates to be closer to prior expectation. Defult value 1 scales \code{levels.prior} to a suggested scale that works empirically.
#* @param trunc.min Minimum possible value for estimated conditional probability table. Default to be 0.0001
#* @param trunc.max Maximum possible value for estimated conditional probability table. Default to be 0.9999
#* @param subpop This could be the column name of the variable in data that is to be used as sub-population indicator, or a list of column names if more than one  variable are to be used. Or it could be a vector of sub-population assignments  of the same length of death records. It could be numerical indicators or character  vectors of names.
#* @param java_option Option to initialize java JVM. Default to ``-Xmx1g'', which sets the maximum heap size to be 1GB. If R produces ``java.lang.OutOfMemoryError: Java heap space'' error message, consider increasing heap size using this option, or one of the following: (1) decreasing \code{Nsim}, (2) increasing \code{thin}, or (3) disabling \code{auto.length}.
#* @param seed Seed used for initializing sampler. The algorithm will produce the same outcome with the same seed in each machine.
#* @param phy.code A matrix of physician assigned cause distribution. The physician assigned causes need not be the same as the list of causes used in InSilicoVA and InterVA-4. The cause list used could be a higher level aggregation of the InSilicoVA causes. See \code{phy.cat} for more detail. The first column of \code{phy.code} should be death ID that could be matched to the symptom dataset, the following columns are the probabilities of each cause category used by physicians.
#* @param phy.cat A two column matrix describing the correspondence between InSilicoVA causes and the physician assigned causes. Note each InSilicoVA cause (see \code{causetext}) could only correspond to one physician assigned cause. See \code{SampleCategory} for an example. 'Unknown' category should not be included in this matrix.
#* @param phy.unknown The name of the physician assigned cause that correspond to unknown COD.
#* @param phy.external The name of the physician assigned cause that correspond to external causes. This will only be used if \code{external.sep} is set to TRUE. In that case, all external causes should be grouped together, as they are assigned deterministically by the corresponding symptoms.
#* @param phy.debias Fitted object from physician coding debias function (see \code{\link{physician_debias}}) that overwrites \code{phy.code}.
#* @param exclude.impossible.cause option to exclude impossible causes at the individual level. The following rules are implemented: \code{subset}: Causes with 0 probability given the age group and gender of the observation, according to the InterVA conditional probabilities, are removed; \code{subset2}: In addition to the same rules as \code{subset}, also remove Prematurity for baby born during at least 37 weeks of pregnancy and remove Birth asphyxia for baby  not born during at least 37 weeks of pregnancy; \code{all}: Causes with 0 probability given any symptom of the observation, according to the InterVA conditional probabilities, are removed; \code{interVA}:  Causes with 0 probability given any positive indicators according to the InterVA conditional probabilities, are removed; and \code{none}: no causes are removed. \code{subset2} is the default.
#* @param no.is.missing logical indicator to treat all absence of symptoms as missing. Default to FALSE. If set to TRUE, InSilicoVA will perform calculations similar to InterVA-4 w.r.t treating absent symptoms. It is highly recommended to set this argument to FALSE.
#* @param indiv.CI credible interval for individual probabilities. If set to NULL, individual COD distributions will not be calculated to accelerate model fitting time. See \code{\link{get.indiv}} for details of updating the C.I. later after fitting the model.
#* @param groupcode logical indicator of including the group code in the output causes
#* @post /insilicova


insilico <- function(data, data.type = c("WHO2012", "WHO2016")[1], sci = NULL, isNumeric = FALSE,
  updateCondProb = TRUE, keepProbbase.level = TRUE,
  CondProb = NULL, CondProbNum = NULL, datacheck = TRUE, datacheck.missing = TRUE,
  warning.write = FALSE, directory = NULL, external.sep = TRUE, Nsim = 4000, thin = 10, burnin = 2000,
  auto.length = TRUE, conv.csmf = 0.02, jump.scale = 0.1,
  levels.prior = NULL, levels.strength = 1, trunc.min = 0.0001, trunc.max = 0.9999,
  subpop = NULL, java_option = "-Xmx1g", seed = 1,
  phy.code = NULL, phy.cat = NULL, phy.unknown = NULL, phy.external = NULL,
  phy.debias = NULL, exclude.impossible.cause = c("subset2", "subset", "all", "InterVA", "none")[1],
  no.is.missing = FALSE, indiv.CI = NULL, groupcode=FALSE, ...){

	# handling changes throughout time
	  args <- as.list(match.call())
	  if(!is.null(args$length.sim)){
	  	Nsim <- args$length.sim
	  	message("length.sim argument is replaced with Nsim argument, will remove in later versions.\n")
	  }

	  if(methods::is(exclude.impossible.cause, "logical")){
	  	stop(paste0("exclude.impossible.cause now has more options. Please check out the package documentation. Rest to '", exclude.impossible.cause, "'."))
	  }

	fit <- insilico.fit(data = data,
						data.type = data.type,
						sci = sci,
						isNumeric = isNumeric,
						updateCondProb = updateCondProb,
						keepProbbase.level = keepProbbase.level,
						CondProb = CondProb,
						CondProbNum = CondProbNum,
						datacheck = datacheck,
						datacheck.missing = datacheck.missing,
						warning.write = warning.write,
						directory = directory,
						external.sep = external.sep,
						Nsim = Nsim,
						thin = thin,
						burnin = burnin,
						auto.length = auto.length,
						conv.csmf = conv.csmf,
						jump.scale = jump.scale,
						levels.prior = levels.prior,
						levels.strength = levels.strength,
						trunc.min = trunc.min,
						trunc.max = trunc.max,
						subpop = subpop,
						java_option = java_option,
						seed = seed,
						phy.code = phy.code,
						phy.cat = phy.cat,
						phy.unknown = phy.unknown,
						phy.external = phy.external,
						phy.debias = phy.debias,
						exclude.impossible.cause = exclude.impossible.cause,
						no.is.missing = no.is.missing,
						indiv.CI = indiv.CI,
						groupcode=groupcode)
	return(fit)
}

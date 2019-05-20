## insilicosystemsargs class ----

#' Constructor function for the \code{insilicosystemsargs} class.
#'
#' Constructor function for the \code{insilicosystemsargs} class, with default values for the parameters if not provided by the user.
#'
#' For the protein-coding (and non-coding) biological function ratios (i.e. PC.TC.p, PC.TL.p, etc): if none of the ratios are provided,
#' then they are set to their default values. Otherwise, if at least one value among the 6 (5 for noncoding genes) is set by the user:
#' \itemize{
#' \item if the sum of the provided values is 1 or more: the non-specified values are set to 0, and the specified values are normalised
#' such that their sum is 1.
#' \item if the sum of the provided values is less than 1: the non-specified values are set such that the sum of all ratios equals 1.
#' }
#' Example: if the user sets PC.TC.p to 1 and PC.TL.p to 0.6, but does not provide any values for the other ratios,
#' then PC.TC.p is set to 1/(1+0.6)=0.625, PC.TL.p to 0.6/(1+0.6)=0.375, and PC.RD.p, PC.PD.p, PC.PTM.p and PC.MR.p are all set to 0.
#' Accordingly, if the user only sets NC.TC.p to 0.6, then NC.TL.p, NC.RD.p, NC.PD.p and NC.PTM.p are all set to 0.1.
#'
#' @param G Integer. Number of genes in the system. Default value is 10.
#' @param PC.p Numeric. Probability of each gene to be a protein-coding gene. Default value is 0.7.
#' @param PC.TC.p Numeric. Probability of a protein-coding gene to be a regulator of transcription. Default value is 0.4 (see details).
#' @param PC.TL.p Numeric. Probability of a protein-coding gene to be a regulator of translation. Default value is 0.3 (see details).
#' @param PC.RD.p Numeric. Probability of a protein-coding gene to be a regulator of RNA decay. Default value is 0.1 (see details).
#' @param PC.PD.p Numeric. Probability of a protein-coding gene to be a regulator of protein decay. Default value is 0.1 (see details).
#' @param PC.PTM.p Numeric. Probability of a protein-coding gene to be a regulator of protein post-translational modification. Default value is 0.05 (see details).
#' @param PC.MR.p Numeric. Probability of a protein-coding gene to be a metabolic enzyme. Default value is 0.05 (see details).
#' @param NC.TC.p Numeric. Probability of a noncoding gene to be a regulator of transcription. Default value is 0.3 (see details).
#' @param NC.TL.p Numeric. Probability of a noncoding gene to be a regulator of translation. Default value is 0.3 (see details).
#' @param NC.RD.p Numeric. Probability of a noncoding gene to be a regulator of RNA decay. Default value is 0.3 (see details).
#' @param NC.PD.p Numeric. Probability of a noncoding gene to be a regulator of protein decay. Default value is 0.05 (see details).
#' @param NC.PTM.p Numeric. Probability of a noncoding gene to be a regulator of protein post-translational modification. Default value is 0.05 (see details).
#' @param TC.pos.p Numeric. Probability of a regulation targeting gene transcription to be positive. Default value is 0.5.
#' @param TL.pos.p Numeric. Probability of a regulation targeting gene translation to be positive. Default value is 0.5.
#' @param PTM.pos.p Numeric. Probability of a regulation targeting protein post-translational modification to be positive (i.e the targeted protein
#' is transformed into its modified form, as opposed to the modified protein being transformed back into its original form). Default value is 0.5.
#' @param basal_transcription_rate_samplingfct Function from which the transcription rates of genes are sampled (input x is the required sample size). Default value is
#' a function returning \code{(10^v)/3600}, with \code{v} a vector of size x sampled from a normal distribution with mean of 3 and sd of 0.5.
#' @param basal_translation_rate_samplingfct Function from which the translation rates of genes are sampled (input x is the required sample size). Default value is
#' a function returning \code{(10^v)/3600}, with \code{v} a vector of size x sampled from a normal distribution with mean of 2.146 and sd of 0.7.
#' @param basal_RNAlifetime_samplingfct Function from which the transcript lifetimes are sampled (input x is the required sample size). Default value is
#' a function returning \code{(10^v)*3600}, with \code{v} a vector of size x sampled from a normal distribution with mean of 0.95 and sd of 0.2.
#' @param basal_protlifetime_samplingfct Function from which the protein lifetime are sampled (input x is the required sample size). Default value is
#' a function returning \code{(10^v)*3600}, with \code{v} a vector of size x sampled from a normal distribution with mean of 1.3 and sd of 0.4.
#' @param TC.PC.outdeg.distr Form of the distribution of the number of targets (out-degree) of protein regulators in the transcription regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TC.NC.outdeg.distr Form of the distribution of the number of targets (out-degree) of noncoding regulators in the transcription regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TC.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the transcription regulation graph. Default value is 3.
#' @param TC.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the transcription regulation graph. Default value is 5.
#' @param TC.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the transcription regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TC.NC.indeg.distr Type of preferential attachment for the targets of noncoding regulators in the transcription regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TC.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation in the transcription regulation graph. Default value is 0.2.
#' @param TC.NC.autoregproba Numeric. Probability of noncoding regulators to perform autoregulation in the transcription regulation graph. Default value is 0.
#' @param TC.PC.twonodesloop Logical. Are 2-nodes loops authorised in the transcription regulation graph with protein regulators? Default value is FALSE.
#' @param TC.NC.twonodesloop Logical. Are 2-nodes loops authorised in the transcription regulation graph with noncoding regulators? Default value is FALSE.
#' @param TCbindingrate_samplingfct Function from which the binding rates of transcription regulators on their targets are sampled (input \code{means} is a vector of length equal to the required sample size, giving for each
#' edge (regulatory interaction) for which a binding rate is being sampled the value of the sampled unbinding rate divided by the steady-state abundance of the regulator
#' in absence of any regulation in the system). Default value is a function returning \code{10^v}, where \code{v} is a vector with the same length as \code{means} whose elements are sampled
#' from a truncated normal distribution with mean equal to the log10 of the corresponding element in \code{means}, and sd = 0.1, the minimum authorised value being the log10 of the corresponding element in \code{means}.
#' @param TCunbindingrate_samplingfct Function from which the unbinding rates of transcription regulators from their target are sampled (input x is the required sample size). Default value is
#' a function returning \code{10^v}, with \code{v} a vector of size x sampled from a normal distribution with mean of -3 and sd of 0.2.
#' @param TCfoldchange_samplingfct Function from which the transcription fold change induced by a bound regulator is sampled (input x is the required sample size). Default value is
#' a truncated normal distribution with a mean of 3, sd of 10 and minimum authorised value of 1.5.
#' @param TL.PC.outdeg.distr Form of the distribution of the number of targets (out-degree) of protein regulators in the translation regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TL.NC.outdeg.distr Form of the the distribution of the number of targets (out-degree) of noncoding regulators in the translation regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TL.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the translation regulation graph. Default value is 4.
#' @param TL.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the translation regulation graph. Default value is 6.
#' @param TL.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the translation regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TL.NC.indeg.distr Type of preferential attachment for the targets of noncoding regulators in the translation regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TL.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation in the translation regulation graph. Default value is 0.2.
#' @param TL.NC.autoregproba Numeric. Probability of noncoding regulators to perform autoregulation in the translation regulation graph. Default value is 0.
#' @param TL.PC.twonodesloop Logical. Are 2-nodes loops authorised in the translation regulation graph with protein regulators? Default value is FALSE.
#' @param TL.NC.twonodesloop Logical. Are 2-nodes loops authorised in the translation regulation graph with noncoding regulators? Default value is FALSE.
#' @param TLbindingrate_samplingfct Function from which the binding rate of translation regulators on target are sampled (input \code{means} is a vector of length equal to the required sample size, giving for each
#' edge (regulatory interaction) for which a binding rate is being sampled the value of the sampled unbinding rate divided by the steady-state abundance of the regulator
#' in absence of any regulation in the system). Default value is a function returning \code{10^v}, where \code{v} is a vector with the same length as \code{means} whose elements are sampled
#' from a truncated normal distribution with mean equal to the log10 of the corresponding element in \code{means}, and sd = 0.1, the minimum authorised value being the log10 of the corresponding element in \code{means}.
#' @param TLunbindingrate_samplingfct Function from which the unbinding rate of translation regulators from target are sampled (input x is the required sample size). Default value is
#' a function returning \code{10^v}, with \code{v} a vector of size x sampled from a normal distribution with mean of -3 and sd of 0.2.
#' @param TLfoldchange_samplingfct Function from which the translation fold change induced by a bound regulator are sampled (input x is the required sample size). Default value is
#' a truncated normal distribution with a mean of 3, sd of 10 and minimum authorised value of 1.5.
#' @param RD.PC.outdeg.distr Form of the distribution of the number of targets (out-degree) of protein regulators in the RNA decay regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param RD.NC.outdeg.distr Form of the the distribution of the number of targets (out-degree) of noncoding regulators in the RNA decay regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param RD.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the RNA decay regulation graph. Default value is 4.
#' @param RD.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the RNA decay regulation graph. Default value is 6.
#' @param RD.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the RNA decay graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param RD.NC.indeg.distr Type of preferential attachment for the targets of noncoding regulators in the RNA decay graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param RD.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation in the RNA decay regulation graph. Default value is 0.2.
#' @param RD.NC.autoregproba Numeric. Probability of noncoding regulators to perform autoregulation in the RNA decay regulation graph. Default value is 0.
#' @param RD.PC.twonodesloop Logical. Are 2-nodes loops authorised in the RNA decay regulation graph with protein regulators? Default value is FALSE.
#' @param RD.NC.twonodesloop Logical. Are 2-nodes loops authorised in the RNA decay regulation graph with noncoding regulators? Default value is FALSE.
#' @param RDregrate_samplingfct Function from which the RNA decay rates of targets of RNA decay regulators are sampled (input x is the required sample size). Default value is
#' a function returning \code{10^v}, with \code{v} a vector of size x sampled from a normal distribution with mean of -5 and sd of 1.5.
#' @param PD.PC.outdeg.distr Form of the distribution of the number of targets (out-degree) of protein regulators in the protein decay regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PD.NC.outdeg.distr Form of the the distribution of the number of targets (out-degree) of noncoding regulators in the protein decay regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PD.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the protein decay regulation graph. Default value is 4.
#' @param PD.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the protein decay regulation graph. Default value is 6.
#' @param PD.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the protein decay regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PD.NC.indeg.distr Type of preferential attachment for the targets of noncoding regulators in the protein decay graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PD.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation in the protein decay regulation graph. Default value is 0.2.
#' @param PD.NC.autoregproba Numeric. Probability of noncoding regulators to perform autoregulation in the protein decay regulation graph. Default value is 0.
#' @param PD.PC.twonodesloop Logical. Are 2-nodes loops authorised in the protein decay graph with protein regulators in the protein decay regulation graph? Default value is FALSE.
#' @param PD.NC.twonodesloop Logical. Are 2-nodes loops authorised in the protein decay graph with noncoding regulators in the protein decay regulation graph? Default value is FALSE.
#' @param PDregrate_samplingfct Function from which the protein decay rates of targets of protein decay regulators are sampled (input x is the required sample size).  Default value is
#' a function returning \code{10^v}, with \code{v} a vector of size x sampled from a normal distribution with mean of -5 and sd of 1.5.
#' @param PTM.PC.outdeg.distr Form of the distribution of the number of targets (out-degree) of protein regulators in the post-translational modification regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PTM.NC.outdeg.distr Form of the the distribution of the number of targets (out-degree) of noncoding regulators in the post-translational modification regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PTM.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the protein post-translational modification graph. Default value is 4.
#' @param PTM.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the protein post-translational modification graph. Default value is 6.
#' @param PTM.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the protein post-translational modification graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PTM.NC.indeg.distr Type of preferential attachment for the targets of noncoding regulators in the protein post-translational modification graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PTM.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation. Default value is 0.2.
#' @param PTM.NC.autoregproba Numeric. Probability of noncoding regulators to perform autoregulation. Default value is 0.
#' @param PTM.PC.twonodesloop Logical. Are 2-nodes loops authorised in the protein post-translational modification graph with protein regulators? Default value is FALSE.
#' @param PTM.NC.twonodesloop Logical. Are 2-nodes loops authorised in the protein post-translational modification graph with noncoding regulators? Default value is FALSE.
#' @param PTMregrate_samplingfct Function from which the protein transformation rates of targets of post-translational modification regulators are sampled (input x is the required sample size). Default value is
#' a function returning \code{10^v}, with \code{v} a vector of size x sampled from a normal distribution with mean of -5 and sd of 1.5.
#' @param regcomplexes Can the regulators controlling a common target form regulatory complexes in the different regulatory graphs? Can be 'none', 'prot' (only protein can form regulatory complexes) or 'both'
#' (both regulatory RNAs and proteins can form regulatory complexes). Default value is "prot".
#' @param regcomplexes.p Numeric. Probability that regulators controlling a common target form regulatory complexes; ignore if \code{regcomplexes} = 'none'. Default value is 0.3.
#' @param regcomplexes.size Integer. Number of components of a regulatory complex; ignore if \code{regcomplexes} = 'none'. Default value is 2.
#' @param complexesformationrate_samplingfct Function from which the formation rate of regulatory complexes are sampled (input x is the required sample size). Default value is
#' a function returning \code{10^v}, with \code{v} a vector of size x sampled from a normal distribution with mean of -3 and sd of 0.7.
#' @param complexesdissociationrate_samplingfct Function from which the dissociation rate of regulatory complexes are sampled (input x is the required sample size). Default value is
#' a function returning \code{10^v}, with \code{v} a vector of size x sampled from a normal distribution with mean of 3 and sd of 0.9.
#' @return An object of the class \code{insilicosystemargs}, that is a named list of the different parameters.
#' @examples
#' sysargs = insilicosystemargs(G = 15, PC.p = 0.2,
#'  basal_transcription_rate_samplingfct = function(x){runif(x, 0.1, 0.2)})
#' @export
insilicosystemargs <- function(
  G = 10,
  PC.p = 0.7,
  PC.TC.p = NULL,
  PC.TL.p = NULL,
  PC.RD.p = NULL,
  PC.PD.p = NULL,
  PC.PTM.p = NULL,
  PC.MR.p = NULL,
  NC.TC.p = NULL,
  NC.TL.p = NULL,
  NC.RD.p = NULL,
  NC.PD.p = NULL,
  NC.PTM.p = NULL,
  TC.pos.p = 0.5,
  TL.pos.p = 0.5,
  PTM.pos.p = 0.5,
  #  TC.PC.pos.p = 0.5,
  #  TC.NC.pos.p = 0.5,
  #  TL.PC.pos.p = 0.5,
  #  TL.NC.pos.p = 0.5,
  #  PTM.PC.pos.p = 0.5,
  #  PTM.NC.pos.p = 0.5,
  basal_transcription_rate_samplingfct = NULL,
  basal_translation_rate_samplingfct = NULL,
  basal_RNAlifetime_samplingfct = NULL,
  basal_protlifetime_samplingfct = NULL,
  TC.PC.outdeg.distr = "powerlaw",
  TC.NC.outdeg.distr = "powerlaw",
  TC.PC.outdeg.exp = 3,
  TC.NC.outdeg.exp = 5,
  TC.PC.indeg.distr = "powerlaw",
  TC.NC.indeg.distr = "powerlaw",
  TC.PC.autoregproba = 0.2,
  TC.NC.autoregproba = 0,
  TC.PC.twonodesloop = FALSE,
  TC.NC.twonodesloop = FALSE,
  TCbindingrate_samplingfct = NULL,
  TCunbindingrate_samplingfct = NULL,
  TCfoldchange_samplingfct = NULL,
  TL.PC.outdeg.distr = "powerlaw",
  TL.NC.outdeg.distr = "powerlaw",
  TL.PC.outdeg.exp = 4,
  TL.NC.outdeg.exp = 6,
  TL.PC.indeg.distr = "powerlaw",
  TL.NC.indeg.distr = "powerlaw",
  TL.PC.autoregproba = 0.2,
  TL.NC.autoregproba = 0,
  TL.PC.twonodesloop = FALSE,
  TL.NC.twonodesloop = FALSE,
  TLbindingrate_samplingfct = NULL,
  TLunbindingrate_samplingfct = NULL,
  TLfoldchange_samplingfct = NULL,
  RD.PC.outdeg.distr = "powerlaw",
  RD.NC.outdeg.distr = "powerlaw",
  RD.PC.outdeg.exp = 4,
  RD.NC.outdeg.exp = 6,
  RD.PC.indeg.distr = "powerlaw",
  RD.NC.indeg.distr = "powerlaw",
  RD.PC.autoregproba = 0.2,
  RD.NC.autoregproba = 0,
  RD.PC.twonodesloop = FALSE,
  RD.NC.twonodesloop = FALSE,
  RDregrate_samplingfct = NULL,
  PD.PC.outdeg.distr = "powerlaw",
  PD.NC.outdeg.distr = "powerlaw",
  PD.PC.outdeg.exp = 4,
  PD.NC.outdeg.exp = 6,
  PD.PC.indeg.distr = "powerlaw",
  PD.NC.indeg.distr = "powerlaw",
  PD.PC.autoregproba = 0.2,
  PD.NC.autoregproba = 0,
  PD.PC.twonodesloop = FALSE,
  PD.NC.twonodesloop = FALSE,
  PDregrate_samplingfct = NULL,
  PTM.PC.outdeg.distr = "powerlaw",
  PTM.NC.outdeg.distr = "powerlaw",
  PTM.PC.outdeg.exp = 4,
  PTM.NC.outdeg.exp = 6,
  PTM.PC.indeg.distr = "powerlaw",
  PTM.NC.indeg.distr = "powerlaw",
  PTM.PC.autoregproba = 0.2,
  PTM.NC.autoregproba = 0,
  PTM.PC.twonodesloop = FALSE,
  PTM.NC.twonodesloop = FALSE,
  PTMregrate_samplingfct = NULL,
  regcomplexes = "prot",
  regcomplexes.p = 0.3,
  regcomplexes.size = 2,
  complexesformationrate_samplingfct = NULL,
  complexesdissociationrate_samplingfct = NULL
){

  if(is.null(basal_transcription_rate_samplingfct)) basal_transcription_rate_samplingfct = function(x){ logval = rnorm(x, mean = -0.92, sd = 0.35); val = 10^logval; return(val/60) }
  if(is.null(basal_translation_rate_samplingfct)) basal_translation_rate_samplingfct = function(x){ logval = rnorm(x, mean = 2.146, sd = 0.7); val = 10^logval; return(val/3600) }
  if(is.null(basal_RNAlifetime_samplingfct)) basal_RNAlifetime_samplingfct = function(x){ logval = rnorm(x, mean = 1.36, sd = 0.2); val = 10^logval; return(val*60) }
  if(is.null(basal_protlifetime_samplingfct)) basal_protlifetime_samplingfct = function(x){ logval = rnorm(x, mean = 5.43, sd = 1); val = 2^logval; return(val*60) }
  if(is.null(TCbindingrate_samplingfct)) TCbindingrate_samplingfct = function(means){ logval = truncnorm::rtruncnorm(length(means),a = log10(means), mean = log10(means), sd = 0.1); return(10^logval) }
  if(is.null(TCunbindingrate_samplingfct)) TCunbindingrate_samplingfct = function(x){ logval = rnorm(x, mean = -3, sd = 0.2); return(10^(logval)) }
  if(is.null(TCfoldchange_samplingfct)) TCfoldchange_samplingfct = function(x){ truncnorm::rtruncnorm(x, a = 1.5, mean = 3, sd = 10) }
  if(is.null(TLbindingrate_samplingfct)) TLbindingrate_samplingfct = function(means){ logval = truncnorm::rtruncnorm(length(means),a = log10(means), mean = log10(means), sd = 0.1); return(10^logval)  }
  if(is.null(TLunbindingrate_samplingfct)) TLunbindingrate_samplingfct = function(x){ logval = rnorm(x, mean = -3, sd = 0.2); return(10^(logval)) }
  if(is.null(TLfoldchange_samplingfct)) TLfoldchange_samplingfct = function(x){ truncnorm::rtruncnorm(x, a = 1.5, mean = 3, sd = 10)  }
  if(is.null(RDregrate_samplingfct)) RDregrate_samplingfct = function(x){ logval = rnorm(x, mean = -5, sd = 1.5); return(10^logval) }
  if(is.null(PDregrate_samplingfct)) PDregrate_samplingfct = function(x){ logval = rnorm(x, mean = -5, sd = 1.5); return(10^logval) }
  if(is.null(PTMregrate_samplingfct)) PTMregrate_samplingfct = function(x){ logval = rnorm(x, mean = -8, sd = 1.5); return(10^logval) }
  if(is.null(complexesformationrate_samplingfct)) complexesformationrate_samplingfct = function(x){ logval = rnorm(x, mean = -3, sd = 0.7); return(10^(logval)) }
  if(is.null(complexesdissociationrate_samplingfct)) complexesdissociationrate_samplingfct = function(x){ logval = rnorm(x, mean = 2, sd = 1.2); return(10^(logval)) }

  NC.p = 1 - PC.p

  ## Probability of each protein-coding gene biological function
  temp = c(PC.TC.p, PC.TL.p, PC.RD.p, PC.PD.p, PC.PTM.p, PC.MR.p)

  if(is.null(temp)){ ## if no values are provided, give default values
    PC.TC.p = 0.4
    PC.TL.p = 0.3
    PC.RD.p = 0.1
    PC.PD.p = 0.1
    PC.PTM.p = 0.05
    PC.MR.p = 0.05
  } else if(sum(temp) >= 1 | length(temp) == 6){ ## if at least one value is provided, and their sum is >=1, normalise the given values and set to 0 the others
    for(v in c("PC.TC.p", "PC.TL.p", "PC.RD.p", "PC.PD.p", "PC.PTM.p", "PC.MR.p")){
      assign(v, ifelse(is.null(get(v)), 0, get(v)/sum(temp)))
    }
  } else if(sum(temp) < 1){ ## if at least one value is provided but don't sum up to 1, assign the remaining proba such that the sum of all proba is 1
    residual = (1-sum(temp))/(6-length(temp))
    for(v in c("PC.TC.p", "PC.TL.p", "PC.RD.p", "PC.PD.p", "PC.PTM.p", "PC.MR.p")){
      assign(v, ifelse(is.null(get(v)), residual, get(v)))
    }
  }

  ## Probability of each noncoding gene biological function
  temp = c(NC.TC.p, NC.TL.p, NC.RD.p, NC.PD.p, NC.PTM.p)

  if(is.null(temp)){ ## if no values are provided, give default values
    NC.TC.p = 0.3
    NC.TL.p = 0.3
    NC.RD.p = 0.3
    NC.PD.p = 0.05
    NC.PTM.p = 0.05
  } else if(sum(temp) >= 1 | length(temp) == 5){ ## if at least one value is provided, and their sum is >=1, normalise the given values and set to 0 the others
    for(v in c("NC.TC.p", "NC.TL.p", "NC.RD.p", "NC.PD.p", "NC.PTM.p")){
      assign(v, ifelse(is.null(get(v)), 0, get(v)/sum(temp)))
    }
  } else if(sum(temp) < 1){ ## if at least one value is provided but don't sum up to 1, assign the remaining proba such that the sum of all proba is 1
    residual = (1-sum(temp))/(6-length(temp))
    for(v in c("NC.TC.p", "NC.TL.p", "NC.RD.p", "NC.PD.p", "NC.PTM.p")){
      assign(v, ifelse(is.null(get(v)), residual, get(v)))
    }
  }




  ## There cannot be negative regulation for transcript and protein decay
  RD.pos.p = 1 ## RD.PC.pos.p: probability that the RNA decay is positively regulated by protein regulators (faster decay)
  PD.pos.p = 1 ## PD.PC.pos.p: probability that the protein decay is positively regulated by protein regulators (faster decay)
  #  RD.PC.pos.p = 1 ## RD.PC.pos.p: probability that the RNA decay is positively regulated by protein regulators (faster decay)
  #  RD.NC.pos.p = 1 ## RD.NC.pos.p: probability that the RNA decay is positively regulated by noncoding regulators (faster decay)
  #  PD.PC.pos.p = 1 ## PD.PC.pos.p: probability that the protein decay is positively regulated by protein regulators (faster decay)
  #  PD.NC.pos.p = 1 ## PD.NC.pos.p: probability that the protein decay is positively regulated by noncoding regulators (faster decay)

  value = list(  "G" = G,
                 "PC.p" = PC.p,
                 "PC.TC.p" = PC.TC.p,
                 "PC.TL.p" = PC.TL.p,
                 "PC.RD.p" = PC.RD.p,
                 "PC.PD.p" = PC.PD.p,
                 "PC.PTM.p" = PC.PTM.p,
                 "PC.MR.p" = PC.MR.p,
                 "NC.p" = NC.p,
                 "NC.TC.p" = NC.TC.p,
                 "NC.TL.p" = NC.TL.p,
                 "NC.RD.p" = NC.RD.p,
                 "NC.PD.p" = NC.PD.p,
                 "NC.PTM.p" = NC.PTM.p,
                 "TC.pos.p" = TC.pos.p,
                 "TL.pos.p" = TL.pos.p,
                 "RD.pos.p" = RD.pos.p,
                 "PD.pos.p" = PD.pos.p,
                 "PTM.pos.p" = PTM.pos.p,
                 #                 "TC.PC.pos.p" = TC.PC.pos.p,
                 #                 "TC.NC.pos.p" = TC.NC.pos.p,
                 #                 "TL.PC.pos.p" = TL.PC.pos.p,
                 #                 "TL.NC.pos.p" = TL.NC.pos.p,
                 #                 "RD.PC.pos.p" = RD.PC.pos.p,
                 #                 "RD.NC.pos.p" = RD.NC.pos.p,
                 #                 "PD.PC.pos.p" = PD.PC.pos.p,
                 #                 "PD.NC.pos.p" = PD.NC.pos.p,
                 #                 "PTM.PC.pos.p" = PTM.PC.pos.p,
                 #                 "PTM.NC.pos.p" = PTM.NC.pos.p,
                 "basal_transcription_rate_samplingfct" = basal_transcription_rate_samplingfct,
                 "basal_translation_rate_samplingfct" = basal_translation_rate_samplingfct,
                 "basal_RNAlifetime_samplingfct" = basal_RNAlifetime_samplingfct,
                 "basal_protlifetime_samplingfct" = basal_protlifetime_samplingfct,
                 "TC.PC.outdeg.distr" = TC.PC.outdeg.distr,
                 "TC.NC.outdeg.distr" = TC.NC.outdeg.distr,
                 "TC.PC.outdeg.exp" = TC.PC.outdeg.exp,
                 "TC.NC.outdeg.exp" = TC.NC.outdeg.exp,
                 "TC.PC.indeg.distr" = TC.PC.indeg.distr,
                 "TC.NC.indeg.distr" = TC.NC.indeg.distr,
                 "TC.PC.autoregproba" = TC.PC.autoregproba,
                 "TC.NC.autoregproba" = TC.NC.autoregproba,
                 "TC.PC.twonodesloop" = TC.PC.twonodesloop,
                 "TC.NC.twonodesloop" = TC.NC.twonodesloop,
                 "TCbindingrate_samplingfct" = TCbindingrate_samplingfct,
                 "TCunbindingrate_samplingfct" = TCunbindingrate_samplingfct,
                 "TCfoldchange_samplingfct" = TCfoldchange_samplingfct,
                 "TL.PC.outdeg.distr" = TL.PC.outdeg.distr,
                 "TL.NC.outdeg.distr" = TL.NC.outdeg.distr,
                 "TL.PC.outdeg.exp" = TL.PC.outdeg.exp,
                 "TL.NC.outdeg.exp" = TL.NC.outdeg.exp,
                 "TL.PC.indeg.distr" = TL.PC.indeg.distr,
                 "TL.NC.indeg.distr" = TL.NC.indeg.distr,
                 "TL.PC.autoregproba" = TL.PC.autoregproba,
                 "TL.NC.autoregproba" = TL.NC.autoregproba,
                 "TL.PC.twonodesloop" = TL.PC.twonodesloop,
                 "TL.NC.twonodesloop" = TL.NC.twonodesloop,
                 "TLbindingrate_samplingfct" = TLbindingrate_samplingfct,
                 "TLunbindingrate_samplingfct" = TLunbindingrate_samplingfct,
                 "TLfoldchange_samplingfct" = TLfoldchange_samplingfct,
                 "RD.PC.outdeg.distr" = RD.PC.outdeg.distr,
                 "RD.NC.outdeg.distr" = RD.NC.outdeg.distr,
                 "RD.PC.outdeg.exp" = RD.PC.outdeg.exp,
                 "RD.NC.outdeg.exp" = RD.NC.outdeg.exp,
                 "RD.PC.indeg.distr" = RD.PC.indeg.distr,
                 "RD.NC.indeg.distr" = RD.NC.indeg.distr,
                 "RD.PC.autoregproba" = RD.PC.autoregproba,
                 "RD.NC.autoregproba" = RD.NC.autoregproba,
                 "RD.PC.twonodesloop" = RD.PC.twonodesloop,
                 "RD.NC.twonodesloop" = RD.NC.twonodesloop,
                 "RDregrate_samplingfct" = RDregrate_samplingfct,
                 "PD.PC.outdeg.distr" = PD.PC.outdeg.distr,
                 "PD.NC.outdeg.distr" = PD.NC.outdeg.distr,
                 "PD.PC.outdeg.exp" = PD.PC.outdeg.exp,
                 "PD.NC.outdeg.exp" = PD.NC.outdeg.exp,
                 "PD.PC.indeg.distr" = PD.PC.indeg.distr,
                 "PD.NC.indeg.distr" = PD.NC.indeg.distr,
                 "PD.PC.autoregproba" = PD.PC.autoregproba,
                 "PD.NC.autoregproba" = PD.NC.autoregproba,
                 "PD.PC.twonodesloop" = PD.PC.twonodesloop,
                 "PD.NC.twonodesloop" = PD.NC.twonodesloop,
                 "PDregrate_samplingfct" = PDregrate_samplingfct,
                 "PTM.PC.outdeg.distr" = PTM.PC.outdeg.distr,
                 "PTM.NC.outdeg.distr" = PTM.NC.outdeg.distr,
                 "PTM.PC.outdeg.exp" = PTM.PC.outdeg.exp,
                 "PTM.NC.outdeg.exp" = PTM.NC.outdeg.exp,
                 "PTM.PC.indeg.distr" = PTM.PC.indeg.distr,
                 "PTM.NC.indeg.distr" = PTM.NC.indeg.distr,
                 "PTM.PC.autoregproba" = PTM.PC.autoregproba,
                 "PTM.NC.autoregproba" = PTM.NC.autoregproba,
                 "PTM.PC.twonodesloop" = PTM.PC.twonodesloop,
                 "PTM.NC.twonodesloop" = PTM.NC.twonodesloop,
                 "PTMregrate_samplingfct" = PTMregrate_samplingfct,
                 "regcomplexes" = regcomplexes,
                 "regcomplexes.p" = regcomplexes.p,
                 "regcomplexes.size" = regcomplexes.size ,
                 "complexesformationrate_samplingfct" = complexesformationrate_samplingfct,
                 "complexesdissociationrate_samplingfct" = complexesdissociationrate_samplingfct)

  attr(value, "class") = "insilicosystemargs"
  return(value)
}


## insilicoindividualargs class ----

#' Constructor function for the \code{insilicoindividualargs} class.
#'
#' Constructor function for the \code{insilicoindividualargs} class, with default values for the parameters if not provided by the user.
#'
#' @param ploidy Integer. Number of alleles for each gene (ploidy of the individuals). Default value is 2.
#' @param ngenevariants Integer. Number of alleles existing for each gene and segregating in the in silico population. Default value is 5.
#' @param qtleffect_samplingfct Function from which is sampled the value of a QTL effect coefficient (input x is the required sample size). Default value is a truncated normal distribution with mean 1 and sd 0.1 (only gives positive values).
#' @param initvar_samplingfct Function from which is sampled the variation of the initial abundance of a species (input x is the required sample size). Default value is a truncated normal distribution with mean 1 and sd 0.1 (only gives positive values).
#' @return An object of the class \code{insilicoindividualargs}, that is a named list of the different parameters.
#' @examples
#' indargs = insilicoindividualargs(ploidy = 4)
#' @export
insilicoindividualargs <- function(
  ploidy = 2,
  ngenevariants = 5,
  qtleffect_samplingfct = function(x){truncnorm::rtruncnorm(x, a = 0, b = Inf, mean = 1, sd = 0.1)},
  initvar_samplingfct = function(x){truncnorm::rtruncnorm(x, a = 0, b = Inf, mean = 1, sd = 0.1)}
){
  gcnList = sapply(1:ploidy, function(x){paste0("GCN", x)})
  ## qtlnames: names of the qtl effect coefficients
  ## The first 5 are the qtl affecting all genes, the last 5 only affect protein coding genes
  qtlnames = c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDregrate", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregrate", "qtlPTMregrate")

  value = list("ploidy" = ploidy,
               "gcnList" = gcnList,
               "ngenevariants" = ngenevariants,
               "qtleffect_samplingfct" = qtleffect_samplingfct,
               "initvar_samplingfct" = initvar_samplingfct,
               "qtlnames" = qtlnames
               )


  attr(value, "class") = "insilicoindividualargs"

  return(value)
}

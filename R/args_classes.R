## insilicosystemsargs class ----

#' Constructor function for the \code{insilicosystemsargs} class.
#'
#' Constructor function for the \code{insilicosystemsargs} class, with default values if not provided by the user.
#'
#' @param G Integer. Number of genes in the system. Default value is 10.
#' @param PC.p Numeric. Probability of each gene to be a protein-coding gene. Default value is 0.7.
#' @param PC.TC.p Numeric. Ratio of regulators of transcription among the protein-coding genes. Default value is 0.4.
#' @param PC.TL.p Numeric. Ratio of regulators of translation among the protein-coding genes. Default value is 0.3.
#' @param PC.RD.p Numeric. Ratio of regulators of RNA decay among the protein-coding genes. Default value is 0.1.
#' @param PC.PD.p Numeric. Ratio of regulators of protein decay among the protein-coding genes. Default value is 0.1.
#' @param PC.PTM.p Numeric. Ratio of regulators of protein post-translational modification among the protein-coding genes. Default value is 0.05.
#' @param PC.MR.p Numeric. Ratio of metabolic enzymes among the protein-coding genes. Default value is 0.05.
#' @param NC.TC.p Numeric. Ratio of regulators of transcription among the noncoding genes. Default value is 0.3.
#' @param NC.TL.p Numeric. Ratio of regulators of translation among the noncoding genes. Default value is 0.3.
#' @param NC.RD.p Numeric. Ratio of regulators of RNA decay among the noncoding genes. Default value is 0.3.
#' @param NC.PD.p Numeric. Ratio of regulators of protein decay among the noncoding genes. Default value is 0.05.
#' @param NC.PTM.p Numeric. Ratio of regulators of protein post-translational modification among the noncoding genes. Default value is 0.05.
#' @param TC.pos.p Numeric. Probability that the transcription is positively regulated by regulators. Default value is 0.5.
#' @param TL.pos.p Numeric. Probability that the translation is positively regulated by regulators. Default value is 0.5.
#' @param PTM.pos.p Numeric. Probability that the regulators transform the original protein into its modified form
#' (as opposed to transforming the modified protein back into its original form). Default value is 0.5.
# #' @param TC.PC.pos.p Numeric. Probability that the transcription is positively regulated by protein regulators. Default value is 0.5.
# #' @param TC.NC.pos.p Numeric. Probability that the transcription is positively regulated by noncoding regulators. Default value is 0.5.
# #' @param TL.PC.pos.p Numeric. Probability that the translation is positively regulated by protein regulators. Default value is 0.5.
# #' @param TL.NC.pos.p Numeric. Probability that the translation is positively regulated by noncoding regulators. Default value is 0.5.
# #' @param PTM.PC.pos.p Numeric. Probability that the protein regulators transform the original protein into its modified form
# #' (as opposed to transforming the modified protein back into its original form). Default value is 0.5.
# #' @param PTM.NC.pos.p Numeric. Probability that the noncoding regulators transform the original protein into its modified form
# #' (as opposed to transforming the modified protein back into its original form). Default value is 0.5.
#' @param basal_transcription_rate_samplingfct Function from which the transcription rates of genes are sampled (input x is the required sample size). Default value is
#' Values from Schwanhausser et al., 2013: transcription rate distribution log-normal, from 0.1 to 100 mRNA/hour -> we want the transcription rate in seconds. Default value is
#' @param basal_translation_rate_samplingfct Function from which the translation rates of genes are sampled (input x is the required sample size). Default value is
#' Values from Schwanhausser et al., 2013: translation rate distribution log-normal, from 0.1 to 10^5 protein/mRNA/hour -> we want the transcription rate in seconds. Default value is
#' @param basal_RNAlifetime_samplingfct Function from which the transcript lifetimes are sampled (input x is the required sample size). Default value is
#' Values from Schwanhausser et al., 2013: RNAs half-life distribution log-normal (chose boundary values to be from 1 min to 100 hours) -> we want the half-life in seconds. Default value is
#' @param basal_protlifetime_samplingfct Function from which the protein lifetime are sampled (input x is the required sample size). Default value is
#' Values from Schwanhausser et al., 2013: proteins half-life distribution log-normal (chose boundary values to be from 1 hour to 1000 hours) -> we want the half-life in seconds. Default value is
#' @param TC.PC.outdeg.distr Form of the distribution of the number of targets of transcription factors; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TC.NC.outdeg.distr Form of the distribution of the number of targets of noncoding RNAs regulating transcription; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TC.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the transcription regulation graph. Default value is 3.
#' @param TC.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the transcription regulation graph. Default value is 1.
#' @param TC.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the transcription regulation graph; can be either "powerlaw" or "exponential"). Default value is "powerlaw".
#' @param TC.NC.indeg.distr Type of preferential attachment for the targets of ncRNAs in the transcription regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TC.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation in the transcription regulation graph. Default value is 0.2.
#' @param TC.NC.autoregproba Numeric. Probability of ncRNAs to perform autoregulation in the transcription regulation graph. Default value is 0.
#' @param TC.PC.twonodesloop Logical. Are 2-nodes loops authorised in the transcription regulation graph with protein regulators? Default value is FALSE.
#' @param TC.NC.twonodesloop Logical. Are 2-nodes loops authorised in the transcription regulation graph with noncoding regulators? Default value is FALSE.
#' @param TCbindingrate_samplingfct Function from which the binding rate of transcription regulators on target are sampled (input x is the required sample size). Default value is
#' @param TCunbindingrate_samplingfct Function from which the unbinding rate of transcription regulators from target are sampled (input x is the required sample size). Default value is
#' @param TCfoldchange_samplingfct Function from which the transcription fold change induced by a bound regulator are sampled (input x is the required sample size). Default value is
#' @param TL.PC.outdeg.distr Form of the distribution of the number of targets of translation factors; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TL.NC.outdeg.distr Form of the the distribution of the number of targets of noncoding RNAs regulating translation; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TL.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the translation regulation graph. Default value is 3.
#' @param TL.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the translation regulation graph. Default value is 1.
#' @param TL.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the translation regulation graph can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TL.NC.indeg.distr Type of preferential attachment for the targets of ncRNAs in the translation regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param TL.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation in the translation regulation graph. Default value is 0.2.
#' @param TL.NC.autoregproba Numeric. Probability of ncRNAs to perform autoregulation in the translation regulation graph. Default value is 0.
#' @param TL.PC.twonodesloop Logical. Are 2-nodes loops authorised in the translation regulation graph with protein regulators? Default value is FALSE.
#' @param TL.NC.twonodesloop Logical. Are 2-nodes loops authorised in the translation regulation graph with noncoding regulators? Default value is FALSE.
#' @param TLbindingrate_samplingfct Function from which the binding rate of translation regulators on target are sampled (input x is the required sample size). Default value is
#' @param TLunbindingrate_samplingfct Function from which the unbinding rate of translation regulators from target are sampled (input x is the required sample size). Default value is
#' @param TLfoldchange_samplingfct Function from which the translation fold change induced by a bound regulator are sampled (input x is the required sample size). Default value is
#' @param RD.PC.outdeg.distr Form of the distribution of the number of targets of RNA decay factors; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param RD.NC.outdeg.distr Form of the the distribution of the number of targets of noncoding RNAs regulating RNA decay; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param RD.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the RNA decay regulation graph. Default value is 3.
#' @param RD.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the RNA decay regulation graph. Default value is 1.
#' @param RD.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the RNA decay graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param RD.NC.indeg.distr Type of preferential attachment for the targets of ncRNAs in the RNA decay graph;can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param RD.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation in the RNA decay regulation graph. Default value is 0.2.
#' @param RD.NC.autoregproba Numeric. Probability of ncRNAs to perform autoregulation in the RNA decay regulation graph. Default value is 0.
#' @param RD.PC.twonodesloop Logical. Are 2-nodes loops authorised in the RNA decay regulation graph with protein regulators? Default value is FALSE.
#' @param RD.NC.twonodesloop Logical. Are 2-nodes loops authorised in the RNA decay regulation graph with noncoding regulators? Default value is FALSE.
#' @param RDregrate_samplingfct Function from which the RNA decay rates of targets of RNA decay regulators are sampled (input x is the required sample size). Default value is
#' @param PD.PC.outdeg.distr Form of the distribution of the number of targets of protein decay factors; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PD.NC.outdeg.distr Form of the the distribution of the number of targets of noncoding RNAs regulating protein decay; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PD.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the protein decay regulation graph. Default value is 3.
#' @param PD.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the protein decay regulation graph. Default value is 1.
#' @param PD.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the protein decay regulation graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PD.NC.indeg.distr Type of preferential attachment for the targets of ncRNAs in the protein decay graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PD.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation in the protein decay regulation graph. Default value is 0.2.
#' @param PD.NC.autoregproba Numeric. Probability of ncRNAs to perform autoregulation in the protein decay regulation graph. Default value is 0.
#' @param PD.PC.twonodesloop Logical. Are 2-nodes loops authorised in the protein decay graph with protein regulators in the protein decay regulation graph? Default value is FALSE.
#' @param PD.NC.twonodesloop Logical. Are 2-nodes loops authorised in the protein decay graph with noncoding regulators in the protein decay regulation graph? Default value is FALSE.
#' @param PDregrate_samplingfct Function from which the protein decay rates of targets of protein decay regulators are sampled (input x is the required sample size). Default value is
#' @param PTM.PC.outdeg.distr Form of the distribution of the number of targets of protein post-translational modification regulators; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PTM.NC.outdeg.distr Form of the the distribution of the number of targets of noncoding RNAs regulating protein post-translational modification; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PTM.PC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the protein regulators in the protein post-translational modification graph. Default value is 3.
#' @param PTM.NC.outdeg.exp Numeric. Exponent of the distribution for the out-degree of the noncoding regulators in the protein post-translational modification graph. Default value is 1.
#' @param PTM.PC.indeg.distr Type of preferential attachment for the targets of protein regulators in the protein post-translational modification graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PTM.NC.indeg.distr Type of preferential attachment for the targets of ncRNAs in the protein post-translational modification graph; can be either "powerlaw" or "exponential". Default value is "powerlaw".
#' @param PTM.PC.autoregproba Numeric. Probability of protein regulators to perform autoregulation. Default value is 0.2.
#' @param PTM.NC.autoregproba Numeric. Probability of ncRNAs to perform autoregulation. Default value is 0.
#' @param PTM.PC.twonodesloop Logical. Are 2-nodes loops authorised in the protein post-translational modification graph with protein regulators? Default value is FALSE.
#' @param PTM.NC.twonodesloop Logical. Are 2-nodes loops authorised in the protein post-translational modification graph with noncoding regulators? Default value is FALSE.
#' @param PTMregrate_samplingfct Function from which the protein transformation rates of targets of post-translational modification regulators are sampled (input x is the required sample size). Default value is
#' @param regcomplexes Can the regulators controlling a common target form regulatory complexes in the protein post-translational modification graph? Can be 'none', 'prot' (only protein can form regulatory complexes) or 'both'
#' (both RNAs and proteins can form regulatory complexes). Default value is "prot".
#' @param regcomplexes.p Numeric. Probability that regulators controlling a common target form regulatory complexes; ignore if regcomplexes = 'none'. Default value is 0.3.
#' @param regcomplexes.size Integer. Number of components of a regulatory complex; ignore if regcomplexes = 'none'. Default value is 2.
#' @param complexesformationrate_samplingfct Function from which the formation rate of regulatory complexes are sampled (input x is the required sample size). Default value is
#' @param complexesdissociationrate_samplingfct Function from which the dissociation rate of regulatory complexes are sampled (input x is the required sample size). Default value is
#' @param mycolsCS Named vector of colour names or code. Colours used in the plots to represent portein-coding and noncoding genes. Default value is c("PC" = "#e03616",  "NC" = "#58355e", "Tot" = "#31161F").
#' @param mycolsGF Colours used in the plots to represent the different gene expression steps (transcription, translation, etc). Default value is c("TC" = "#FF7F11", "TL" = "#FF963C", "RD" = "#5AB7A4", "PD" = "#78C4B4", "PTM" = "#FF1B1C", "MR" = "#FF6D6E").
#' @param mycolsPosNeg Colours used in plots to represent positive/negative regulatory interactions. Default value is c("1" = "#D63230", "-1" = "#69BAF4").
#' @return An object of the class \code{insilicosystemargs}, that is a list of the different parameters.
#' @export
insilicosystemargs <- function(
  G = 10,
  PC.p = 0.7,
  PC.TC.p = 0.4,
  PC.TL.p = 0.3,
  PC.RD.p = 0.1,
  PC.PD.p = 0.1,
  PC.PTM.p = 0.05,
  PC.MR.p = 0.05,
  NC.TC.p = 0.3,
  NC.TL.p = 0.3,
  NC.RD.p = 0.3,
  NC.PD.p = 0.05,
  NC.PTM.p = 0.05,
  TC.pos.p = 0.5,
  TL.pos.p = 0.5,
  PTM.pos.p = 0.5,
#  TC.PC.pos.p = 0.5,
#  TC.NC.pos.p = 0.5,
#  TL.PC.pos.p = 0.5,
#  TL.NC.pos.p = 0.5,
#  PTM.PC.pos.p = 0.5,
#  PTM.NC.pos.p = 0.5,
  basal_transcription_rate_samplingfct = function(x){ logval = rnorm(x, mean = 0.5, sd = 0.5); val = 10^logval; return(val/3600) },
  basal_translation_rate_samplingfct = function(x){ logval = rnorm(x, mean = 2.5, sd = 0.8); val = 10^logval; return(val/3600) },
  basal_RNAlifetime_samplingfct = function(x){ logval = rnorm(x, mean = 0, sd = 0.5); val = 10^logval; return(val*3600) },
  basal_protlifetime_samplingfct = function(x){ logval = rnorm(x, mean = 1.75, sd = 0.5); val = 10^logval; return(val*3600) },
  TC.PC.outdeg.distr = "powerlaw",
  TC.NC.outdeg.distr = "powerlaw",
  TC.PC.outdeg.exp = 3,
  TC.NC.outdeg.exp = 1,
  TC.PC.indeg.distr = "powerlaw",
  TC.NC.indeg.distr = "powerlaw",
  TC.PC.autoregproba = 0.2,
  TC.NC.autoregproba = 0,
  TC.PC.twonodesloop = FALSE,
  TC.NC.twonodesloop = FALSE,
  TCbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  TCunbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  TCfoldchange_samplingfct = function(x){ truncnorm::rtruncnorm(x, a = 1, mean = 3, sd = 10) },
  TL.PC.outdeg.distr = "powerlaw",
  TL.NC.outdeg.distr = "powerlaw",
  TL.PC.outdeg.exp = 3,
  TL.NC.outdeg.exp = 1,
  TL.PC.indeg.distr = "powerlaw",
  TL.NC.indeg.distr = "powerlaw",
  TL.PC.autoregproba = 0.2,
  TL.NC.autoregproba = 0,
  TL.PC.twonodesloop = FALSE,
  TL.NC.twonodesloop = FALSE,
  TLbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  TLunbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  TLfoldchange_samplingfct = function(x){ sample(2:30, x, replace = T)  },
  RD.PC.outdeg.distr = "powerlaw",
  RD.NC.outdeg.distr = "powerlaw",
  RD.PC.outdeg.exp = 3,
  RD.NC.outdeg.exp = 1,
  RD.PC.indeg.distr = "powerlaw",
  RD.NC.indeg.distr = "powerlaw",
  RD.PC.autoregproba = 0.2,
  RD.NC.autoregproba = 0,
  RD.PC.twonodesloop = FALSE,
  RD.NC.twonodesloop = FALSE,
  RDregrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  PD.PC.outdeg.distr = "powerlaw",
  PD.NC.outdeg.distr = "powerlaw",
  PD.PC.outdeg.exp = 3,
  PD.NC.outdeg.exp = 1,
  PD.PC.indeg.distr = "powerlaw",
  PD.NC.indeg.distr = "powerlaw",
  PD.PC.autoregproba = 0.2,
  PD.NC.autoregproba = 0,
  PD.PC.twonodesloop = FALSE,
  PD.NC.twonodesloop = FALSE,
  PDregrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  PTM.PC.outdeg.distr = "powerlaw",
  PTM.NC.outdeg.distr = "powerlaw",
  PTM.PC.outdeg.exp = 3,
  PTM.NC.outdeg.exp = 1,
  PTM.PC.indeg.distr = "powerlaw",
  PTM.NC.indeg.distr = "powerlaw",
  PTM.PC.autoregproba = 0.2,
  PTM.NC.autoregproba = 0,
  PTM.PC.twonodesloop = FALSE,
  PTM.NC.twonodesloop = FALSE,
  PTMregrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  regcomplexes = "prot",
  regcomplexes.p = 0.3,
  regcomplexes.size = 2,
  complexesformationrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  complexesdissociationrate_samplingfct = function(x){ runif(x, 0.001, 0.01) },
  mycolsCS = c("PC" = "#e03616",  "NC" = "#58355e", "Tot" = "#31161F"),
  mycolsGF = c("TC" = "#FF7F11", "TL" = "#FF963C", "RD" = "#5AB7A4", "PD" = "#78C4B4", "PTM" = "#FF1B1C", "MR" = "#FF6D6E"),
  mycolsPosNeg = c("1" = "#D63230", "-1" = "#69BAF4")
){
  NC.p = 1 - PC.p

  temp = sum(PC.TC.p + PC.TL.p + PC.RD.p + PC.PTM.p + PC.MR.p)
  PC.TC.p = PC.TC.p/temp
  PC.TL.p = PC.TL.p/temp
  PC.RD.p = PC.RD.p/temp
  PC.PD.p = PC.PD.p/temp
  PC.PTM.p = PC.PTM.p/temp
  PC.MR.p = PC.MR.p/temp

  temp = sum(NC.TC.p + NC.TL.p + NC.RD.p + NC.PTM.p)
  NC.TC.p = NC.TC.p/temp
  NC.TL.p = NC.TL.p/temp
  NC.RD.p = NC.RD.p/temp
  NC.PD.p = NC.PD.p/temp
  NC.PTM.p = NC.PTM.p/temp

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
                 "complexesdissociationrate_samplingfct" = complexesdissociationrate_samplingfct,
                 "mycolsCS" = mycolsCS,
                 "mycolsGF" = mycolsGF,
                 "mycolsPosNeg" = mycolsPosNeg)

  attr(value, "class") = "insilicosystemargs"
  return(value)
}


## insilicoindividualargs class ----

#' Constructor function for the \code{insilicoindividualargs} class.
#'
#' Constructor function for the \code{insilicoindividualargs} class, with default values if not provided by the user.
#'
#' @param ploidy Integer. Number of alleles for each gene. Default value is 2.
#' @param ngenevariants Integer. Number of alleles existing for each gene. Default value is 5.
#' @param qtleffect_samplingfct Function from which is sampled the effect of a QTL (input x is the required sample size). Default value is a truncated normal distribution with mean 1 and sd 0.1 (only gives positive values).
#' @param initvar_samplingfct Function from which is sampled the variation in initial abundance of a species (input x is the required sample size). Default value is a truncated normal distribution with mean 1 and sd 0.1 (only gives positive values).
#' @return An object of the class \code{insilicoindividualargs}, that is a list of the different parameters.
#' @export
insilicoindividualargs <- function(
  ploidy = 2,
  ngenevariants = 5,
  qtleffect_samplingfct = function(x){truncnorm::rtruncnorm(x, a = 0, b = Inf, mean = 1, sd = 0.1)},
  initvar_samplingfct = function(x){truncnorm::rtruncnorm(x, a = 0, b = Inf, mean = 1, sd = 0.1)}
){
  gcnList = sapply(1:ploidy, function(x){paste0("GCN", x)})
  ## qtlnames: names of the qtl effect coefficients
  ## The first 5 are the qtl affecting all genes, the last 4 only affect protein coding genes
  qtlnames = c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind")

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

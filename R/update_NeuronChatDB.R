#' Update NeuronChatDB using user defined new entry information
#'
#' Update NeuronChatDB using user defined new entry information
#'
#' @param DB current interation database (stored in slot 'DB' of NeuronChat object, i.e., object@DB)
#' @param interaction_name a character to denote the name of the interaction to be added
#' @param lig_contributor a vector of gene symbols, contributing to ligand abundance
#' @param receptor_subunit a vector of gene symbols, contributing to receptor abundance
#' @param interaction_type a character to denote the interaction type
#' @param ligand_type a character to denote the ligand_type
#' @param lig_contributor_group a vector of integers to indicate the biological function group of genes in `lig_contributor`;
#' the length of `lig_contributor_group` is the same as `lig_contributor`
#' @param lig_contributor_coeff a numerical vector to indicate the stoichiometric coefficients of ligand-contribuing groups, to calculate ligand abundance;
#' the length is the same as the unique of `lig_contributor_group`
#' @param receptor_subunit_group a vector of integers to indicate the biological function group of genes in `receptor_subunit`;
#' the length of `receptor_subunit_group` is the same as `receptor_subunit`
#' @param receptor_subunit_coeff a numerical vector to indicate the stoichiometric coefficients of receptor subunit groups, to calculate receptor abundance;
#' the length is the same as the unique of `receptor_subunit_group`
#' @return a updated interaction database (a list)
#' @export
update_interactionDB <- function(DB,interaction_name,lig_contributor,receptor_subunit,interaction_type='user_defined',ligand_type='user_defined',
                                 lig_contributor_group=NULL,lig_contributor_coeff=NULL,receptor_subunit_group=NULL,receptor_subunit_coeff=NULL){
DB_new <-DB
ll <- length(DB_new)
if(interaction_name %in% names(DB)){
  stop("The input `interaction_name` already exists in the current DB")
} else {
  interaction_name <- interaction_name
  lig_contributor <- lig_contributor
  receptor_subunit <- receptor_subunit
  if(is.null(lig_contributor_group)){
    lig_contributor_group <- rep(1,length(lig_contributor))
    lig_contributor_coeff <- 1
  }
  if(is.null(receptor_subunit_group)){
    receptor_subunit_group <- rep(1,length(receptor_subunit))
    receptor_subunit_coeff <- 1
  }
  if(length(lig_contributor_group)!=length(lig_contributor) | length(receptor_subunit_group)!=length(receptor_subunit)){
    stop("The input `lig_contributor_group` (or `receptor_subunit_group`) doesn't match the length of `lig_contributor` (or `receptor_subunit`)")
  }
  if(length(lig_contributor_coeff)!=length(unique(lig_contributor_group)) | length(receptor_subunit_coeff)!=length(unique(receptor_subunit_group))){
    stop("The input `lig_contributor_coeff` (or `receptor_subunit_coeff`) doesn't match the length of `unique(lig_contributor_group)` (or `unique(receptor_subunit_group)`")
  }
  interaction_tmp <- list(interaction_name=interaction_name,lig_contributor=lig_contributor,receptor_subunit=receptor_subunit,
                          lig_contributor_group=lig_contributor_group,lig_contributor_coeff=lig_contributor_coeff,
                          receptor_subunit_group=receptor_subunit_group,receptor_subunit_coeff=receptor_subunit_coeff,
                          targets_up=character(0), targets_down=character(0), activator=character(0),inhibitor=character(0), interactors=character(0),
                          interaction_type=interaction_type,ligand_type=ligand_type)
}
DB_new[[ll+1]] <- interaction_tmp
names(DB_new)[ll+1] <- interaction_name
return(DB_new)
}

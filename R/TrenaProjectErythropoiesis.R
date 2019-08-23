#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title TrenaProjectErythropoiesis-class
#'
#' @name TrenaProjectErythropoiesis-class
#' @rdname TrenaProjectErythropoiesis-class
#' @aliases TrenaProjectErythropoiesis
#' @exportClass TrenaProjectErythropoiesis
#'

.TrenaProjectErythropoiesis <- setClass("TrenaProjectErythropoiesis",
                                  contains="TrenaProjectHG38")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectErythropoiesis
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectErythropoiesis-class
#'
#' @export
#'
#' @return An object of the TrenaProjectErythropoiesis class
#'

TrenaProjectErythropoiesis <- function(quiet=TRUE)

{
   directory <- system.file(package="TrenaProjectErythropoiesis", "extdata", "geneSets")
   geneSet.files <- list.files(directory)
   geneSets <- list()

   for(file in geneSet.files){
      full.path <- file.path(directory, file)
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE)
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   footprintDatabaseNames <- NA_character_;
   footprintDatabaseHost <- NA_character_;
   footprintDatabasePort <- NA_integer_;
   dataDirectory <- system.file(package="TrenaProjectErythropoiesis", "extdata")

   covariatesFile <- NA_character_;

   stopifnot(file.exists(dataDirectory))

   .TrenaProjectErythropoiesis(TrenaProjectHG38(projectName="TrenaProjectErythropoiesis",
                                                supportedGenes=geneSets[[1]],
                                                footprintDatabaseHost=footprintDatabaseHost,
                                                footprintDatabasePort=footprintDatabasePort,
                                                footprintDatabaseNames=footprintDatabaseNames,
                                                packageDataDirectory=dataDirectory,
                                                quiet=quiet
                                                ))

} # TrenaProjectErythropoiesis, the constructor
#----------------------------------------------------------------------------------------------------

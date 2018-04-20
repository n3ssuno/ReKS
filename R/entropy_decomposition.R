#' @name
#' entropy_decomposition
#'
#' @title
#' Entropy Decoposition Theorem
#'
#' @description
#' The function will return the information entropy of a given vector of
#' frequences decomposed in two parts: a between-groups and a within-groups one.
#' Moreover, it provides you also the probability of each group and the
#' entropy of each of the groups.
#'
#' @details The total (undecomposed) entropy is equal to the sum of the two
#' decomposed components. The within-groups component is equal to the weighted
#' average of the entropy of each of the groups. So the between-groups component
#' is equal to the residual entropy.
#' See:
#' \itemize{
#' \item{Theil (1972) \emph{Statistical Decomposition Analysis}, North-Holland;}
#' \item{Zadjenweber (1972) "Une Application de la Th{\'e}orie de l’Information
#' {\`a} l'{\'E}conomie: La Mesure de la Concentration, \emph{Revue d'Economie
#' Politique}, 82, 486–-510;}
#' \item{Attaran (1985) "Industrial Diversity and Economic Performance in U.S.
#' Areas, \emph{The Annals of Regional Science}, 20, 44--54;}
#' \item{Frenken (2007) "Entropy statistics and information theory", in
#' Hanusch and Pyka (Eds.) \emph{Elgar Companion to Neo-Schumpeterian
#' Economics}, Edward Elgar.}
#' \item{Frenken, van Oort and Verburg (2007) "Related Variety, Unrelated
#' Variety and Regional Economic Growth", \emph{Regional Studies}, 41,
#' 685--697;}
#' \item{Quatraro (2010) 'Knowledge Coherence, Variety and Economic Growth:
#' Manufacturing Evidence from Italian Regions', \emph{Research Policy}, 39,
#' 1289--1302;}
#' \item{Rocchetta and Mina (2017), "Technological Coherence and the Adaptive
#' Resilience of Regional Economies", \emph{Papers in Evolutionary Economic
#' Geography}, Utrecht University.}
#' }
#'
#' @param data Vector of elements.
#' @param groups Vector that describes how the elements in the data can be
#' aggregated in groups.
#' @return A list with four elements. The information entropy level of the vector.
#' @examples
#' etpy_decomp <- entropy_decomposition(data, grps)
#' etpy_decomp <- entropy_decomposition(table(data), grps)

entropy_decomposition <- function(data, groups) {
    Pg <- by(.get_freqs(data), groups, sum)
    BG <- RKS::entropy(Pg)
    WG <- RKS::entropy(data) - BG
    by_group <- log2(Pg) + 1/Pg * by(data, groups, RKS::entropy)
    etp_dcp <- list(BG = BG,
                    WG = WG,
                    by_group = by_group,
                    Pg = Pg)
    return(etp_dcp)
}

# TODO
# 1. It takes inputs different than entropy_decomposition_panel
#     because in this case it wants two arrays, while in the other case
#     it wants a database and a number of column names
# 2. Also the output of the two functions is different. Here it returns
#     a list of objects, while in the other case a dataframe

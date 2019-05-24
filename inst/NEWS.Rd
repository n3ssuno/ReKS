\name{NEWS}
\title{ReKS News}
\section{Changes in ReKS version 0.5.0 (2019-05-24)}{
    %\subsection{BUG FIXES}{
    %    \itemize{
    %        \item
    %    }
    %}
    \subsection{NEW FEATURES}{
        \itemize{
            \item Added two completely new functions:
                \itemize{
                    \item \code{coi()} that computes the Complexity Outlook
                    Index for each region considered.
                    \item \code{cog()} that computes the Complexity Outlook Gain
                    for each technological domain in each region considered.
                }
                The documentation for these two functions will be added soon.
        }
    }
    \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
        \itemize{
            \item All the functions has been re-written. Follow the example in
            the helps to understand the basic logic of the package.
        }
    }
    \subsection{DEPRECATED AND DEFUNCT}{
        \itemize{
            \item \code{choropleth_map()}
            \item \code{complexity_hh()}
            \item \code{entropy()}
            \item \code{entropy_decomposition()}
            \item \code{entropy_decomposition_panel()}
            \item \code{fitness_tccgp()}
            \item \code{image_biadjacency_matrix()}
            \item \code{jaccard_index()}
            \item \code{occurrence_matrix()}

        }
    }
}
\section{Changes in ReKS version 0.3.0 (2018-07-12)}{
    \subsection{BUG FIXES}{
        \itemize{
            \item \code{complexity_hh()} and \code{complexity_hh_panel()} now
            add a NA when the eigenvector of the second largest eigenvalue of
            the adjacency matrix (thanks to M. Tubiana for letting me know).
        }
    }
    \subsection{NEW FEATURES}{
        \itemize{
            \item Added two new couples of functions:
                \itemize{
                    \item \code{proximity_hkbh()} and
                    \code{proximity_hkbh_panel()} compute the proximity between
                    the different pieces of regional knowledge following the
                    methodology proposed by Hidalgo \emph{et al.} (2007);
                    \item \code{proximity_zctp()} and
                    \code{proximity_zctp_panel()} do the same, but following the
                    methodology developed by Zaccaria \emph{et al.} (2014).
                }
            \item Added a function to plot the
            \bold{maximum_spanning_tree}-class objects obtained from the
            two new proximity functions above said.
        }
    }
    \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
        \itemize{
            \item The functions \code{complexity()} and \code{fitness()} (now
            deprecated, see below) have been splitted in two to fit the
            standards of the package, and the names have been modified to
            describe more properly their content. So now you have
            \code{complexity_hh()} and \code{complexity_hh_panel()}; and
            \code{fitness_tccgp()} and \code{fitness_tccgp_panel()}. They do,
            respectively, the same things of the deprecated ones.
            \item Besides other information, now you can find the
            ``diversification'' and the ``ubiquity'' of each element (region or
            patent class) in each year in the attributes of the objects of class
            \bold{reks_hh_complexity}, \bold{reks_fitness_tccgp},
            \bold{reks_proximity_hkbh}, and \bold{reks_proximity_zctp}. See
            \code{attr(<object>, 'diversification')} and
            \code{attr(<object>, 'ubiquity')}. In the attributes of
            the he objects of class \bold{reks_fitness_tccgp} you can also find
            the number of iterations that the iterative process took to
            converge, and if it converged.
        }
    }
    \subsection{DEPRECATED AND DEFUNCT}{
        \itemize{
            \item \code{complexity()} and \code{fitness()} are deprecated (with
            a warning).
        }
    }
}

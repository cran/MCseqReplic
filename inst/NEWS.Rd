\name{NEWS}
\title{NEWS file of the MCseqReplic package}

\section{Changes in version 1.0.0 [2026-04-29] }{

  \subsection{Changes in existing functions}{
    \itemize{
      \item
      \code{seqdistMCSE} argument \code{silent} replaced by an argument \code{verbose}
    }
  }
  \subsection{Misc}{
    \itemize{
      \item
      Package released on CRAN.
      \item
      CITATION file.
    }
  }

}

\section{Changes in version 0.3.4 [2026-02-21] }{
  \subsection{New function}{
    \itemize{
      \item
      \code{MCpj()} Generates a distribution Pj of the timing error probability that complies with expected size of timing errors and the probability of no error.
    }
  }
}


\section{Changes in version 0.3.3 [2025-10-10] }{
  \subsection{Changes in existing functions}{
    \itemize{
      \item
      \code{MCseqdistSE()} The standards error of each distance (\code{MC.se}) is now the square root of the sum of the variance of the simulated distances and the squared bias. The standard deviation of the distances is returned as \code{MC.sd}. Plot and summary methods have been adjusted conformably.
      \item
      \code{MCratios()} z-ratios now take account of the new definition of \code{MC.se}
    }
  }
}

\section{Changes in version 0.3.2 [2025-03-18] }{
    \subsection{Changes in existing functions}{
      \itemize{
        \item
        \code{MCdisslist()} new logical argument \code{use.udiss}, which when set as \code{TRUE}, tells the function to derive MC-distances from the distances between unique merged replicated sequences.
        \item
        \code{MCudist()} adapted for distances to reference sequences, which can be passed with new argument \code{seqref}.
      }

    }

    \subsection{Bug fixes}{
      \itemize{
        \item
        \code{MCclustqual} option clustmeth="ward.D" (or any other \code{hclust} method) did not work.
        \item
        \code{MCnunique} wrong result of \code{check = TRUE}.
      }
    }
}

\section{Changes in version 0.3.1 [2025-03-06] }{
    \subsection{New functions}{
      \itemize{
        \item
        \code{MCdisscorr()} Correlations between observed and MC-simulated distances.
        \item
        \code{MCclustqual()} Cluster quality for a range of partition sizes per MC-set.
      }

    }
    \subsection{Changes in existing functions}{
      \itemize{
        \item
        \code{MCdisslist()} gains argument \code{seqref} to specify reference sequences to analyze distances to these references.
        \item
        \code{MCseqdistSE()} and associated \code{print} and \code{summary} methods now supports list of matrices of distances to reference sequences.
        \item
        \code{MCratios()} and associated \code{print} and \code{summary} method now also apply to outcome of \code{MCseqdistSE()} (class \code{distMC)} for distances to reference sequences.
        \item
        \code{summary.MCratios()} gains a threshold argument \code{thresh} and returns the number and proportion of \code{diss.z} ratios lower than this threshold.
        \item
        \code{MCclustcomp} now returns all the indexes provided by package \code{aricode}  for each MC-set.
      }
    }
    \subsection{Bug fixes}{
        \itemize{
            \item
            \code{MCseqReplicate()} error when all dss have same length.
        }
    }
}

\section{Changes in version 0.3.0 [2025-01-10] }{
    \subsection{Misc.}{
        \itemize{
            \item
            Packaged created from seqdistMCSE 0.2.0
            \item
            Implements the generation of replicated sets of sequences with timing changes.
        }
    }
}


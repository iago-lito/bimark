# Here is R script to mock object organization of the bimark package

# convenience string for the scripts:
bmclass <- "BimarkModel"

#' Create the data object that'll endure all this package processing
#'
#' R's S4 object opportunities sound nice but.. terrible to use. Furthermore, we
#' do not need such a sophisticated object architecture. As an alternative, we
#' are falling back on this plain named-and-typed list which will gather
#' information as the analysis goes by. Items are simili-slots. Their values are
#' set to NULL if they are not defined.
#'
#' @slot L raw latent histories matrix if available (available if simulated)
#' @slot M raw observation matrix
#'
#' @seealso print.BimarkModel
#'
#' @examples
#' m <- createBimarkModel()
#'
#' @export

createBimarkModel <- function() {
  model <- list(
                L=NULL,
                M=NULL
                )
  class(model) <- bmclass
  return(model)
}

#' Print BimarkModel pseudo-object to the console
#'
#' @seealso createBimarkModel
#'
#' @examples
#' m <- createBimarkModel()
#' m
#' print(m)
#'
#' @export

print.BimarkModel <- function(m, ...) {
  # Consider it as empty if there are no raw histories matrices:
  if (is.null(m$L) && is.null(m$M))
    cat("empty", bmclass, '\n')
}


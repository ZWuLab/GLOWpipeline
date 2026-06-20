# ==============================================================================
# Internal Utilities
# ==============================================================================
#
# Package-internal helpers with no public surface.
#
# INTERNAL HELPER FUNCTIONS:
#   - `%||%`   NULL-coalescing operator (not exported)


#################### INTERNAL HELPER FUNCTIONS ####################

# Internal: NULL-coalescing operator. Not exported.
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

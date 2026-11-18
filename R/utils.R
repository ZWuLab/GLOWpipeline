# ==============================================================================
# Internal Utilities
# ==============================================================================
#
# Package-internal helpers with no public surface.
#
# INTERNAL HELPER FUNCTIONS:
#   - `%||%`          NULL-coalescing operator (not exported)
#   - .thin_cells()   Grid resolution of the display thinning, from the marker size
#   - .thin_points()  One drawn point per grid cell of the current plot, the highest kept


#################### INTERNAL HELPER FUNCTIONS ####################

# Internal: NULL-coalescing operator. Not exported.
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


# Internal: grid resolution of the display thinning. A copy of GLOWr's internal .thin_cells().
#
# A marker of pch 20 at size `cex` is drawn 0.035 inch wide at cex 0.5 and 0.060 inch at 1
# (GLOWr's measurement on a 1200 dpi raster), taken as 0.010 + 0.050 * cex inch. The cell is
# half that diameter, so two points in one cell overlap. The markers of the per-family QQ
# page (pch 16, 17, 15, 18, 8, 4) are at least as wide as pch 20, so the cell is fine enough
# for them.
#
# @param cells_per_inch Numeric or NULL. An explicit resolution, or NULL for the marker's.
# @param cex Numeric marker size.
# @return Numeric, cells per inch.
#' @noRd
.thin_cells <- function(cells_per_inch, cex) {
  if (!is.null(cells_per_inch)) {
    stopifnot(is.numeric(cells_per_inch), length(cells_per_inch) == 1L,
              is.finite(cells_per_inch), cells_per_inch > 0)
    return(cells_per_inch)
  }
  diam_in <- 0.010 + 0.050 * cex
  2 / diam_in
}


# Internal: one drawn point per grid cell, the highest kept. A copy of GLOWr's internal
# .thin_points().
#
# Maps the points to the cells of a grid over the plot region of the current frame
# (par("usr") and par("pin")) and keeps, per cell, the point with the largest y. A logarithmic
# axis is mapped as par("usr") is, and a point outside the plot region is kept without taking
# a cell, so a cropped axis hides nothing the full set would show.
#
# @param x,y Numeric user coordinates (untransformed; a log axis is handled here).
# @param cells_per_inch Numeric grid resolution.
# @return Logical vector, TRUE for the points to draw.
#' @noRd
.thin_points <- function(x, y, cells_per_inch) {
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  # par("usr") is on the log10 scale for a logarithmic axis; map the coordinates the same way
  if (isTRUE(graphics::par("xlog"))) x <- log10(x)
  if (isTRUE(graphics::par("ylog"))) y <- log10(y)
  nx <- max(1, ceiling(pin[1L] * cells_per_inch))
  ny <- max(1, ceiling(pin[2L] * cells_per_inch))
  # A point outside the plot region (a cropped axis) takes no cell and is kept: the device
  # clips it as it would without thinning, and it cannot hide a visible point in its cell.
  inside <- is.finite(x) & is.finite(y) &
            x >= usr[1L] & x <= usr[2L] & y >= usr[3L] & y <= usr[4L]
  keep <- !inside
  ix <- floor((x[inside] - usr[1L]) / (usr[2L] - usr[1L]) * nx)
  iy <- floor((y[inside] - usr[3L]) / (usr[4L] - usr[3L]) * ny)
  ix <- pmin(ix, nx - 1)                   # the upper edge belongs to the last cell
  iy <- pmin(iy, ny - 1)
  cell <- ix * ny + iy                     # one number per (ix, iy); exact in double
  yi <- y[inside]
  o <- order(yi, decreasing = TRUE)        # the highest point of a cell is the one kept
  k <- logical(length(yi))
  k[o[!duplicated(cell[o])]] <- TRUE
  keep[inside] <- k
  keep
}

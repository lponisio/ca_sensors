

## ------------------------------------------------------------
## Function: letter_to_number()
## ------------------------------------------------------------
##
## Purpose:
##   Converts row letters into numeric row positions.
##
##   For example:
##     A  -> 1
##     B  -> 2
##     R  -> 18
##     AA -> 27
##
## Inputs:
##   x: character vector of grid row labels.
##      Usually letters such as "A", "B", "C", ..., "R".
##
## Output:
##   A numeric vector giving the row number corresponding to each letter.
##
## Why this is needed:
##   GRTS needs spatial coordinates. Grid labels are categorical,
##   so we convert letters into row numbers and use those as y-coordinates.

letter_to_number <- function(x) {
  x <- toupper(as.character(x))

  vapply(strsplit(x, ""), function(chars) {
    vals <- match(chars, LETTERS)

    if (any(is.na(vals))) {
      stop("Grid row labels must contain letters only.")
    }

    sum(vals * 26^rev(seq_along(vals) - 1))
  }, numeric(1))
}

## ------------------------------------------------------------
## Function: parse_grid_label()
## ------------------------------------------------------------
##
## Purpose:
##   Splits grid cell labels into their numeric and letter parts.
##
##   Accepts either:
##     "1A", "2B", "9R"
##   or:
##     "A1", "B2", "R9"
##
## Inputs:
##   ids: character vector of grid cell labels.
##
## Output:
##   A data frame with:
##     grid_col_number: numeric grid column, e.g. 1, 2, 3
##     grid_row_letter: lettered grid row, e.g. A, B, C
##     grid_row_number: numeric row position, e.g. A = 1, B = 2
##
## Why this is needed:
##   We need to turn labels like "4C" into x/y coordinates.
##   Here, the number becomes the x-coordinate and the letter becomes
##   the y-coordinate.

parse_grid_label <- function(ids) {
  ids <- trimws(as.character(ids))

  # Check whether labels are in number-letter or letter-number format.
  number_letter <- grepl("^[0-9]+[A-Za-z]+$", ids)
  letter_number <- grepl("^[A-Za-z]+[0-9]+$", ids)

  if (any(!(number_letter | letter_number))) {
    stop("Expected grid labels like 1A, 2B, 9R, or A1, B2.")
  }

  grid_number <- rep(NA_integer_, length(ids))
  grid_letter <- rep(NA_character_, length(ids))

  # Parse labels such as 1A, 2B, 9R.
  grid_number[number_letter] <- as.integer(
    sub("^([0-9]+)([A-Za-z]+)$", "\\1", ids[number_letter])
  )

  grid_letter[number_letter] <- toupper(
    sub("^([0-9]+)([A-Za-z]+)$", "\\2", ids[number_letter])
  )

  # Parse labels such as A1, B2, R9.
  grid_letter[letter_number] <- toupper(
    sub("^([A-Za-z]+)([0-9]+)$", "\\1", ids[letter_number])
  )

  grid_number[letter_number] <- as.integer(
    sub("^([A-Za-z]+)([0-9]+)$", "\\2", ids[letter_number])
  )

  data.frame(
    grid_col_number = grid_number,
    grid_row_letter = grid_letter,
    grid_row_number = letter_to_number(grid_letter)
  )
}


## ------------------------------------------------------------
## Function: largest_remainder_capped()
## ------------------------------------------------------------
##
## Purpose:
##   Converts fractional sample-size targets into whole numbers,
##   while making sure the total sample size equals the requested n.
##
## Example:
##   Suppose we want 10 cells from a location and the vegetation
##   proportions imply:
##
##     forest = 4.4
##     shrub  = 3.3
##     meadow = 2.3
##
##   We cannot sample fractional cells, so this function rounds them
##   to whole numbers while preserving the total of 10.
##
## Inputs:
##   target: numeric vector of desired fractional sample sizes.
##
##   available: numeric vector of how many cells are actually available
##              in each group.
##
##   n: total number of cells requested.
##
## Output:
##   Integer vector of final sample sizes.
##
## Why this is needed:
##   If sampling proportionally across vegetation categories, the exact
##   proportional allocation often produces decimals.

largest_remainder_capped <- function(target, available, n) {
  # Start by taking the floor of each target.
  # Example: 4.4 becomes 4.
  counts <- pmin(floor(target), available)

  # Add one cell at a time until the total equals n.
  while (sum(counts) < n) {
    eligible <- which(counts < available)

    if (length(eligible) == 0) {
      stop("Requested more cells than are available.")
    }

    remainder <- target - floor(target)

    # Give the next cell to the group with the largest leftover fraction.
    # runif() breaks ties randomly.
    chosen <- eligible[order(-remainder[eligible], runif(length(eligible)))][1]

    counts[chosen] <- counts[chosen] + 1L
  }

  as.integer(counts)
}


## ------------------------------------------------------------
## Function: allocate_within_location()
## ------------------------------------------------------------
##
## Purpose:
##   Decides how many cells to sample from each vegetation category
##   within one location.
##
## Inputs:
##   strata: vector of vegetation categories for the cells in one location.
##
##   n: total number of cells to sample from that location.
##
##   allocation: how to allocate cells across vegetation categories.
##      Options:
##        "proportional" = sample vegetation types in proportion to availability
##        "equal"        = sample roughly the same number from each type
##        named vector   = manually specify sample size per vegetation type
##
##      Example named vector:
##        c("forest" = 4, "shrub" = 3, "meadow" = 3)
##
## Output:
##   Named integer vector giving the number of cells to sample from
##   each vegetation category.
##
## Why this is needed:
##   spsurvey::grts() can do stratified GRTS sampling, but it needs
##   the desired sample size for each stratum.

allocate_within_location <- function(strata, n, allocation = "proportional") {
  strata <- as.character(strata)

  available_tab <- table(strata)
  available <- as.integer(available_tab)
  names(available) <- names(available_tab)

  # If the user provides a named vector, use it directly.
  if (is.numeric(allocation)) {
    if (is.null(names(allocation))) {
      stop("Numeric allocation must be a named vector.")
    }

    out <- setNames(integer(length(available)), names(available))
    out[names(allocation)] <- as.integer(allocation)

    if (sum(out) != n) {
      stop("Named allocation must sum to the requested sample size for that location.")
    }

    if (any(out > available)) {
      stop("Requested more cells than available in at least one vegetation category.")
    }

    return(out)
  }

  allocation <- match.arg(allocation, c("proportional", "equal"))

  if (allocation == "proportional") {
    # Allocate cells according to vegetation availability.
    target <- n * available / sum(available)
  } else {
    # Allocate cells as evenly as possible across vegetation categories.
    target <- rep(n / length(available), length(available))
    names(target) <- names(available)
  }

  counts <- largest_remainder_capped(target, available, n)
  names(counts) <- names(available)

  counts
}


## ------------------------------------------------------------
## Function: choose_grts_cells_by_location() 
## ------------------------------------------------------------

## Purpose:
##   Performs separate GRTS draws for South Fork and Emerald Queen.
##
## Inputs:
##   cell_table:
##     Data frame with at least:
##       - grid cell label column, e.g. "cell"
##       - location column, e.g. "location"
##       - vegetation column, e.g. "vegetation"
##
##   n_by_location:
##     Named numeric vector giving the number of cells to sample per location.
##
##   cell_col:
##     Name of the column containing grid cell labels.
##     Default: "cell"
##
##   veg_col:
##     Name of the column containing vegetation categories.
##     Default: "vegetation"
##
##   location_col:
##     Name of the column containing location names.
##     Default: "location"
##
##   stratify_by_vegetation:
##     TRUE/FALSE.
##     If TRUE, the GRTS draw within each location is stratified by vegetation.
##     If FALSE, the function ignores vegetation and samples spatially balanced
##     cells within each location (likely we won't use this)
##
##   vegetation_allocation:
##     How to allocate samples across vegetation categories within each location.
##
##     Options:
##       "proportional" (to area/number of cells with this catagory)
##       "equal"
##       named vector of sample sizes
##
##   seed:
##     Optional random seed for reproducibility.
##
##   n_over:
##     Number of replacement/oversample cells to draw per location.
##     These can be used if selected cells are inaccessible temporarily.
##
## Output:
##   A list with four elements:
##
##   selected:
##     Data frame of primary selected cells to survey.
##
##   replacements:
##     Data frame of backup cells.
##     NULL if n_over = 0.
##
##   allocation:
##     For each location, shows how many cells were allocated to each
##     vegetation category, if stratifying.
##
##   design:
##     Raw GRTS design objects from spsurvey, one per location.
##
## Important assumptions:
##   1. Grid labels are formatted like 1A or A1.
##   2. South Fork cells are within A-R x 1-9.
##   3. Emerald Queen cells are within A-F x 1-9.
##   4. Location names are exactly:
##        "South Fork"
##        "Emerald Queen"


choose_grts_cells_by_location <- function(cell_table,
                                          n_by_location,
                                          cell_col = "cell",
                                          veg_col = "vegetation",
                                          location_col = "location",
                                          stratify_by_vegetation = TRUE,
                                          vegetation_allocation = "proportional",
                                          seed = NULL,
                                          n_over = 0) {

  # Optional seed makes the random draw reproducible (don't leave this
  # the same between days or it will draw the same cells!!).
  if (!is.null(seed)) set.seed(seed)

  # Convert input to a plain data frame.
  dat <- as.data.frame(cell_table)

  # Check for required columns.
  required_cols <- c(cell_col, location_col)

  if (stratify_by_vegetation) {
    required_cols <- c(required_cols, veg_col)
  }

  missing_cols <- setdiff(required_cols, names(dat))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # n_by_location must have names so the function knows which n
  # belongs to which location.
  if (is.null(names(n_by_location))) {
    stop(
      "n_by_location must be a named vector, e.g. ",
      'c("South Fork" = 10, "Emerald Queen" = 6).'
    )
  }

  # Parse grid labels into numeric x/y-like coordinates.
  #
  # This is the only grid-label validation retained.
  # Labels must be parseable as number-letter or letter-number, e.g.
  #   1A, 2B, 9R
  # or
  #   A1, B2, R9
  #
  # The function no longer checks whether the labels fall inside a
  # pre-defined South Fork or Emerald Queen grid extent.
  parsed <- parse_grid_label(dat[[cell_col]])

  dat$grid_col_number <- parsed$grid_col_number
  dat$grid_row_letter <- parsed$grid_row_letter
  dat$grid_row_number <- parsed$grid_row_number

  # Create artificial spatial coordinates from the grid. The numbered
  # part of the cell label becomes x.  The lettered part of the cell
  # label becomes y.

  # These coordinates are only used to make the GRTS draw spatially balanced.
  dat$grts_x <- (dat$grid_col_number - 0.5) 
  dat$grts_y <- (dat$grid_row_number - 0.5) 

  # Storage objects for the outputs from each location.
  selected_all <- list()
  replacement_all <- list()
  allocation_all <- list()
  design_all <- list()

  # Loop over locations and run a separate GRTS draw in each one.
  for (loc in names(n_by_location)) {

    # Number of cells requested for this location.
    n_loc <- as.integer(n_by_location[[loc]])

    # Subset to one location.
    dat_loc <- dat[dat[[location_col]] == loc, , drop = FALSE]

    if (n_loc < 1) {
      stop("Requested sample size must be at least 1 for ", loc, ".")
    }

    if (n_loc > nrow(dat_loc)) {
      stop(
        "Requested ", n_loc, " cells from ", loc,
        ", but only ", nrow(dat_loc), " cells are available."
      )
    }

    # Convert the location-specific data frame into an sf object.
    # spsurvey::grts() expects a spatial sampling frame.
    sframe_loc <- sf::st_as_sf(
      dat_loc,
      coords = c("grts_x", "grts_y"),
      crs = NA,
      remove = FALSE
    )

    if (stratify_by_vegetation) {

      # Decide how many cells to sample from each vegetation category
      # within this location.
      n_by_veg <- allocate_within_location(
        strata = dat_loc[[veg_col]],
        n = n_loc,
        allocation = vegetation_allocation
      )

      # Remove vegetation categories assigned zero samples.
      n_by_veg <- n_by_veg[n_by_veg > 0]

      # Keep only cells belonging to vegetation categories being sampled.
      sframe_loc <- sframe_loc[
        as.character(sframe_loc[[veg_col]]) %in% names(n_by_veg),
      ]

      # Arguments for a stratified GRTS draw.
      grts_args <- list(
        sframe = sframe_loc,
        n_base = n_by_veg,
        stratum_var = veg_col,
        DesignID = loc,
        sep = "-",
        projcrs_check = FALSE
      )

      allocation_all[[loc]] <- n_by_veg

    } else {

      # Arguments for an unstratified GRTS draw.
      grts_args <- list(
        sframe = sframe_loc,
        n_base = n_loc,
        DesignID = loc,
        sep = "-",
        projcrs_check = FALSE
      )

      allocation_all[[loc]] <- n_loc
    }

    # Optional oversample/replacement cells.
    #
    # These are useful if a selected cell cannot be surveyed because it is
    # inaccessible, unsafe, underwater, misclassified, etc.
    if (!is.null(n_over) && n_over > 0) {
      grts_args$n_over <- n_over
    }

    # Run the GRTS draw for this location.
    design_loc <- do.call(spsurvey::grts, grts_args)

    # Extract the primary selected sites and remove the sf geometry column.
    selected <- sf::st_drop_geometry(design_loc$sites_base)

    # Add an explicit order within location.
    selected$survey_order_within_location <- seq_len(nrow(selected))

    # Extract replacement sites, if requested.
    replacements <- NULL

    if (!is.null(design_loc$sites_over)) {
      replacements <- sf::st_drop_geometry(design_loc$sites_over)
      replacements$replacement_order_within_location <- seq_len(nrow(replacements))
    }

    # Store outputs for this location.
    selected_all[[loc]] <- selected
    replacement_all[[loc]] <- replacements
    design_all[[loc]] <- design_loc
  }

  # Combine selected cells from all locations into one data frame.
  selected_all <- do.call(rbind, selected_all)

  # Add an overall order across both locations.
  selected_all$survey_order_overall <- seq_len(nrow(selected_all))

  # Combine replacement cells, if any were created.
  replacements_nonnull <- replacement_all[
    !vapply(replacement_all, is.null, logical(1))
  ]

  replacements_all <- if (length(replacements_nonnull) > 0) {
    out <- do.call(rbind, replacements_nonnull)
    out$replacement_order_overall <- seq_len(nrow(out))
    out
  } else {
    NULL
  }

  # Return all useful pieces.
  list(
    selected = selected_all,
    replacements = replacements_all,
    allocation = allocation_all,
    design = design_all
  )
}


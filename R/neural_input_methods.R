#' @rdname neural_input
#' @export
neural_input.Reg <- function(x, start = 0, end = NULL, resolution = 0.33, ...) {
  # Determine end time if not provided
  if (is.null(end)) {
    end <- max(x$onsets + x$duration) + 10  # Add 10 seconds buffer
  }
  
  # Create time vector
  time <- seq(from = start, to = end, by = resolution)
  
  # Initialize neural input vector
  neural <- numeric(length(time))
  
  # For each event, add the amplitude during its duration
  for (i in seq_along(x$onsets)) {
    event_start <- x$onsets[i]
    event_end <- event_start + x$duration[i]
    
    # Find time indices within event duration
    event_indices <- which(time >= event_start & time < event_end)
    
    # Add amplitude for this event
    neural[event_indices] <- neural[event_indices] + x$amplitude[i]
  }
  
  # Return list with time and neural input
  list(time = time, neural_input = neural)
}
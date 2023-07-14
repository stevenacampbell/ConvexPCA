
transparent_col <- function(color_name, p) {
  # Adjusts color transparency
  
  # color_name: color name
  # p: transparency percentage
  
  # Get RGB values
  col <- col2rgb(color_name)
  
  # Apply transparency
  new_col <- rgb(col[1], col[2], col[3], maxColorValue = 255, alpha = 255*(1-p))

  return(new_col)
}

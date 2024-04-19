#' custum colors
cb_palette_8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_palette <- c("#193cbcff", "#1473afff", "#589acfff", "#89c3efff", "#ea594eff", "#e5b039ff", "#ede65aff")
cb_palette_4 <- c(viridis(4)[1:3], cb_palette[6])

get_viridis_color <- function(x, opt="D", begin=0, end=1, reverse=FALSE) {
  x <- x * (end - begin) + begin
  cmap <- viridisLite::viridis.map[viridisLite::viridis.map$opt == opt,]
  if (reverse) cmap <- cmap[rev(seq_len(nrow(cmap))),]
  map_rgbs <- grDevices::rgb(cmap$R, cmap$G, cmap$B)
  ramp <- grDevices::colorRamp(map_rgbs, space="Lab", interpolate="spline")
  out_rgbs <- ramp(x) / 255
  grDevices::rgb(out_rgbs[,1], out_rgbs[,2], out_rgbs[,3])
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%.%10^", l)
  l <- gsub("[+]", "", l)
  # return this as an expression
  parse(text=l)
}

# scientific <- function(x){
#   ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
# }

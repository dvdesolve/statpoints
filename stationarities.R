### stationarities.R
##
## This script finds all stationary points for function of two variables given as grid (e. g., file with "x" "y" "z" columns)
## Grid bin size (dx, dy) can be arbitrary and vary along any given range. NB: for step-fixed grids numerical differentiation can be performed in much more efficient way

# stationary points classifier
stat_type <- function(dx2, dy2, dxdy) {
    if (is.na(dx2) | is.na(dy2) | is.na(dxdy)) return("unk")
    if (dx2 * dy2 - dxdy^2 < 0) return("saddle")
    if (dx2 * dy2 - dxdy^2 > 0) {
        if (dx2 < 0 & dy2 < 0) return("max")
        else if (dx2 > 0 & dy2 > 0) return("min")
        else return("unk")
    }
    else return("unk")
}

# test function (can be replaced with xyz dataframe below in script)
z_func <- function(x, y) {
    return(2 * x^3 + 6 * x * y^2 - 3 * y^3 - 150 * x)
}

# partial first- and second-order analytical derivatives of test function (just for comparison and debugging purposes)
z_func_dx <- function(x, y) {
    return(6 * x^2 + 6 * y^2 - 150)
}

z_func_dy <- function(x, y) {
    return(12 * x * y - 9 * y^2)
}

z_func_dx2 <- function(x, y) {
    return(12 * x)
}

z_func_dy2 <- function(x, y) {
    return(12 * x - 18 * y)
}

z_func_dxdy <- function(x, y) {
    return(12 * y)
}

# partial first- and second-order numerical derivatives of user function/dataframe
# if grid sizes along x and y axes doesn't change then it's better to calculate dx and dy before differentiation and use them directly in points selection
z_dx <- function(x, y) {
    pos_x <- which(df_x == x)
    
    if (pos_x == 1 | pos_x == length(df_x)) return(NA)
    
    xf <- df_x[pos_x + 1]
    xb <- df_x[pos_x - 1]

    return((df[which(df$x == xf & df$y == y), "z"] - df[which(df$x == xb & df$y == y), "z"]) / (xf - xb))
}

z_dy <- function(x, y) {
    pos_y <- which(df_y == y)
    
    if (pos_y == 1 | pos_y == length(df_y)) return(NA)
    
    yf <- df_y[pos_y + 1]
    yb <- df_y[pos_y - 1]
    
    return((df[which(df$y == yf & df$x == x), "z"] - df[which(df$y == yb & df$x == x), "z"]) / (yf - yb))
}

z_dx2 <- function(x, y) {
    pos_x <- which(df_x == x)
    
    if (pos_x <= 2 | pos_x >= (length(df_x) - 1)) return(NA)
    
    xf <- df_x[pos_x + 1]
    xb <- df_x[pos_x - 1]
    
    return((df[which(df$x == xf & df$y == y), "dx"] - df[which(df$x == xb & df$y == y), "dx"]) / (xf - xb))
}

z_dy2 <- function(x, y) {
    pos_y <- which(df_y == y)
    
    if (pos_y <= 2 | pos_y >= (length(df_y) - 1)) return(NA)
    
    yf <- df_y[pos_y + 1]
    yb <- df_y[pos_y - 1]
    
    return((df[which(df$y == yf & df$x == x), "dy"] - df[which(df$y == yb & df$x == x), "dy"]) / (yf - yb))
}

z_dxdy <- function(x, y) {
    pos_y <- which(df_y == y)
    
    if (pos_y <= 2 | pos_y >= (length(df_y) - 1)) return(NA)
    
    yf <- df_y[pos_y + 1]
    yb <- df_y[pos_y - 1]
    
    return((df[which(df$y == yf & df$x == x), "dx"] - df[which(df$y == yb & df$x == x), "dx"]) / (yf - yb))
}

# just for test purposes; one can read necessary dataframe directly instead of two executing following lines
df <- expand.grid("x" = seq(-10, 10, 0.2), "y" = seq(-10, 10, 0.2))
df$z <- apply(df, 1, function(x) z_func(x["x"], x["y"]))

df$dx_a <- apply(df, 1, function(x) z_func_dx(x["x"], x["y"]))
df$dy_a <- apply(df, 1, function(x) z_func_dy(x["x"], x["y"]))
df$dx2_a <- apply(df, 1, function(x) z_func_dx2(x["x"], x["y"]))
df$dy2_a <- apply(df, 1, function(x) z_func_dy2(x["x"], x["y"]))
df$dxdy_a <- apply(df, 1, function(x) z_func_dxdy(x["x"], x["y"]))

# store all unique abscissae and ordinates
df_x <- unique(df$x)
df_y <- unique(df$y)

# calculate first-order numerical derivatives
df$dx <- apply(df, 1, function(x) z_dx(x["x"], x["y"]))
df$dy <- apply(df, 1, function(x) z_dy(x["x"], x["y"]))

# one can adjust tolerance of selection criteria (in units of [dependent variable]/[independent variable])
tolerance_x <- 1.0
tolerance_y <- 1.0

# filter out any non-stationary points
df$stat <- apply(df, 1, function(x) ifelse(is.na(x["dx"]) | is.na(x["dy"]), FALSE, ifelse(abs(x["dx"]) <= tolerance_x & abs(x["dy"]) <= tolerance_y, TRUE, FALSE)))
stationaries <- df[which(df$stat == TRUE), ]

# calculate second-order numerical derivatives
stationaries$dx2 <- apply(stationaries, 1, function(x) z_dx2(x["x"], x["y"]))
stationaries$dy2 <- apply(stationaries, 1, function(x) z_dy2(x["x"], x["y"]))
stationaries$dxdy <- apply(stationaries, 1, function(x) z_dxdy(x["x"], x["y"]))

# classify stationary points into main classes
stationaries$type <- apply(stationaries, 1, function(x) stat_type(x["dx2"], x["dy2"], x["dxdy"]))
minima <- stationaries[which(stationaries$type == "min"), ]
maxima <- stationaries[which(stationaries$type == "max"), ]
saddles <- stationaries[which(stationaries$type == "saddle"), ]
unks <- stationaries[which(stationaries$type == "unk"), ]

# plot function and stationary points
library(ggplot2)
ggplot(df, aes(x, y, z = z)) + stat_contour(geom = "polygon", aes(fill = ..level..)) +
    geom_raster(aes(fill = z), interpolate = TRUE) + stat_contour(binwidth = 100.0, size = 0.075, colour = "black") + xlab("x") + ylab("y") +
    guides(fill = guide_colorbar(title = "z", barwidth = 15)) + scale_fill_gradientn(colours = c("blue", "steelblue3", "green3", "yellow3", "orange1", "firebrick3")) + theme_bw() + theme(legend.position = "bottom", axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_point(inherit.aes = FALSE, data = stationaries, aes(x, y, colour = type), pch = 20, size = 3.5) + scale_colour_manual(labels = c("max" = "Max", "min" = "Min", "saddle" = "Saddle", "unk" = "Unknown"), values = c("max" = "orangered", "min" = "olivedrab2", "saddle" = "purple1", "unk" = "black")) + guides(colour = guide_legend(title = "Type"))

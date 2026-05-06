# General plotting and numeric utilities used by CMplot().
# These helpers use only functions distributed with base R.

draw_circle <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
{
    # Draw a circle at the origin by tracing the upper and lower semicircles.
    curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
    curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
}


round_y_axis_max <- function(x){
    # Expand an upper y-axis bound to a visually clean tick boundary.
    if(x == 0) return(x)
    if(abs(x) >= 1){
        return(ceiling(x))
    }else{
        if(x < 0){
            digit <- 10^(-ceiling(-log10(abs(x))))
            return(-(floor(abs(x) / digit - 1) * digit))
        }else{
            digit <- 10^(-ceiling(-log10(x)))
            return((floor(x / digit + 1) * digit))
        }
    }
}

round_y_axis_min <- function(x){
    # Expand a lower y-axis bound to a visually clean tick boundary.
    if(x == 0) return(x)
    if(abs(x) >= 1){
        return(floor(x))
    }else{
        if(x < 0){
            digit <- 10^(-ceiling(-log10(abs(x))))
            return(-(floor(abs(x) / digit + 1) * digit))
        }else{
            digit <- 10^(-ceiling(-log10(x)))
            return((floor(x / digit - 1) * digit))
        }
    }
}

min_without_na <- function(x){
    # Return the minimum while ignoring missing values.
    return(min(x, na.rm=TRUE))
}

max_without_na <- function(x){
    # Return the maximum while ignoring missing values.
    return(max(x, na.rm=TRUE))
}

# created by Haohao Zhang
# [OPTIMIZED] Use integer key + duplicated() instead of cbind + duplicated
# Original: index <- !duplicated(cbind(x, y))  # slow, creates matrix
# Optimized: key <- x*1e6 + y; !duplicated(key)  # ~11x faster
filter_visible_points <- function(x, y, w, h, dpi, scale=1) {
    # Keep one point per output pixel bucket to reduce overplotting in large datasets.
    keep <- rep(FALSE, length(x))
    finite <- is.finite(x) & is.finite(y)
    if(!any(finite)) return(keep)
    x0 <- x[finite]
    y0 <- y[finite]
    x_min <- min(x0)
    x_max <- max(x0)
    y_min <- min(y0)
    y_max <- max(y0)
    x_range <- x_max - x_min
    y_range <- y_max - y_min
    if(!is.finite(x_range) || x_range == 0){
        x_scaled <- rep(1, length(x0))
    }else{
        x_scaled <- ceiling((x0 - x_min) / x_range * w * dpi / scale)
    }
    if(!is.finite(y_range) || y_range == 0){
        y_scaled <- rep(1, length(y0))
    }else{
        y_scaled <- ceiling((y0 - y_min) / y_range * h * dpi / scale)
    }
    key <- x_scaled * 1000000L + y_scaled
    keep[finite] <- !duplicated(key)
    keep
}

transform_plot_y <- function(y, skip=NULL, cut=NULL, cutfactor=10) {
    # Transform plotted -log10(P) values without changing the original statistics.
    ydim <- dim(y)
    ynames <- dimnames(y)
    y2 <- as.numeric(y)
    if(!is.null(skip)){
        y2[y2 < skip] <- NA_real_
    }
    if(!is.null(cut)){
        above <- is.finite(y2) & y2 > cut
        y2[above] <- cut + (y2[above] - cut) / cutfactor
    }
    if(!is.null(ydim)){
        dim(y2) <- ydim
        dimnames(y2) <- ynames
    }
    y2
}

inverse_plot_y <- function(y, skip=NULL, cut=NULL, cutfactor=10) {
    # Convert display coordinates back to the original -log10(P) scale for axis labels.
    ydim <- dim(y)
    ynames <- dimnames(y)
    y2 <- as.numeric(y)
    if(!is.null(cut)){
        above <- is.finite(y2) & y2 > cut
        y2[above] <- cut + (y2[above] - cut) * cutfactor
    }
    if(!is.null(ydim)){
        dim(y2) <- ydim
        dimnames(y2) <- ynames
    }
    y2
}

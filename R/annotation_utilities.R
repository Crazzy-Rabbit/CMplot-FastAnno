# Annotation and highlighted-label utilities used by CMplot().
# These functions parse user annotation tables and place labels on Manhattan plots.

# [FEATURE] Centralized highlighted label renderer for the new top/nearby annotation modes.
draw_highlight_labels <- function(
    x,
    y,
    words=NULL,
    point.cex=1,
    text.cex=1,
    pch=19,
    type = "p",
    point.col = "red",
    text.col = "black",
    text.font=3,
    top.cex=1,
    mode="top",
    side="auto",
    nearby.offset=0.012,
    nearby.step=0.04,
    line.col="grey85",
    line.lwd=1,
    line.lty=1,
    line.mode="auto",
    line.bend=0.08,
    arrow=TRUE,
    arrow.length=0.06,
    top.space=0.06,
    top.inside=FALSE,
    optimize=TRUE,
    lanes=1,
    lane.gap=0.055,
    min.gap=0.004,
    xlim=c(-Inf, Inf),
    ylim=c(-Inf, Inf)
)
{
    # Draw highlighted SNP points and optional labels in scatter, top, or nearby modes.
    overlap <- function(x1, y1, sw1, sh1, boxes) {
        if (length(boxes) == 0) return(FALSE)
        for (i in c(1:length(boxes))) {
            bnds <- boxes[[i]]
            x2 <- bnds[1]
            y2 <- bnds[2]
            sw2 <- bnds[3]
            sh2 <- bnds[4]

            if (x1 < x2)
                overlap <- x1 + sw1 > x2
            else
                overlap <- x2 + sw2 > x1

            if (y1 < y2)
                overlap <- overlap && (y1 + sh1 > y2)
            else
                overlap <- overlap && (y2 + sh2 > y1)

            if (overlap) {
                return(TRUE)
            }
        }
        return(FALSE)
    }

    layout <- function(x, y, words, cex=1, xlim=c(-Inf, Inf), ylim=c(-Inf, Inf)) {
        sdx <- sd(x, na.rm=TRUE)
        sdy <- sd(y, na.rm=TRUE)
        if (sdx == 0) sdx <- 1
        if (sdy == 0) sdy <- 1
        # [OPTIMIZED] Pre-compute strwidth/strheight once
        wid_all <- strwidth(words, cex=cex)
        ht_all <- strheight(words, cex=cex)
        boxes <- list()
        for(i in seq_along(words)){
            wid <- wid_all[i]
            ht <- ht_all[i]
            if(i <= (length(words) / 2)){
                boxes[[length(boxes) + 1]] <- c(x[i]-0.5*wid, y[i]-0.5*ht, wid, ht)
            }else{
                xupdt <- xrot <- x[i]
                yupdt <- yrot <- y[i]
                r <- 0
                theta <- runif(1, 0, 2 * pi)
                ht <- 1.5 * ht
                isOverlaped <- TRUE
                while(isOverlaped){
                    if(
                        !overlap(xupdt-0.5*wid, yupdt-0.5*ht, wid, ht, boxes) &&
                        (xupdt-0.5*wid) > xlim[1] &&
                        (yupdt-0.5*ht) > ylim[1] &&
                        (xupdt+0.5*wid) < xlim[2] &&
                        (yupdt+0.5*ht) < ylim[2]
                    ){
                        boxes[[length(boxes) + 1]] <- c(xupdt-0.5*wid, yupdt-0.5*ht, wid, ht)
                        isOverlaped <- FALSE
                    }else{
                        theta <- theta + 0.1
                        r <- r + 0.001 / (2 * base::pi)
                        xupdt <- xrot + 0.1 * sdx * r * cos(theta)
                        yupdt <- yrot + sdy * r * sin(theta)
                    }
                }
            }
        }
        result <- do.call(rbind, boxes)
        colnames(result) <- c("x", "y", "width", "ht")
        rownames(result) <- words
        result
    }

    draw_highlight_points <- function(x, y, pch, type, point.col, point.cex, ylim) {
        if(type=="h"){
            points(x,y,pch=pch,type="h",col=point.col, lwd=point.cex+1)
            points(x,y,pch=pch,type="p",col=point.col, cex=point.cex)
        }else if(type=="l"){
            segments(x, ylim[1], x, ylim[2], col=point.col, lwd=point.cex, lty=2)
        }else{
            points(x,y,pch=pch,type=type,col=point.col,cex=point.cex)
        }
    }

    # [FEATURE] Normalize user-facing side options and choose automatic left/right placement.
    normalize_label_side <- function(side, n, x=NULL, xlim=c(-Inf, Inf)) {
        side <- rep(side, length.out=n)
        side <- tolower(side)
        side[side %in% c("l", "left")] <- "left"
        side[side %in% c("r", "right")] <- "right"
        side[side %in% c("alt", "alternate")] <- "alternate"
        side[!side %in% c("auto", "left", "right", "alternate")] <- "auto"
        if(any(side == "alternate")){
            alt <- rep(c("left", "right"), length.out=sum(side == "alternate"))
            side[side == "alternate"] <- alt
        }
        if(!is.null(x) && any(side == "auto")){
            midpoint <- mean(xlim)
            if(!is.finite(midpoint)) midpoint <- stats::median(x, na.rm=TRUE)
            side[side == "auto"] <- ifelse(x[side == "auto"] > midpoint, "left", "right")
        }
        side
    }

    label_boxes_overlap <- function(box, boxes) {
        if(!length(boxes)) return(FALSE)
        for(b in boxes){
            if(box[1] < b[2] && box[2] > b[1] && box[3] < b[4] && box[4] > b[3]){
                return(TRUE)
            }
        }
        FALSE
    }

    # [FEATURE] Spread top labels in target order so connector arms do not cross unnecessarily.
    spread_label_x <- function(x, words, cex, xlim, anchor.x=NULL, min.gap=0.004) {
        n <- length(x)
        if(n <= 1L || any(!is.finite(xlim))) return(pmin(pmax(x, xlim[1]), xlim[2]))
        if(is.null(anchor.x)) anchor.x <- x
        cex <- rep(cex, length.out=n)
        half_width <- suppressWarnings(strheight(words, cex=cex) / 2)
        half_width[!is.finite(half_width)] <- 0
        xrange <- diff(xlim)
        auto_half_gap <- xrange / max(n * 2.4, 80)
        pad <- max(xrange * max(min.gap, 0), auto_half_gap, max(half_width, na.rm=TRUE) * 0.4, na.rm=TRUE)
        half_width <- half_width + pad

        ord <- order(anchor.x, x)
        pos <- x[ord]
        half_ord <- half_width[ord]
        left <- xlim[1] + half_ord
        right <- xlim[2] - half_ord
        if(any(left > right)) return(pmin(pmax(x, xlim[1]), xlim[2]))

        pos <- pmin(pmax(pos, left), right)
        for(i in 2:n){
            pos[i] <- max(pos[i], pos[i-1] + half_ord[i-1] + half_ord[i])
        }
        overflow <- pos[n] - right[n]
        if(is.finite(overflow) && overflow > 0) pos <- pos - overflow
        pos <- pmin(pmax(pos, left), right)
        for(i in (n-1):1){
            pos[i] <- min(pos[i], pos[i+1] - half_ord[i+1] - half_ord[i])
        }
        underflow <- left[1] - pos[1]
        if(is.finite(underflow) && underflow > 0) pos <- pos + underflow
        pos <- pmin(pmax(pos, left), right)

        out <- numeric(n)
        out[ord] <- pos
        out
    }

    # [FEATURE] Assign top labels to non-overlapping lanes while preserving x-order within each lane.
    assign_top_label_lanes <- function(anchor.x, preferred.x, words, cex, xlim, lanes=1L, min.gap=0.004, optimize=TRUE) {
        n <- length(anchor.x)
        lanes <- max(1L, min(as.integer(lanes[1]), max(1L, n)))
        if(!isTRUE(optimize) || n <= 1L) lanes <- 1L
        label.lane <- rep(1L, n)
        if(lanes > 1L && all(is.finite(xlim))){
            xrange <- diff(xlim)
            span <- suppressWarnings(strheight(words, cex=cex))
            span[!is.finite(span)] <- 0
            span <- span + max(xrange * max(min.gap, 0), 0)
            ord <- order(anchor.x, preferred.x)
            lane.right <- rep(-Inf, lanes)
            lane.count <- integer(lanes)
            for(ii in ord){
                costs <- rep(Inf, lanes)
                for(ln in seq_len(lanes)){
                    last.right <- lane.right[ln]
                    candidate <- preferred.x[ii]
                    if(is.finite(last.right)){
                        candidate <- max(candidate, last.right + span[ii])
                    }
                    overflow <- max(0, candidate + span[ii] / 2 - xlim[2])
                    costs[ln] <- abs(candidate - preferred.x[ii]) + overflow * 1000 +
                        lane.count[ln] * xrange * 0.0005 + (ln - 1L) * xrange * 0.0002
                }
                ln <- which.min(costs)
                label.lane[ii] <- ln
                lane.count[ln] <- lane.count[ln] + 1L
                lane.right[ln] <- max(lane.right[ln], preferred.x[ii] + span[ii] / 2, na.rm=TRUE)
            }
        }
        if(lanes > 1L){
            # [FEATURE] Cross-lane labels are still globally spaced on x because vertical labels
            # can visually collide even when their lane baselines differ.
            label.x <- spread_label_x(preferred.x, as.character(words), cex,
                xlim=xlim, anchor.x=anchor.x, min.gap=min.gap)
        }else{
            label.x <- preferred.x
            for(ln in sort(unique(label.lane))){
                idx <- which(label.lane == ln)
                label.x[idx] <- spread_label_x(preferred.x[idx], as.character(words[idx]), cex[idx],
                    xlim=xlim, anchor.x=anchor.x[idx], min.gap=min.gap)
            }
        }
        list(x=label.x, lane=label.lane, lanes=lanes)
    }

    if(!is.null(words)){
        if(length(x) != length(words))  stop("length of highlighted labels is not equal to the highlighted SNPs.")
        indx <- order(y, decreasing=TRUE)
        x <- x[indx]
        y <- y[indx]
        words <- words[indx]
        if(length(point.cex)!=1){if(length(point.cex)==length(x)){point.cex=point.cex[indx]}else{stop("unequal length of 'cex' for highlighted points.")}}else{point.cex=rep(point.cex,length(x))}
        if(length(pch)!=1){if(length(pch)==length(x)){pch=pch[indx]}else{stop("unequal length of 'pch' for highlighted points.")}}else{pch=rep(pch,length(x))}
        if(length(point.col)!=1){if(length(point.col)==length(x)){point.col=point.col[indx]}else{stop("unequal length of 'col' for highlighted points.")}}else{point.col=rep(point.col,length(x))}
        if(length(text.col)!=1){if(length(text.col)==length(x)){text.col=text.col[indx]}else{stop("unequal length of 'col' for highlighted text.")}}else{text.col=rep(text.col,length(x))}
        if(length(text.cex)!=1){if(length(text.cex)==length(x)){text.cex=text.cex[indx]}else{stop("unequal length of 'cex' for highlighted text.")}}else{text.cex=rep(text.cex,length(x))}
        
        words_ety <- words[words == "" | is.na(words)]
        if(length(words_ety)){
            logical_idx <- words == "" | is.na(words)
            if(type=="h"){
                points(x[logical_idx],y[logical_idx],pch=pch[logical_idx],type="h",col=point.col[logical_idx], lwd=point.cex[logical_idx]+1)
                points(x[logical_idx],y[logical_idx],pch=pch[logical_idx],type="p",col=point.col[logical_idx], cex=point.cex[logical_idx])
            }else if(type=="l"){
                segments(x[logical_idx], ylim[1], x[logical_idx], ylim[2], col=point.col[logical_idx], lwd=point.cex[logical_idx], lty=2)
            }else{
                points(x[logical_idx],y[logical_idx],pch=pch[logical_idx],type="p",col=point.col[logical_idx],cex=point.cex[logical_idx])
            }
            words <- words[!logical_idx]
            x <- x[!logical_idx]
            y <- y[!logical_idx]
            point.cex <- point.cex[!logical_idx]
            pch <- pch[!logical_idx]
            point.col <- point.col[!logical_idx]
            text.col <- text.col[!logical_idx]
            text.cex <- text.cex[!logical_idx]
        }

        if(length(words) == 0) return(invisible(NULL))

        if(mode == "nearby"){
            # [FEATURE] GWASLab-inspired side annotation: labels are placed next to target points.
            x1 <- x
            y1 <- y
            yrange <- diff(ylim)
            xrange <- diff(xlim)
            if(!is.finite(yrange) || yrange <= 0) yrange <- 1
            if(!is.finite(xrange) || xrange <= 0) xrange <- 1
            nearby.offset <- max(nearby.offset, 0)
            nearby.step <- max(nearby.step, 0)
            side.vec <- normalize_label_side(side, length(x1), x=x1, xlim=xlim)
            dir <- ifelse(side.vec == "left", -1, 1)
            label.x <- x1 + dir * xrange * nearby.offset
            word.width <- strwidth(as.character(words), cex=text.cex)
            word.height <- strheight(as.character(words), cex=text.cex)
            left.edge <- ifelse(dir > 0, label.x, label.x - word.width)
            right.edge <- ifelse(dir > 0, label.x + word.width, label.x)
            flip <- (left.edge < xlim[1] | right.edge > xlim[2]) & rep(tolower(side), length.out=length(x1)) %in% c("auto", "alternate")
            if(any(flip)){
                dir[flip] <- -dir[flip]
                label.x[flip] <- x1[flip] + dir[flip] * xrange * nearby.offset
            }
            label.x <- pmin(pmax(label.x, xlim[1]), xlim[2])
            label.y <- y1 + yrange * nearby.step
            boxes <- list()
            ord <- order(y1, decreasing=TRUE)
            for(ii in ord){
                h <- max(word.height[ii], yrange * 0.015, na.rm=TRUE)
                w <- word.width[ii]
                try.y <- label.y[ii]
                for(step.i in 0:80){
                    left <- if(dir[ii] > 0) label.x[ii] else label.x[ii] - w
                    right <- if(dir[ii] > 0) label.x[ii] + w else label.x[ii]
                    box <- c(left, right, try.y - h * 0.6, try.y + h * 0.6)
                    if(!label_boxes_overlap(box, boxes)) break
                    try.y <- try.y + h * 1.35
                }
                label.y[ii] <- try.y
                boxes[[length(boxes) + 1L]] <- box
            }
            if(line.mode != "none"){
                segments(x1, y1, label.x, label.y,
                    col=rep(line.col, length.out=length(x1)),
                    lwd=rep(line.lwd, length.out=length(x1)),
                    lty=rep(line.lty, length.out=length(x1)), xpd=NA)
            }
            draw_highlight_points(x1, y1, pch, type, point.col, point.cex, ylim)
            text(label.x, label.y, as.character(words),
                adj=c(ifelse(dir > 0, 0, 1), 0.5), xpd=NA,
                cex=text.cex, col=text.col, font=text.font)
            return(invisible(NULL))
        }

        if(mode == "top"){
            # [FEATURE] Top annotation: keep diagonal arms above the plotting area, then drop vertical lines to targets.
            x1 <- x
            y1 <- y
            top.inside <- isTRUE(top.inside)
            text.cex <- text.cex * top.cex
            yrange <- diff(ylim)
            if(!is.finite(yrange) || yrange <= 0) yrange <- 1
            xrange <- diff(xlim)
            if(!is.finite(xrange) || xrange <= 0) xrange <- 1
            top.space <- max(top.space, 0)
            lane.gap <- max(lane.gap, 0)
            min.gap <- max(min.gap, 0)
            lanes <- max(1L, as.integer(lanes[1]))
            side.vec <- normalize_label_side(side, length(x1), x=x1, xlim=xlim)
            preferred.x <- x1
            preferred.x[side.vec == "left"] <- x1[side.vec == "left"] - xrange * 0.035
            preferred.x[side.vec == "right"] <- x1[side.vec == "right"] + xrange * 0.035
            if(any(!is.finite(xlim))) preferred.x <- x1
            lane.layout <- assign_top_label_lanes(x1, preferred.x, as.character(words), text.cex,
                xlim=xlim, lanes=lanes, min.gap=min.gap, optimize=optimize)
            label.x <- lane.layout$x
            label.lane <- lane.layout$lane
            lane.offset <- (label.lane - 1L) * lane.gap
            if(top.inside){
                line.top <- ylim[2] - yrange * (top.space * 0.65 + lane.offset)
                text.y <- ylim[2] - yrange * (top.space * 0.2 + lane.offset)
            }else{
                line.top <- ylim[2] + yrange * (top.space * 0.9 + lane.offset)
                text.y <- ylim[2] + yrange * (top.space * 1.45 + lane.offset)
            }
            line.col <- rep(line.col, length.out=length(x1))
            line.lwd <- rep(line.lwd, length.out=length(x1))
            line.lty <- rep(line.lty, length.out=length(x1))
            line.bend <- max(line.bend, 0)
            if(top.inside){
                arm.y <- pmin(line.top, pmax(y1, ylim[2] - yrange * (top.space * 1.2 + lane.offset * 0.35)))
            }else{
                arm.y <- ylim[2] + yrange * (top.space * 0.25 + lane.offset * 0.35)
                arm.y <- pmin(line.top, pmax(ylim[2], arm.y))
            }
            straight.tol <- xrange * 0.0025
            if(!is.finite(straight.tol) || straight.tol <= 0) straight.tol <- 0
            use.straight <- abs(label.x - x1) <= straight.tol
            if(line.mode == "straight") use.straight <- rep(TRUE, length(x1))
            if(line.mode == "elbow") use.straight <- rep(FALSE, length(x1))

            elbow.idx <- which(!use.straight)
            if(line.mode == "none"){
                # no connector lines
            }else if(isTRUE(arrow)){
                straight.idx <- which(use.straight)
                if(length(straight.idx)){
                    arrows(label.x[straight.idx], line.top[straight.idx],
                        x1[straight.idx], y1[straight.idx],
                        length=arrow.length, angle=15, code=2,
                        col=line.col[straight.idx], lwd=line.lwd[straight.idx],
                        lty=line.lty[straight.idx], xpd=TRUE)
                }
                if(length(elbow.idx)){
                    segments(label.x[elbow.idx], line.top[elbow.idx],
                        x1[elbow.idx], arm.y[elbow.idx],
                        col=line.col[elbow.idx], lwd=line.lwd[elbow.idx],
                        lty=line.lty[elbow.idx], xpd=TRUE)
                    arrows(x1[elbow.idx], arm.y[elbow.idx],
                        x1[elbow.idx], y1[elbow.idx],
                        length=arrow.length, angle=15, code=2,
                        col=line.col[elbow.idx], lwd=line.lwd[elbow.idx],
                        lty=line.lty[elbow.idx], xpd=TRUE)
                }
            }else{
                straight.idx <- which(use.straight)
                if(length(straight.idx)){
                    segments(label.x[straight.idx], line.top[straight.idx],
                        x1[straight.idx], y1[straight.idx],
                        col=line.col[straight.idx], lwd=line.lwd[straight.idx],
                        lty=line.lty[straight.idx], xpd=TRUE)
                }
                if(length(elbow.idx)){
                    segments(label.x[elbow.idx], line.top[elbow.idx],
                        x1[elbow.idx], arm.y[elbow.idx],
                        col=line.col[elbow.idx], lwd=line.lwd[elbow.idx],
                        lty=line.lty[elbow.idx], xpd=TRUE)
                    segments(x1[elbow.idx], arm.y[elbow.idx],
                        x1[elbow.idx], y1[elbow.idx],
                        col=line.col[elbow.idx], lwd=line.lwd[elbow.idx],
                        lty=line.lty[elbow.idx], xpd=TRUE)
                }
            }
            draw_highlight_points(x1, y1, pch, type, point.col, point.cex, ylim)
            text(label.x, text.y, as.character(words),
                srt=90, adj=c(ifelse(top.inside, 1, 0), 0.5), xpd=TRUE, cex=text.cex,
                col=text.col, font=text.font)
            return(invisible(NULL))
        }

        x1 <- x
        y1 <- y
        xadj <- sample(c(1.5, 0, -0.5), size=length(x), replace=TRUE)
        # xadj <- rep(c(1.5, 0, -0.5), length=max(3, length(x)))
        # xadj <- sort(xadj)[1:length(x)]
        # xadj[order(x)] <- xadj
        yadj <- rep(c(1.5, 0, -0.5), length=max(3, length(x)))
        yadj <- sort(yadj)[1:length(x)]
        # [OPTIMIZED] Pre-compute strwidth/strheight once before the loop
        word_widths <- strwidth(words, cex=text.cex)
        word_heights <- strheight(words, cex=text.cex)
        for(i in seq_along(x)){
            wht <- word_heights[i]
            wwd <- word_widths[i]
            if(xadj[i] == 0){
                if(yadj[i] == -0.5){
                    if((y[i] + 2*wht) > max(ylim)){
                        y[i] = y[i] - 1.5*wht
                    }else{
                        y[i] = y[i] + 1.5*wht
                    }
                }
                if(yadj[i] == 1.5)  y[i] = y[i] - 1.5*wht
            }else{
                if(yadj[i] == -0.5){
                    if((y[i] + 1.5*wht) > max(ylim)){
                        y[i] = y[i] - wht
                    }else{
                        y[i] = y[i] + wht
                    }
                }
                if(yadj[i] == -0.5) y[i] = y[i] + wht
                if(yadj[i] == 1.5)  y[i] = y[i] - wht
            }
            if(xadj[i] == 1.5){
                if((x[i] - 1.2*wwd) < min(xlim)){
                    x[i] = x[i] + 0.6*wwd
                }else{
                    x[i] = x[i] - 0.6*wwd
                }
            }
            if(xadj[i] == -0.5){
                if((x[i] + 1.2*wwd) > max(xlim)){
                    x[i] = x[i] - 0.6*wwd
                }else{
                    x[i] = x[i] + 0.6*wwd
                }
            }
        }

        x <- c(x1,x)
        y <- c(y1,y)
        words <- c(rep("OO", length(words)), as.character(words))
        lay <- layout(x=x,y=y,words=words,cex=c(rep(text.cex[1],length(x1)),text.cex),xlim=xlim,ylim=ylim)
        n <- length(x1)
        indd <- (n+1):length(x)
        for(i in indd){
            ii <- i - n
            xl <- lay[i,1]
            yl <- lay[i,2]
            w <- lay[i,3]
            h <- lay[i,4]
            nx <- xl + 0.5 * w
            ny <- yl + 0.5 * h
            if((nx + 0.5 * word_widths[ii]) < x1[ii]){
                nx=nx + 0.5 * word_widths[ii]
            }else if((nx - 0.5 * word_widths[ii]) > x1[ii]){
                nx=nx - 0.5 * word_widths[ii]
            }
            if((ny + word_heights[ii]) < y1[ii]){
                ny=ny + 0.5 * word_heights[ii]
            }else if((ny - word_heights[ii]) > y1[ii]){
                ny=ny - 0.5 * word_heights[ii]
            }
            # arrows(x1[i-n], y1[i-n], nx, ny, length=.08, angle=15, code=2, col="grey", lwd=2)
            segments(x1[ii], y1[ii], nx, ny, col="black", lwd=text.cex[ii])
        }
        if(type=="h"){
            points(x1,y1,pch=pch,type="h",col=point.col, lwd=point.cex+1)
            points(x1,y1,pch=pch,type="p",col=point.col, cex=point.cex)
        }else if(type=="l"){
            segments(x1, ylim[1], x1, ylim[2], col=point.col, lwd=point.cex, lty=2)
            # points(x1,y1,pch=pch,type="p",col=point.col, cex=point.cex)
        }else{
            points(x1,y1,pch=pch,type=type,col=point.col,cex=point.cex)
        }
        text(lay[indd,1]+0.5*lay[indd,3],lay[indd,2]+0.5*lay[indd,4],words[indd],xpd=TRUE,cex=text.cex,col=text.col,font=text.font)
    }else{
        if(type=="h"){
            points(x,y,pch=pch,type="h",col=point.col, lwd=point.cex+1)
            points(x,y,pch=pch,type="p",col=point.col, cex=point.cex)
        }else if(type=="l"){
            segments(x, ylim[1], x, ylim[2], col=point.col, lwd=point.cex, lty=2)
            # points(x,y,pch=pch,type="p",col=point.col, cex=point.cex)
        }else{
            points(x,y,pch=pch,type=type,col=point.col,cex=point.cex)
        }
    }
}

# [FEATURE] Read user-supplied annotation targets from CSV/TSV/data.frame.
read_annotation_targets <- function(annotation.file, sep=NULL) {
    # Read annotation targets from a data.frame, CSV file, or tab-delimited file.
    if(is.data.frame(annotation.file)) return(annotation.file)
    if(!is.character(annotation.file) || length(annotation.file) != 1L){
        stop("'annotation.file' should be a file path or a data.frame.")
    }
    if(!file.exists(annotation.file)) stop(paste("annotation.file does not exist:", annotation.file))
    if(is.null(sep)){
        ext <- tolower(sub("^.*\\.([^.]+)$", "\\1", basename(annotation.file)))
        if(identical(ext, basename(annotation.file))) ext <- ""
        sep <- ifelse(ext == "csv", ",", "\t")
    }
    utils::read.table(annotation.file, header=TRUE, sep=sep, stringsAsFactors=FALSE,
        check.names=FALSE, comment.char="", quote="\"", fill=TRUE)
}

# [FEATURE] Resolve annotation file column names while allowing common SNP/gene aliases.
resolve_annotation_column <- function(dat, requested=NULL, candidates=NULL, required=TRUE, fallback=NULL) {
    # Resolve a requested annotation column by index, exact name, or common aliases.
    nms <- names(dat)
    if(!length(nms)) stop("annotation.file must contain a header row.")
    if(!is.null(requested)){
        if(is.numeric(requested)){
            if(requested < 1 || requested > ncol(dat)) stop("annotation column index is out of range.")
            return(requested)
        }
        exact <- match(requested, nms)
        if(!is.na(exact)) return(exact)
        folded <- match(tolower(requested), tolower(nms))
        if(!is.na(folded)) return(folded)
        if(required) stop(paste("Cannot find annotation column:", requested))
        return(fallback)
    }
    if(length(candidates)){
        folded <- match(tolower(candidates), tolower(nms))
        folded <- folded[!is.na(folded)]
        if(length(folded)) return(folded[1])
    }
    if(required) stop(paste("Cannot find annotation column. Tried:", paste(candidates, collapse=", ")))
    fallback
}

# [FEATURE] Convert annotation table rows into CMplot highlight/highlight.text lists.
build_highlight_from_annotations <- function(annotation.file, trait, R, sep=NULL, snp.col=NULL, label.col=NULL, trait.col=NULL) {
    # Convert annotation rows into per-trait highlight and label lists.
    anno <- read_annotation_targets(annotation.file, sep=sep)
    snp.idx <- resolve_annotation_column(anno, snp.col,
        candidates=c("SNP", "SNPID", "Marker", "MarkerName", "ID", "rsID", "variant", "variant_id"),
        required=TRUE)
    label.idx <- resolve_annotation_column(anno, label.col,
        candidates=c("Label", "label", "Gene", "gene", "GeneName", "Annotation", "annotation", "Name", "name"),
        required=FALSE, fallback=snp.idx)
    trait.idx <- resolve_annotation_column(anno, trait.col,
        candidates=c("Trait", "trait", "Phenotype", "phenotype"),
        required=FALSE, fallback=NULL)

    snp <- as.character(anno[[snp.idx]])
    label <- as.character(anno[[label.idx]])
    label[is.na(label) | label == ""] <- snp[is.na(label) | label == ""]
    keep <- !is.na(snp) & snp != ""
    snp <- snp[keep]
    label <- label[keep]
    if(!length(snp)) stop("annotation.file does not contain any non-empty SNP IDs.")

    if(is.null(trait.idx)){
        if(R > 1){
            warning("annotation.file has no trait column; the same annotation rows are reused for every trait. Add a Trait column or set annotation.trait.col for per-trait annotations.", call.=FALSE)
        }
        highlight <- rep(list(snp), R)
        highlight.text <- rep(list(label), R)
    }else{
        anno.trait <- as.character(anno[[trait.idx]])[keep]
        highlight <- vector("list", R)
        highlight.text <- vector("list", R)
        for(i in seq_len(R)){
            use <- anno.trait %in% c(trait[i], as.character(i))
            highlight[[i]] <- snp[use]
            highlight.text[[i]] <- label[use]
        }
    }
    list(highlight=highlight, highlight.text=highlight.text)
}

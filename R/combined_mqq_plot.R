# Combined Manhattan + Q-Q plot branch used by CMplot().

draw_combined_mqq_plots <- function(env) {
    evalq({
        requested.mqq <- plot.type[plot.type %in% c("mqq", "qqm")]
        if(!length(requested.mqq)) return(invisible(NULL))
        if(multracks || multraits) stop("'mqq' and 'qqm' combined plots currently support one panel per trait; set multracks=FALSE and multraits=FALSE.")
        if(!is.null(file.name) && length(file.name) != R) stop(paste("please provide a vector containing file names of all", R, "traits."))

        open_combined_device <- function(prefix, name, wh, ht) {
            if(file=="jpg") jpeg(paste(prefix, ".", name, ".jpg", sep=""), width=wh*dpi, height=ht*dpi, res=dpi, quality=100)
            if(file=="pdf") pdf(paste(prefix, ".", name, ".pdf", sep=""), width=wh, height=ht)
            if(file=="tiff") tiff(paste(prefix, ".", name, ".tiff", sep=""), width=wh*dpi, height=ht*dpi, res=dpi)
            if(file=="png") png(paste(prefix, ".", name, ".png", sep=""), width=wh*dpi, height=ht*dpi, res=dpi, bg=NA)
        }

        draw_y_axis <- function(Min, Max, scale=1) {
            ats <- pretty(c(Min, Max), n=6)
            ats <- ats[is.finite(ats) & ats >= Min & ats <= Max]
            if(!length(ats)) ats <- c(Min, Max)
            labs <- inverse_y(ats)
            labs <- ifelse(abs(labs - round(labs)) < 1e-8, as.character(round(labs)), format(round(labs, 2), trim=TRUE))
            axis(2, at=ats, labels=labs, las=1, lwd=axis.lwd*scale, cex.axis=axis.cex*scale, font=lab.font)
            axis(2, at=c(Min, Max), labels=c("", ""), tcl=0, lwd=axis.lwd*scale)
        }

        draw_combined_manhattan_panel <- function(i, panel.wh, ht) {
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            logpvalue <- logpvalueT[,i]
            finite.y <- logpvalue[is.finite(logpvalue)]
            if(!length(finite.y)) stop("No finite values remain for Manhattan plotting after skip/cut filtering.")
            y.candidates <- finite.y
            if(!is.null(threshold) && !is.null(threshold[[i]])) y.candidates <- c(y.candidates, threshold_to_y(threshold[[i]]))
            if(is.null(ylim)){
                Max <- round_y_axis_max(max_without_na(y.candidates))
                Min <- round_y_axis_min(min_without_na(y.candidates))
                if(!is.null(skip)) Min <- max(Min, skip)
                if(Max <= Min) Max <- Min + 1
                plot.ylim <- c(Min, Max)
            }else{
                plot.ylim <- plot_y(ylim[[i]])
                Min <- min_without_na(plot.ylim)
                Max <- max_without_na(plot.ylim)
            }

            par(mar=highlight_mar(c(mar[1], mar[2], mar[3], 1)), xaxs="i", yaxs="i")
            plot(NULL, xlim=c(min_without_na(pvalue.posN)-band, band+max_without_na(pvalue.posN)),
                ylim=plot.ylim, ann=FALSE, axes=FALSE)
            visible <- filter_visible_points(pvalue.posN, logpvalue, panel.wh, ht, dpi=dpi)
            visible <- visible & is.finite(logpvalue)
            point.col <- rep(rep(colx, N[i]), add[[i]])
            points(pvalue.posN[visible], logpvalue[visible], pch=pch[i], type=type,
                lwd=cex[2]+1, cex=cex[2], col=point.col[visible])
            if(chr.border){
                for(b in seq_along(chr.border.pos)){
                    segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=axis.lwd, lty=2)
                }
            }
            if(!is.null(threshold) && !is.null(threshold[[i]])){
                for(thr in seq_along(threshold[[i]])){
                    h <- threshold_to_y(threshold[[i]][thr])
                    if(is.finite(h)) segments(0, h, max_without_na(pvalue.posN), h,
                        col=threshold.col[thr], lwd=threshold.lwd[thr], lty=threshold.lty[thr])
                }
            }
            draw_rect_chr_axis(scale=1, Min=Min, Max=Max, title=TRUE)
            draw_y_axis(Min, Max)
            mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)

            if(!is.null(highlight) && !is.na(highlight_index[[i]][1])){
                hcol <- if(is.null(highlight.col)) point.col[highlight_index[[i]]] else highlight_col[[i]]
                draw_highlight_labels(
                    x=pvalue.posN[highlight_index[[i]]],
                    y=logpvalue[highlight_index[[i]]],
                    xlim=c(min_without_na(pvalue.posN)-band, max_without_na(pvalue.posN)+band),
                    ylim=c(Min, Max),
                    words=highlight.text[[i]],
                    point.cex=highlight.cex,
                    text.cex=highlight.text.cex,
                    pch=highlight.pch,
                    type=highlight.type,
                    point.col=hcol,
                    text.col=highlight.text.col,
                    text.font=highlight.text.font,
                    top.cex=highlight.text.top.cex,
                    mode=highlight.text.mode,
                    side=highlight.text.side,
                    nearby.offset=highlight.text.nearby.offset,
                    nearby.step=highlight.text.nearby.step,
                    line.col=highlight.text.line.col,
                    line.lwd=highlight.text.line.lwd,
                    line.lty=highlight.text.line.lty,
                    line.mode=highlight.text.line.mode,
                    line.bend=highlight.text.line.bend,
                    arrow=highlight.text.arrow,
                    arrow.length=highlight.text.arrow.length,
                    top.space=highlight.text.top.space,
                    top.inside=highlight.text.top.inside,
                    optimize=highlight.text.optimize,
                    lanes=highlight.text.lanes,
                    lane.gap=highlight.text.lane.gap,
                    min.gap=highlight.text.min.gap)
            }
            if(box) box(lwd=axis.lwd)
        }

        draw_combined_qq_panel <- function(i, panel.wh, ht) {
            qdat <- prepare_qq_data(i)
            log.Quantiles <- qdat$log.Quantiles
            log.P.values <- qdat$log.P.values
            N.qq <- qdat$N
            if(N.qq == 0) stop("No finite values remain for Q-Q plotting.")
            if(conf.int){
                xi <- ceiling((10^-log.Quantiles) * N.qq)
                xi[xi == 0] <- 1
                c95 <- qbeta(0.95, xi, N.qq - xi + 1)
                c05 <- qbeta(0.05, xi, N.qq - xi + 1)
                ci.low <- plot_y(-log10(c05))
                ci.high <- plot_y(-log10(c95))
                index <- length(c95):1
            }else{
                ci.low <- ci.high <- NA_real_
            }
            y.candidates <- c(log.P.values, ci.low, ci.high)
            if(is.null(ylim)){
                YlimMax <- round_y_axis_max(max_without_na(y.candidates[is.finite(y.candidates)]))
                YlimMin <- if(!is.null(skip)) skip else 0
                if(YlimMax <= YlimMin) YlimMax <- YlimMin + 1
                plot.ylim <- c(YlimMin, YlimMax)
            }else{
                plot.ylim <- plot_y(ylim[[i]])
                YlimMin <- min_without_na(plot.ylim)
                YlimMax <- max_without_na(plot.ylim)
            }
            x.max <- floor(max_without_na(log.Quantiles)+1)
            par(mar=c(mar[1], max(4.2, mar[2]), mar[3], mar[4]), xpd=TRUE)
            plot(NULL, xlim=c(0, x.max), ylim=plot.ylim, axes=FALSE, xlab="", ylab="")
            axis(1, mgp=c(3,xticks.pos,0), at=pretty(c(0, x.max), n=6),
                lwd=axis.lwd, cex.axis=axis.cex)
            draw_y_axis(YlimMin, YlimMax)
            mtext(side=1, text=expression(Expected~~-log[10](italic(p))), line=ylab.pos+1,
                cex=lab.cex, font=lab.font, xpd=TRUE)
            mtext(side=2, text=expression(Observed~~-log[10](italic(p))), line=ylab.pos,
                cex=lab.cex, font=lab.font, xpd=TRUE)
            if(conf.int){
                qq.col <- if(is.null(conf.int.col)) t(col)[i] else conf.int.col[i]
                polygon(c(log.Quantiles[index], log.Quantiles),
                    c(ci.low[index], ci.high),
                    col=rgb(t(col2rgb(qq.col)), alpha=points.alpha, maxColorValue=255),
                    border=rgb(t(col2rgb(qq.col)), alpha=points.alpha, maxColorValue=255))
            }
            if(!is.null(threshold.col)){
                xline <- seq(0, x.max, length.out=256)
                lines(xline, plot_y(xline), lwd=threshold.lwd[1], lty=threshold.lty[1], col=threshold.col[1])
            }
            visible <- filter_visible_points(log.Quantiles, log.P.values, panel.wh, ht, dpi=dpi)
            visible <- visible & is.finite(log.P.values)
            points(log.Quantiles[visible], log.P.values[visible], col=t(col)[i], pch=19, cex=cex[3])
            if(box) box(lwd=axis.lwd)
        }

        for(mode in requested.mqq){
            for(i in seq_len(R)){
                if(verbose) cat(paste(" Combined ", toupper(mode), " plotting ", trait[i], ".\n", sep=""))
                ht <- ifelse(is.null(height), 6, height)
                wh <- ifelse(is.null(width), 14, width)
                prefix <- ifelse(mode == "mqq", "MQQplot", "QQMplot")
                out.name <- ifelse(is.null(file.name), trait[i], file.name[i])
                if(file.output){
                    open_combined_device(prefix, out.name, wh, ht)
                }else{
                    if(is.null(dev.list())) dev.new(width=wh, height=ht)
                }
                layout(matrix(c(1,2), nrow=1), widths=mqqratio)
                if(mode == "mqq"){
                    draw_combined_manhattan_panel(i, wh * mqqratio[1] / sum(mqqratio), ht)
                    draw_combined_qq_panel(i, wh * mqqratio[2] / sum(mqqratio), ht)
                }else{
                    draw_combined_qq_panel(i, wh * mqqratio[1] / sum(mqqratio), ht)
                    draw_combined_manhattan_panel(i, wh * mqqratio[2] / sum(mqqratio), ht)
                }
                layout(1)
                if(file.output) dev.off()
            }
        }
    }, envir=env)
}

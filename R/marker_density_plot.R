# Marker density plotting utilities used by CMplot().
# The main helper draws chromosome windows colored by marker density.

draw_marker_density_plot <- function(
    chr,
    pos,
    chr.orig.labels,
    col=c("darkgreen", "yellow", "red"),
    main=NULL,
    main.cex=1.2,
    main.font=2,
    chr.labels=NULL, 
    chr.pos.max=FALSE,
    bin=1e6,
    bin.breaks=NULL,
    band=3,
    width=5,
    legend.cex=1,
    legend.y.intersp=1,
    legend.x.intersp=1,
    xticks.pos=1,
    axis.cex=1,
    axis.lwd=1.5,
    plot=TRUE,
    dpi=NULL,
    wh=NULL,
    ht=NULL
)
{
    # Draw or compute the marker density track for each chromosome.
    legend.min <- 1
    legend.max <- NULL
    if(is.null(legend.cex)) legend.cex = 1
    if(!is.null(bin.breaks)){
        bin.breaks <- sort(bin.breaks)
        if(sum(bin.breaks < 0)) stop("breaks should not contain a negative value.")
        if(bin.breaks[1]){
            legend.min <- bin.breaks[1]
        }else{
            bin.breaks <- bin.breaks[-1]
        }
        legend.max <- bin.breaks[length(bin.breaks)]
    }
    if(is.null(col) | length(col) == 1){col=c("darkgreen", "yellow", "red")}
    max.chr <- max(chr)
    chr.num <- unique(chr)
    chorm.maxlen <- max(pos)
    bp <- ifelse(chorm.maxlen < 1e3, 1, ifelse(chorm.maxlen < 1e6, 1e3, 1e6))
    bp_label <- ifelse(bp == 1, "bp", ifelse(bp == 1e3, "Kb", "Mb"))
    if(is.null(main))   main <- paste("The number of SNPs within ", bin / bp, bp_label, " window size", sep="")
    if(plot)    plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chr.num) * band + band), main=main, cex.main=main.cex, font.main=main.font, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
    pos.x <- list()
    chr.pos.max.v <- NULL
    col.index <- list()
    maxbin.num <- NULL
    windinfo <- list()
    for(i in 1 : length(chr.num)){
        pos.x[[i]] <- pos[chr == chr.num[i]]
        maxposindx <- which.max(pos.x[[i]])
        max.pos <- pos.x[[i]][maxposindx]
        chr.pos.max.v <- c(chr.pos.max.v, max.pos)
        cut.breaks <- seq(0, max.pos, bin)
        cut.len <- length(cut.breaks)
        if(cut.breaks[length(cut.breaks)] < max.pos)  cut.breaks <- c(cut.breaks, cut.breaks[length(cut.breaks)] + bin)
        if(chr.pos.max){
            pos.x[[i]] <- pos.x[[i]][-maxposindx]
        }
        if(cut.len <= 1){
            maxbin.num <- c(maxbin.num, length(pos.x[[i]]))
            col.index[[i]] <- rep(length(pos.x[[i]]), length(pos.x[[i]]))
            names(col.index[[i]]) <- 1
        }else{
            # [OPTIMIZED] Integer binning plus tabulate() replaces cut()/table() for marker density counts.
            bin.n <- max(1L, ceiling(max.pos / bin))
            cut.r <- pmin(floor(pos.x[[i]] / bin) + 1L, bin.n)
            eachbin.num <- tabulate(cut.r, nbins=bin.n)
            maxbin.num <- c(maxbin.num, max(eachbin.num))
            col.index[[i]] <- eachbin.num[cut.r]
            names(col.index[[i]]) <- cut.r
        }
        if(plot){
            windinfo <- c(windinfo, tapply(pos.x[[i]], as.numeric(names(col.index[[i]])), function(x){
                return(c(ifelse(!is.null(chr.labels), chr.labels[i], chr.orig.labels[i]),
                    min(x),max(x),length(x)))})
            )
        }
    }
    if(plot){
        windinfo <- as.data.frame(do.call(rbind, windinfo))
        colnames(windinfo) <- c("Chr", "Start", "End", "Num")
        rownames(windinfo) <- NULL
        for(i in 2:ncol(windinfo)){windinfo[, i]<-as.numeric(windinfo[, i])}
    }
    Maxbin.num <- max(maxbin.num)
    maxbin.num <- Maxbin.num
    if(!is.null(legend.max)){
        maxbin.num <- legend.max
    }
    if(Maxbin.num < legend.min)    stop("the maximum number of markers in windows is smaller than the lower boundary of breaks.")
    col=colorRampPalette(col)(maxbin.num - legend.min + 1)
    col.seg=NULL
    for(i in 1 : length(chr.num)){
        if(plot){
            polygon(c(0, 0, chr.pos.max.v[i], chr.pos.max.v[i]), 
            c(-width/5 - band * (i - length(chr.num) - 1), width/5 - band * (i - length(chr.num) - 1), 
            width/5 - band * (i - length(chr.num) - 1), -width/5 - band * (i - length(chr.num) - 1)), col="grey95", border="grey95")
            rect(xleft=0, ybottom = -width/5 - band * (i - length(chr.num) - 1), xright=chr.pos.max.v[i], ytop=width/5 - band * (i - length(chr.num) - 1), border="grey80")
        }
        if(!is.null(legend.max)){
            if(legend.max < Maxbin.num){
                col.index[[i]][col.index[[i]] > legend.max] <- legend.max
            }
        }
        col.index[[i]][col.index[[i]] < legend.min] <- legend.min
        if(!plot)   col.seg <- c(col.seg, col[col.index[[i]] - legend.min + 1])
        if(!is.null(ht) && !is.null(wh) && !is.null(dpi)){
            is_visable <-  filter_visible_points(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), wh * (max(pos.x[[i]])/chorm.maxlen), ht, dpi=dpi)
            if(plot)    segments(pos.x[[i]][is_visable], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]][is_visable], width/5 - band * (i - length(chr.num) - 1), 
                            col=col[col.index[[i]][is_visable] - legend.min + 1], lwd=1)
        }else{
            if(plot)    segments(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]], width/5 - band * (i - length(chr.num) - 1), 
                            col=col[col.index[[i]] - legend.min + 1], lwd=1)
        }
    }
    
    chr.num <- rev(chr.orig.labels)
    if(plot){
        if(!is.null(chr.labels)){
            mtext(at=seq(band, length(chr.num) * band, band), text=chr.labels, side=2, las=2, font=1, cex=axis.cex*0.6, line=0.2, xpd=TRUE)
        }else{
            if(max.chr == 0)    mtext(at=seq(band, length(chr.num) * band, band), text=chr.num, side=2, las=2, font=1, cex=axis.cex*0.6, line=0.2, xpd=TRUE)
            if(max.chr != 0)    mtext(at=seq(band, length(chr.num) * band, band), text=paste("Chr", chr.num, sep=""), side=2, las=2, font=1, cex=axis.cex*0.6, line=0.2, xpd=TRUE)
        }
    }
    if(plot){
        xticks=seq(0, chorm.maxlen / bp, length=10)
        
        if(round(xticks[2]) <= 10){
            xticks=seq(0, chorm.maxlen / bp, round(xticks[2], 1))
        }else{
            xticks=seq(0, chorm.maxlen / bp, round(xticks[2]))    
        }
        
        if((chorm.maxlen/bp - max(xticks)) > 0.5*xticks[2]){
            xticks=c(xticks, round(chorm.maxlen / bp))
        }
        axis(3, mgp=c(3,xticks.pos,0), at=xticks*bp, labels=paste(xticks, bp_label, sep=""), font=1, cex.axis=axis.cex*0.8, tck=0.01, lwd=axis.lwd, padj=1.2)
        axis(3, at=c(0, chorm.maxlen), labels=c("",""), tcl=0, lwd=axis.lwd)
    }

    if(is.null(bin.breaks)){
        legend.len <- 10
        if(maxbin.num <= legend.len)    legend.len <- maxbin.num
        legend.y <- round(seq(0, maxbin.num, length=legend.len + 1))
        legend.y <- unique(legend.y)
        len <- ifelse(length(legend.y)==1, 1, legend.y[2])
        legend.y <- seq(legend.y[2], maxbin.num, len)
    }else{
        legend.y <- bin.breaks
    }
    
    if(!is.null(bin.breaks)){
        if(legend.max < Maxbin.num){
            legend.y[length(legend.y)] <- paste(">=", maxbin.num, sep="")
            legend.y.col <- c(legend.y[c(-length(legend.y))], maxbin.num)
        }else{
            legend.y.col <- legend.y
        }
    }else{
        legend.y.col <- legend.y
    }
    if(legend.min != 1){
        legend.y[1] <- paste("<=", legend.min, sep="")
    }
    legend.y <- c("0", legend.y)
    legend.y.col <- as.numeric(legend.y.col)
    legend.col <- c("grey95", col[legend.y.col - legend.min + 1])
    if(plot){
        legend(x=(chorm.maxlen + chorm.maxlen/50), y=(-width/2.5 + band), title="", legend=legend.y, pch=15, pt.cex=legend.cex*3, col=legend.col,
        cex=legend.cex, bty="n", y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, yjust=0, xjust=0, xpd=TRUE)
        return(windinfo)
    }else{
        return(list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y))
    }
}

draw_marker_density_branch <- function(env) {
    # Run the plot.type == "d" branch inside the CMplot() runtime environment.
    evalq({
        #SNP-Density plot
        wind_snp_num <- NULL
        if("d" %in% plot.type){
            if(verbose) cat(" Marker density plotting.\n")
            if(file.output){
                ht=ifelse(is.null(height), 6, height)
                wh=ifelse(is.null(width), 9, width)
                if(file=="jpg") jpeg(paste("Marker_Density.",ifelse(is.null(file.name),taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                if(file=="pdf") pdf(paste("Marker_Density.",ifelse(is.null(file.name),taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht)
                if(file=="tiff")    tiff(paste("Marker_Density.",ifelse(is.null(file.name),taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                if(file=="png") png(paste("Marker_Density.",ifelse(is.null(file.name),taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                # par(xpd=TRUE)
                par(mar=c(mar[1]-2, mar[2]-1, mar[3]+1, mar[4]))
            }else{
                ht=ifelse(is.null(height), 6, height)
                wh=ifelse(is.null(width), 9, width)
                if(is.null(dev.list())) dev.new(width=wh,height=ht)
                # par(xpd=TRUE)
            }
            wind_snp_num <- draw_marker_density_plot(Pmap[, 1], Pmap[, 2], chr.ori, chr.pos.max=chr.pos.max, dpi=dpi, wh=wh, ht=ht, chr.labels=chr.labels, col=chr.den.col, bin=bin.size, bin.breaks=bin.breaks, main=main[1], main.cex=main.cex, main.font=main.font, legend.cex=legend.cex, xticks.pos=xticks.pos, axis.cex=axis.cex, axis.lwd=axis.lwd)
            if(file.output) dev.off()
        }
    }, envir=env)
}


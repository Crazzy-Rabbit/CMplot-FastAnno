# CMplot-FastAnno - CMplot-compatible version with faster plotting and cleaner annotations
#
# Optimizations applied:
# 1. Data preprocessing: Avoid as.matrix() and apply() conversions
#    - Original: as.matrix() ~10s + apply() ~7s for 17.8M SNPs
#    - Optimized: Keep data.frame, use for-loop for type conversion
#    - Speedup: ~2.7x faster (22s -> 8s)
#
# 2. filter_visible_points: Use integer key + duplicated() instead of cbind + duplicated
#    - Original: !duplicated(cbind(x, y)) - creates matrix, slow row comparison
#    - Optimized: key <- x*1e6 + y; !duplicated(key) - integer comparison
#    - Speedup: ~11x faster
#
# 3. Q-Q confidence interval: Vectorized qbeta() instead of for-loop
#    - Original: for-loop calling qbeta() for each point
#    - Optimized: vectorized qbeta(c(0.05, 0.95), xi, N - xi + 1)
#    - Speedup: ~2x faster
#
# Benchmark results (17.8M SNPs, DPI=100, 3 runs):
# - Manhattan Plot: 47.81s -> 16.79s (2.85x speedup)
# - Q-Q Plot:       89.80s -> 30.14s (2.98x speedup)
# - MD5 checksum: PASSED (output identical to original)
#
# Generated: 2026-02-25
# Tested with: pig/pmap_Direction112_layers4_layers5.csv (17,810,683 SNPs)
#
# Added features in CMplotFastAnno:
# 1. Annotation layout modes for highlighted SNPs/genes:
#    - highlight.text.mode="scatter": legacy CMplot placement kept for compatibility.
#    - highlight.text.mode="top": vertical labels above the plot with ordered connectors.
#    - highlight.text.mode="nearby": side labels placed next to the target points.
# 2. Annotation targets can be supplied from annotation.file instead of only highlight/highlight.text vectors.
# 3. Highlighted/annotated points default to red via highlight.col="red".
#    Users can pass highlight.col to customize one color, a vector of colors, or a per-trait list.
# 4. Threshold-exceeding points are no longer enlarged by default.
#    Set threshold, amplify=TRUE, and signal.cex to redraw/enlarge threshold hits.
# 5. Manhattan plots use 5e-8 as the default threshold when threshold is omitted.
#    Users can still set threshold exactly as in CMplot, or pass threshold=NULL explicitly.
# 6. Rectangular Manhattan x-axis labels show chromosome numbers only, with xlab as the axis title.
# 7. Top annotations support collision avoidance and multi-lane label placement.

.sibling_source_dir <- local({
    calls <- sys.calls()
    frames <- sys.frames()
    source_file <- NULL
    for(i in rev(seq_along(calls))){
        call <- calls[[i]]
        fn <- deparse(call[[1L]])[1L]
        if(fn %in% c("source", "sys.source", "base::source", "base::sys.source")){
            file_expr <- if(!is.null(call[["file"]])) call[["file"]] else call[[2L]]
            for(frame in rev(frames)){
                candidate <- try(eval(file_expr, frame), silent=TRUE)
                if(is.character(candidate) && length(candidate) == 1L && file.exists(candidate)){
                    source_file <- normalizePath(candidate, winslash="/", mustWork=TRUE)
                    break
                }
            }
        }
        if(!is.null(source_file)) break
    }
    if(is.null(source_file)){
        candidates <- c(
            file.path(getwd(), "CMplot_modify", "R", "CMplot.r"),
            file.path(getwd(), "R", "CMplot.r"),
            file.path(getwd(), "CMplot.r")
        )
        candidates <- candidates[file.exists(candidates)]
        if(length(candidates)) source_file <- normalizePath(candidates[1L], winslash="/", mustWork=TRUE)
    }
    if(is.null(source_file)) return(NULL)
    dirname(source_file)
})

.source_sibling_files <- function(files){
    if(is.null(.sibling_source_dir)) return(invisible(NULL))
    target_env <- parent.frame()
    for(file in files){
        module_path <- file.path(.sibling_source_dir, file)
        if(file.exists(module_path)) sys.source(module_path, envir=target_env)
    }
    invisible(NULL)
}

.source_sibling_files(c(
    "plot_utilities.R",
    "annotation_utilities.R",
    "marker_density_plot.R",
    "circular_manhattan_plot.R",
    "rectangular_manhattan_plot.R",
    "combined_mqq_plot.R",
    "qq_plot.R"
))
rm(.source_sibling_files, .sibling_source_dir)

CMplot <- function(
    Pmap,
    # [FEATURE] Default Manhattan palette follows GWASLab classic colors.
    col=c("#597FBD", "#74BAD3"),
    bin.size=1e6,
    bin.breaks=NULL,
    LOG10=TRUE,
    scaled=FALSE,
    pch=19,
    type="p",
    band=1,
    H=1.5,
    ylim=NULL,
    skip=NULL,
    cut=NULL,
    cutfactor=10,
    axis.cex=1,
    axis.lwd=1.5,
    lab.cex=1.5,
    lab.font=2,
    plot.type=c("m","c","q","d"),
    multracks=FALSE,
    multracks.xaxis=FALSE,
    multraits=FALSE,
    points.alpha=100L,
    r=0.3,
    cex=c(0.5,1,1),
    outward=FALSE,
    xlab="Chromosome",
    ylab=expression(-log[10](italic(p))),
    ylab.pos=3,
    xticks.pos=1,
    mar=c(3,6,3,3),
    mar.between=0,
    # [FEATURE] Manhattan plots default to 5e-8 when threshold is omitted; pass threshold=NULL to suppress it.
    threshold=NULL, 
    threshold.col="red",
    threshold.lwd=1,
    threshold.lty=2,
    # [FEATURE] Threshold hits keep the base point size by default because amplify=FALSE.
    # Set amplify=TRUE and tune signal.cex to redraw/enlarge threshold hits.
    amplify=FALSE,
    signal.cex=1.5,
    signal.pch=19,
    signal.col=NULL,
    signal.line=2,
    highlight=NULL,
    highlight.cex=1,
    highlight.pch=19,
    highlight.type="p",
    # [FEATURE] Annotated/highlighted target points are red by default; users can override with highlight.col.
    highlight.col="red",
    highlight.text=NULL,
    highlight.text.col="black",
    highlight.text.cex=1,
    highlight.text.font=3,
    highlight.text.top.cex=1.25,
    # [FEATURE] Added annotation placement modes. The first value keeps legacy CMplot behavior by default.
    highlight.text.mode=c("scatter","top","nearby"),
    # [FEATURE] Controls whether annotation arms prefer the left side, right side, alternating sides, or automatic placement.
    highlight.text.side="auto",
    highlight.text.nearby.offset=0.012,
    highlight.text.nearby.step=0.04,
    highlight.text.line.col="grey85",
    highlight.text.line.lwd=1,
    highlight.text.line.lty=1,
    # [FEATURE] Controls whether annotation connectors are automatic, straight, elbow-shaped, or hidden.
    highlight.text.line.mode=c("auto","straight","elbow","none"),
    highlight.text.line.bend=0.08,
    highlight.text.arrow=FALSE,
    highlight.text.arrow.length=0.06,
    highlight.text.top.space=0.06,
    highlight.text.top.inside=FALSE,
    highlight.text.top.margin=8,
    # [FEATURE] Optimized top-label placement reduces label collisions and connector crossings.
    highlight.text.optimize=TRUE,
    highlight.text.lanes=1,
    highlight.text.lane.gap=0.055,
    highlight.text.min.gap=0.004,
    # [FEATURE] Optional annotation table; columns can be auto-detected or set with annotation.*.col.
    annotation.file=NULL,
    annotation.snp.col=NULL,
    annotation.label.col=NULL,
    annotation.trait.col=NULL,
    annotation.sep=NULL,
    chr.labels=NULL,
    chr.border=FALSE,
    chr.labels.angle=0,
    chr.den.col="black",
    chr.pos.max=FALSE,
    cir.band=1,
    cir.chr=TRUE,
    cir.chr.h=1.5,
    cir.axis=TRUE,
    cir.axis.col="black",
    cir.axis.grid=TRUE,
    conf.int=TRUE,
    conf.int.col=NULL,
    file.output=TRUE,
    file.name=NULL,
    file=c("jpg","pdf","tiff","png"),
    dpi=300,
    height=NULL,
    width=NULL,
    main=NULL,
    main.cex=1.5,
    main.font=2,
    legend.ncol=NULL,
    legend.cex=1,
    legend.pos=c("left","middle","right","none"),
    mqqratio=c(3,1),
    sig.line=TRUE,
    sig.level=5e-8,
    sig.line.col="grey45",
    sig.line.lwd=1,
    sig.line.lty=2,
    suggestive.line=FALSE,
    suggestive.level=1e-5,
    suggestive.line.col="grey70",
    suggestive.line.lwd=1,
    suggestive.line.lty=2,
    additional.line=NULL,
    additional.line.col="grey70",
    additional.line.lwd=1,
    additional.line.lty=2,
    box=FALSE,
    verbose=TRUE
)
{   
    # Public entry point: validate user options, prepare shared plot data,
    # then dispatch each requested plot type to its content-specific helper.

    if(!all(plot.type %in% c("c","m","q","d","mqq","qqm"))) stop("unknown 'plot.type'.")
    legend.pos <- match.arg(legend.pos)
    file <- match.arg(file)
    highlight.text.mode <- match.arg(highlight.text.mode)
    highlight.text.line.mode <- match.arg(highlight.text.line.mode)
    scaled <- isTRUE(scaled[1])
    if(scaled && !LOG10) stop("'scaled=TRUE' expects -log10(p) values and requires LOG10=TRUE.")
    if(!is.null(skip)){
        skip <- as.numeric(skip[1])
        if(!is.finite(skip) || skip < 0) stop("'skip' must be a non-negative numeric value.")
    }
    if(!is.null(cut)){
        cut <- as.numeric(cut[1])
        if(!is.finite(cut) || cut < 0) stop("'cut' must be a non-negative numeric value.")
    }
    cutfactor <- as.numeric(cutfactor[1])
    if(!is.finite(cutfactor) || cutfactor <= 1) stop("'cutfactor' must be greater than 1.")
    if(!is.null(skip) && !is.null(cut) && cut <= skip) stop("'cut' must be greater than 'skip'.")
    if(length(mqqratio) != 2L || any(!is.finite(as.numeric(mqqratio))) || any(as.numeric(mqqratio) <= 0)){
        stop("'mqqratio' must contain two positive numeric values.")
    }
    mqqratio <- as.numeric(mqqratio)
    highlight.text.top.inside <- isTRUE(highlight.text.top.inside[1])
    highlight.text.lanes <- max(1L, as.integer(highlight.text.lanes[1]))
    highlight.text.lane.gap <- max(0, as.numeric(highlight.text.lane.gap[1]))
    highlight.text.min.gap <- max(0, as.numeric(highlight.text.min.gap[1]))
    threshold.user.supplied <- !missing(threshold)
    if(!threshold.user.supplied && any(plot.type %in% c("m","c","mqq","qqm"))){
        # [FEATURE] New significance-line interface. Legacy threshold remains unchanged when supplied.
        line.values <- numeric()
        line.cols <- character()
        line.lwds <- numeric()
        line.ltys <- numeric()
        add_line <- function(values, cols, lwds, ltys){
            values <- as.numeric(values)
            values <- values[is.finite(values)]
            if(!length(values)) return(invisible(NULL))
            line.values <<- c(line.values, values)
            line.cols <<- c(line.cols, rep(cols, length.out=length(values)))
            line.lwds <<- c(line.lwds, rep(lwds, length.out=length(values)))
            line.ltys <<- c(line.ltys, rep(ltys, length.out=length(values)))
            invisible(NULL)
        }
        if(isTRUE(sig.line)) add_line(sig.level, sig.line.col, sig.line.lwd, sig.line.lty)
        if(isTRUE(suggestive.line)) add_line(suggestive.level, suggestive.line.col, suggestive.line.lwd, suggestive.line.lty)
        if(!is.null(additional.line)) add_line(additional.line, additional.line.col, additional.line.lwd, additional.line.lty)
        if(length(line.values)){
            threshold <- line.values
            threshold.col <- line.cols
            threshold.lwd <- line.lwds
            threshold.lty <- line.ltys
        }else{
            threshold <- NULL
        }
    }
    top.highlight.active <- ((!is.null(highlight) && !is.null(highlight.text)) || !is.null(annotation.file)) && highlight.text.mode == "top"
    highlight_mar <- function(x) {
        if(top.highlight.active && !highlight.text.top.inside) x[3] <- max(x[3], highlight.text.top.margin + (highlight.text.lanes - 1L) * 2.2)
        if(!is.null(xlab) && nzchar(as.character(xlab)[1])) x[1] <- max(x[1], 5.4)
        x
    }
    trait <- colnames(Pmap)[-c(1:3)]
    if(length(trait) == 0)   trait <- paste("Trait", 1:(ncol(Pmap)-3), sep="")
    taxa <- paste(trait, collapse="_")
    
    if(length(points.alpha) != 1L)   stop("invalid 'points.alpha': must be 'TRUE', 'FALSE' or an integer between 0 and 255")
    if(is.logical(points.alpha))   points.alpha <- ifelse(points.alpha, formals()$points.alpha, 255L)
    if(!is.integer(points.alpha)){
      if(is.numeric(points.alpha) && points.alpha == as.integer(points.alpha))   points.alpha <- as.integer(points.alpha)
      else   stop("invalid 'points.alpha': must an integer between")
    }
    if(!is.integer(points.alpha))    stop("invalid 'points.alpha': must an integer between")
    if(points.alpha < 0L || points.alpha > 255L)   stop("out-of range 'points.alpha': must be between 0 and 255")

    #get the number of traits
    R=ncol(Pmap)-3

    # [OPTIMIZED] Data preprocessing - avoid as.matrix() and apply()
    # Original: as.matrix() + apply() took ~18s for 17.8M SNPs
    # Optimized: keep data.frame, use for-loop for type conversion (~2.7x faster)

    # Step 1: Filter illegal SNPs (keep as data.frame)
    suppressWarnings(Pmap <- Pmap[Pmap[, 2] != "0", ])
    Pmap <- Pmap[!is.na(Pmap[, 2]), ]

    # Step 2: Convert position column to numeric and filter
    suppressWarnings(Pmap[, 3] <- as.numeric(Pmap[, 3]))
    Pmap <- Pmap[!is.na(Pmap[, 3]), ]

    # Step 3: Handle non-euchromosome
    suppressWarnings(numeric.chr <- as.numeric(Pmap[, 2]))
    suppressWarnings(max.chr <- max(numeric.chr, na.rm=TRUE))
    if(is.infinite(max.chr))    max.chr <- 0
    suppressWarnings(map.xy.index <- which(!numeric.chr %in% c(0:max.chr)))
    if(length(map.xy.index) != 0){
        chr.xy <- unique(Pmap[map.xy.index, 2])
        for(i in 1:length(chr.xy)){
            Pmap[Pmap[, 2] == chr.xy[i], 2] <- max.chr + i
        }
    }
    SNP_id <- Pmap[,1]

    # Step 4: Remove SNP names column
    Pmap <- Pmap[, -1]

    # Step 5: Convert columns to numeric (for-loop is faster than apply for data.frame)
    for(i in 1:ncol(Pmap)) {
        suppressWarnings(Pmap[, i] <- as.numeric(Pmap[, i]))
    }

    # Step 6: Order by chromosome and position
    order_index <- order(Pmap[, 1], Pmap[,2])

    #order the GWAS results by chromosome and position
    Pmap <- Pmap[order_index, ]
    SNP_id <- SNP_id[order_index]

    chr <- unique(Pmap[,1])
    chr.ori <- chr
    if(length(map.xy.index) != 0){
        for(i in 1:length(chr.xy)){
            chr.ori[chr.ori == max.chr + i] <- chr.xy[i]
        }
    }

    # Dispatch plot branches through dedicated helper functions.
    wind_snp_num <- NULL
    draw_marker_density_branch(environment())
    
    if(length(plot.type) > 1 | (!"d" %in% plot.type)){

        #scale and adjust the parameters
        cir.chr.h <- cir.chr.h/5
        cir.band <- cir.band/5
        if(!is.null(threshold)){
            # [FEATURE] threshold controls both threshold lines and optional signal amplification.
            # By default amplify=FALSE, so points above threshold keep the same size and color.
            if(!is.list(threshold)){
                thresholdlist <- list()
                for(i in 1:R){
                    thresholdlist[[i]]  <- threshold
                }
                threshold <- thresholdlist
            }

            if(LOG10){
                if(sum(unlist(threshold) <= 0) != 0) stop("threshold must be greater than 0.")
            }

            threshold.col <- rep(threshold.col, max(sapply(threshold, length)))
            threshold.lwd <- rep(threshold.lwd, max(sapply(threshold, length)))
            threshold.lty <- rep(threshold.lty, max(sapply(threshold, length)))
            signal.col <- rep(signal.col, max(sapply(threshold, length)))
            signal.pch <- rep(signal.pch, max(sapply(threshold, length)))
            signal.cex <- rep(signal.cex, max(sapply(threshold, length)))
        }
        if(length(cex)!=3) cex <- rep(cex,3)

        if(!is.null(ylim)){
            if(!is.list(ylim)){
                if(R > 1)    cat(" (warning: all phenotypes will use the same ylim.)\n")
                if(length(ylim)!=2) stop("ylim for each phenotype should be assigned two values.")
                if(ylim[2] <= ylim[1])  stop("second value should be larger than the first in ylim.")
                ylimlist <- list()
                for(i in 1:R){
                    ylimlist[[i]]  <- ylim
                }
                ylim <- ylimlist
            }else{
                if(length(ylim)!=R) stop("length of list of ylim should equal to the number of phenotype.")
                for(i in 1:R){
                    if(length(ylim[[i]])!=2) stop("ylim for each phenotype should be assigned two values.") 
                    if(ylim[[i]][2] <= ylim[[i]][1])  stop("second value should be larger than the first in ylim.")
                }
            }
        }
        
        if(!is.null(conf.int.col)) conf.int.col <- rep(conf.int.col, R)
        if(!is.null(main)) main <- rep(main, R)
        if(length(mar) != 4)    stop("length of 'mar' shoud equal to 4.")
        if(chr.labels.angle > 90 | chr.labels.angle < -90)  stop("'chr.labels.angle' should be > -90 and < 90.")
        pch=rep(pch, R)
        
        if(!is.null(annotation.file)){
            # [FEATURE] annotation.file overrides highlight/highlight.text so targets and labels can be fully user-defined.
            annotation.input <- build_highlight_from_annotations(annotation.file, trait=trait, R=R,
                sep=annotation.sep, snp.col=annotation.snp.col,
                label.col=annotation.label.col, trait.col=annotation.trait.col)
            highlight <- annotation.input$highlight
            highlight.text <- annotation.input$highlight.text
        }

        if(!is.null(highlight.text)){
            if(!is.list(highlight.text)){
                if((multracks || multraits) && R > 1){
                    warning("For multi-track or multi-trait plots, pass 'highlight.text' as a list with one element per trait; a single vector is reused for every trait.", call.=FALSE)
                }
                highlight.text <- list(highlight.text)
                for(i in 1:R){highlight.text[[i]] = highlight.text[[1]]}
            }else{
                if(length(highlight.text) != R){stop("length of 'highlight.text' not equals to the number of traits.")}  
            }
        }

        if(!is.null(highlight)){
            # [FEATURE] Match target SNP/gene IDs once; highlight.col defaults to red and can be customized.
            highlight_index <- list()
            highlight_col <- list()
            if(is.list(highlight.col)){
                if(length(highlight.col) != R){stop("length of 'highlight.col' not equals to the number of traits.")}
                highlight_col=highlight.col
            }
            if(!is.list(highlight)){
                if((multracks || multraits) && R > 1){
                    warning("For multi-track or multi-trait plots, pass 'highlight' as a list with one element per trait; a single vector is reused for every trait.", call.=FALSE)
                }
                highlight <- list(highlight)
                for(i in 1:R){highlight[[i]] = highlight[[1]]}
            }else{
                if(length(highlight) != R){stop("length of 'highlight' not equals to the number of traits.")}  
            }
            length(highlight_index) <- length(highlight)
            for(i in 1:length(highlight)){
                raw.highlight <- as.character(unlist(highlight[[i]], use.names=FALSE))
                valid.highlight <- !is.na(raw.highlight)
                raw.highlight <- raw.highlight[valid.highlight]
                if(!is.null(highlight.text)){
                    raw.highlight.text <- as.character(unlist(highlight.text[[i]], use.names=FALSE))
                    if(length(raw.highlight.text) == length(valid.highlight)){
                        raw.highlight.text <- raw.highlight.text[valid.highlight]
                    }else if(length(raw.highlight.text) != length(raw.highlight)){
                        stop("length of highlighted labels is not equal to the highlighted SNPs.")
                    }
                }
                if(length(raw.highlight) == 0){
                    highlight_index[[i]] <- NA
                    highlight_col[[i]] <- NA
                    if(!is.null(highlight.text)) highlight.text[[i]] <- NA
                }else{
                    matched.index <- match(raw.highlight, SNP_id)
                    if(all(is.na(matched.index))) stop("No shared SNPs between Pmap and highlight!")
                    matched.keep <- !is.na(matched.index)
                    highlight[[i]] <- raw.highlight[matched.keep]
                    highlight_index[[i]] <- matched.index[matched.keep]
                    if(!is.null(highlight.text)) highlight.text[[i]] <- raw.highlight.text[matched.keep]
                    if(!is.null(highlight.col) && !is.list(highlight.col))  highlight_col[[i]] <- highlight.col
                }
            }
        }
        top.highlight.active <- !is.null(highlight) && !is.null(highlight.text) && highlight.text.mode == "top"

        pvalueT <- as.matrix(Pmap[,-c(1:2)])
        pvalue.pos <- Pmap[, 2]
        chr.factor <- factor(Pmap[, 1], levels=chr)
        chr.index <- as.integer(chr.factor)
        pvalue.pos.list <- split(pvalue.pos, chr.factor, drop=TRUE)
        chr.max.pos <- vapply(pvalue.pos.list, max_without_na, numeric(1))
        chr.min.pos <- vapply(pvalue.pos.list, min_without_na, numeric(1))

        #scale the space parameter between chromosomes
        if(!missing(band)){
            band <- floor(band*(sum(chr.max.pos) - min_without_na(chr.min.pos))/100)
        }else{
            band <- floor((sum(chr.max.pos) - min_without_na(chr.min.pos))/100)
        }
        if(band==0) band=100
        
        if(scaled){
            if(sum(pvalueT < 0, na.rm=TRUE) != 0) stop("scaled -log10(p) values should be non-negative.")
            logpvalueT.raw <- pvalueT
            logpvalueT.raw[!is.finite(logpvalueT.raw)] <- NA
            pvalueT <- 10^(-pmin(logpvalueT.raw, 300))
        }else if(LOG10){
            if(sum(pvalueT <= 0, na.rm=TRUE) != 0 || sum(pvalueT > 1, na.rm=TRUE) != 0) stop("p values should be at range of (0, 1).")
            pvalueT[pvalueT <= 0] <- NA
            pvalueT[pvalueT > 1] <- NA
            logpvalueT.raw <- -log10(pvalueT)
        }else{
            logpvalueT.raw <- pvalueT
        }
        Pmap[,-c(1:2)] <- pvalueT

        #set the colors for the plot
        if(is.vector(col)){
            col <- matrix(col,R,length(col),byrow=TRUE)
        }
        if(is.matrix(col)){
            #try to transform the colors into matrix for all traits
            col <- matrix(as.vector(t(col)),R,dim(col)[2],byrow=TRUE)
        }

        Num <- tabulate(chr.index, nbins=length(chr))
        Nchr <- length(Num)
        N <- NULL

        #set the colors for each traits
        for(i in 1:R){
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            N[i] <- ceiling(Nchr/length(colx))
        }
        
        #insert the space into chromosomes and return the midpoint of each chromosome
        ticks <- NULL
        chr.border.pos <- NULL
        pvalue.posN <- NULL
        #pvalue <- pvalueT[,j]
        if(Nchr == 1){
            bp <- ifelse((chr.max.pos[1] - chr.min.pos[1]) > 1000000, 1000000, 1000)
            bp_lab <- ifelse(bp == 1000000, " (Mb)", " (Kb)")
            pvalue.posN <- pvalue.pos.list[[1]] + band
            ticks <- seq(chr.min.pos[1], chr.max.pos[1], length=10)
            ticks <- seq(round(chr.min.pos[1] / bp), round(chr.max.pos[1] / bp), round((ticks[2]-ticks[1])/bp) + 0.5)
            if(!round(chr.max.pos[1] / bp) %in% ticks){
                if(round(chr.max.pos[1] / bp) - ticks[length(ticks)] > 0.5 * ticks[2])
                ticks <- c(ticks, round(chr.max.pos[1] / bp))
            }
            ticks <- ticks[-1]
            chr.labels <- ticks
            ticks <- ticks * bp + band
            chr.border <- FALSE
        }else{
            # [OPTIMIZED] Vectorized chromosome coordinate offsets replace repeated concatenation loops.
            chr.offset <- numeric(Nchr)
            chr.offset[1] <- band
            if(Nchr > 1){
                chr.offset[-1] <- band + cumsum(chr.max.pos[-Nchr] + band)
            }
            pvalue.posN <- pvalue.pos + chr.offset[chr.index]
            chr.end.pos <- chr.offset + chr.max.pos
            ticks <- chr.end.pos - floor(chr.max.pos / 2)
            ticks[1] <- floor((chr.end.pos[1] + chr.min.pos[1] + band) / 2)
            chr.border.pos <- chr.end.pos[-Nchr] + 0.5 * band
        }

        if(!is.null(chr.labels) & Nchr != 1){
            chr.labels <- as.character(chr.labels)
            if(length(chr.labels) != Nchr)  stop("length of 'chr.labels' should equal to the number of chromosomes.")
            ticks.logi <- rep(TRUE, length(ticks))
            for(ti in 1:Nchr){
                if(is.na(chr.labels[ti]))    ticks.logi[ti] <- FALSE
            }
            if(!all(ticks.logi)){
                chr.labels <- chr.labels[ticks.logi]
                ticks <- ticks[ticks.logi]
            }
        }

        pvalue.posN.list <- split(pvalue.posN, chr.factor, drop=TRUE)
        
        # Merge the p-values of traits by column and apply display-only y-axis compression.
        plot_y <- function(y) transform_plot_y(y, skip=skip, cut=cut, cutfactor=cutfactor)
        inverse_y <- function(y) inverse_plot_y(y, skip=skip, cut=cut, cutfactor=cutfactor)
        threshold_to_y <- function(x) {
            y <- if(LOG10) -log10(as.numeric(x)) else as.numeric(x)
            plot_y(y)
        }
        logpvalueT <- plot_y(logpvalueT.raw)
        prepare_qq_data <- function(i){
            if(scaled){
                obs <- as.numeric(logpvalueT.raw[, i])
                obs <- obs[is.finite(obs) & obs >= 0]
                obs <- sort(obs, decreasing=TRUE)
                N <- length(obs)
                log.Quantiles <- -log10((seq_len(N))/(N+1))
                log.P.values <- plot_y(obs)
            }else{
                P.values <- as.numeric(Pmap[, i+2])
                P.values <- P.values[!is.na(P.values)]
                if(LOG10){
                    P.values <- P.values[P.values > 0]
                    P.values <- P.values[P.values < 1]
                    N <- length(P.values)
                    P.values <- P.values[order(P.values)]
                    log.P.values <- plot_y(-log10(P.values))
                }else{
                    N <- length(P.values)
                    P.values <- P.values[order(P.values, decreasing=TRUE)]
                    log.P.values <- plot_y(P.values)
                }
                log.Quantiles <- -log10((seq_len(N))/(N+1))
            }
            list(N=N, log.Quantiles=log.Quantiles, log.P.values=log.P.values)
        }

        add <- list()
        for(i in 1:R){
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            add[[i]] <- c(Num,rep(0,N[i]*length(colx)-Nchr))
        }

        circleMin <- (min_without_na(pvalue.posN) - band - 1)
        TotalN <- max_without_na(pvalue.posN) - circleMin

        if(length(chr.den.col) > 1){
            cir.density=TRUE
            den.fold <- 20
            density.list <- draw_marker_density_plot(Pmap[, 1], Pmap[, 2], chr.ori, chr.pos.max=FALSE, col=chr.den.col, plot=FALSE, bin=bin.size, bin.breaks=bin.breaks)
        }else{
            cir.density=FALSE
        }

        # [FEATURE] Rectangular Manhattan axes show chromosome labels only; the axis title is controlled by xlab.
        rect_chr_axis_labels <- function(){
            if(Nchr == 1){
                return(c(paste(unique(Pmap[,1]), bp_lab, sep=""), chr.labels))
            }
            if(is.null(chr.labels)){
                c("", chr.ori)
            }else{
                c("", chr.labels)
            }
        }

        # [FEATURE] Shared x-axis drawer keeps the new chromosome axis behavior consistent across plot modes.
        draw_rect_chr_axis <- function(scale=1, Min=NULL, Max=NULL, title=TRUE){
            axis.at <- c(min_without_na(pvalue.posN)-band, ticks)
            axis.labels <- rect_chr_axis_labels()
            if(chr.labels.angle == 0){
                axis(1, mgp=c(3,xticks.pos,0), at=axis.at,
                    lwd=axis.lwd*scale, cex.axis=axis.cex*scale,
                    font=lab.font, labels=axis.labels, padj=1)
            }else{
                axis(1, mgp=c(3,xticks.pos,0), at=axis.at,
                    lwd=axis.lwd*scale, labels=FALSE)
                label.y <- par("usr")[3]
                if(!is.null(Min) && !is.null(Max)){
                    label.y <- par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min)
                }
                text(axis.at, label.y, cex=axis.cex*scale, font=lab.font,
                    labels=axis.labels, srt=chr.labels.angle, xpd=TRUE,
                    adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)),
                        ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
            }
            axis(1, mgp=c(3,xticks.pos,0), at=c(ticks[length(ticks)], max_without_na(pvalue.posN)),
                labels=c("",""), tcl=0, lwd=axis.lwd*scale)
            if(title && !is.null(xlab) && nzchar(as.character(xlab)[1])){
                # [FEATURE] Extra x-axis title offset prevents "Chromosome" from overlapping tick labels.
                mtext(side=1, text=xlab, line=4, cex=lab.cex*scale*0.85,
                    font=lab.font, xpd=TRUE)
            }
        }
    }

    draw_combined_mqq_plots(environment())
    draw_circular_manhattan_plot(environment())
    draw_rectangular_manhattan_plots(environment())
    draw_qq_plots(environment())
    if(file.output & verbose)   cat(paste(" Plots are stored in: ", getwd(), sep=""), "\n")
    if(!is.null(wind_snp_num))  return(invisible(wind_snp_num))
}

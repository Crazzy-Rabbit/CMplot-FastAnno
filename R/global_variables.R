# Global variable declarations for R CMD check.
# Plot branch helpers are evaluated inside the CMplot() runtime environment, so
# these names are deliberately supplied by CMplot() rather than by package scope.
if(getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "H", "LOG10", "N", "Nchr", "Pmap", "R", "TotalN", "add",
        "amplify", "axis.cex", "axis.lwd", "band", "bin.breaks",
        "bin.size", "bp_lab", "cex", "chr", "chr.border",
        "chr.border.pos", "chr.den.col", "chr.labels", "chr.ori",
        "chr.pos.max", "cir.axis", "cir.axis.col", "cir.axis.grid",
        "cir.band", "cir.chr", "cir.chr.h", "cir.density",
        "circleMin", "conf.int", "conf.int.col", "den.fold",
        "density.list", "dpi", "draw_rect_chr_axis", "file",
        "file.name", "file.output", "height", "highlight",
        "highlight.cex", "highlight.col", "highlight.pch",
        "highlight.text", "highlight.text.arrow",
        "highlight.text.arrow.length", "highlight.text.cex",
        "highlight.text.col", "highlight.text.font",    
        "highlight.text.lane.gap", "highlight.text.lanes",
        "highlight.text.line.bend", "highlight.text.line.col",
        "highlight.text.line.lty", "highlight.text.line.lwd",
        "highlight.text.line.mode", "highlight.text.min.gap",
        "highlight.text.mode", "highlight.text.nearby.offset",
        "highlight.text.nearby.step", "highlight.text.optimize",
        "highlight.text.side", "highlight.text.top.cex",
        "highlight.text.top.inside", "highlight.text.top.space",
        "highlight.type", "highlight_col",
        "highlight_index", "highlight_mar", "lab.cex", "lab.font",
        "legend.cex", "legend.ncol", "legend.pos", "logpvalueT",
        "inverse_y", "main", "main.cex", "main.font", "mar", "mar.between",
        "multracks", "multracks.xaxis", "multraits", "outward", "pch",
        "mqqratio", "plot.type", "plot_y", "points.alpha", "prepare_qq_data",
        "pvalue.posN", "pvalue.posN.list", "pvalueT", "r", "signal.cex", "signal.col", "signal.line",
        "signal.pch", "taxa", "threshold", "threshold.col",
        "threshold.lty", "threshold.lwd", "threshold_to_y", "ticks", "trait", "type",
        "verbose", "wh", "width", "wind_snp_num", "xticks.pos", "ylab",
        "ylab.pos", "ylim", "skip"
    ))
}

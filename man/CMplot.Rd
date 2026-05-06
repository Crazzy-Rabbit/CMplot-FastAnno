\name{CMplot}
\alias{CMplot}
\title{Draw Manhattan, Circular Manhattan, Q-Q, and Marker Density Plots}
\description{
Draws CMplot-compatible genomic visualizations from a marker map and one or
more trait columns. Supported plot types include rectangular Manhattan plots,
circular Manhattan plots, Q-Q plots, and marker density plots.
}
\usage{
CMplot(Pmap, col, bin.size, bin.breaks, LOG10, scaled, pch, type, band, H,
    ylim, skip, cut, cutfactor, axis.cex, axis.lwd, lab.cex, lab.font, plot.type, multracks,
    multracks.xaxis, multraits, points.alpha, r, cex, outward, xlab, ylab,
    ylab.pos, xticks.pos, mar, mar.between, threshold, threshold.col,
    threshold.lwd, threshold.lty, amplify, signal.cex, signal.pch,
    signal.col, signal.line, highlight, highlight.cex, highlight.pch,
    highlight.type, highlight.col, highlight.text, highlight.text.col,
    highlight.text.cex, highlight.text.font, highlight.text.top.cex,
    highlight.text.mode, highlight.text.side, highlight.text.nearby.offset,
    highlight.text.nearby.step, highlight.text.line.col,
    highlight.text.line.lwd, highlight.text.line.lty,
    highlight.text.line.mode, highlight.text.line.bend,
    highlight.text.arrow, highlight.text.arrow.length,
    highlight.text.top.space, highlight.text.top.inside,
    highlight.text.top.margin,
    highlight.text.optimize, highlight.text.lanes, highlight.text.lane.gap,
    highlight.text.min.gap, annotation.file, annotation.snp.col,
    annotation.label.col, annotation.trait.col, annotation.sep, chr.labels,
    chr.border, chr.labels.angle, chr.den.col, chr.pos.max, cir.band,
    cir.chr, cir.chr.h, cir.axis, cir.axis.col, cir.axis.grid, conf.int,
    conf.int.col, file.output, file.name, file, dpi, height, width, main,
    main.cex, main.font, legend.ncol, legend.cex, legend.pos, mqqratio,
    sig.line, sig.level, sig.line.col, sig.line.lwd, sig.line.lty,
    suggestive.line, suggestive.level, suggestive.line.col,
    suggestive.line.lwd, suggestive.line.lty, additional.line,
    additional.line.col, additional.line.lwd, additional.line.lty,
    box, verbose)
}
\arguments{
\item{Pmap}{A data frame whose first three columns are marker ID, chromosome,
and position, followed by one or more numeric trait or p-value columns.}
\item{col}{Colors used for Manhattan and Q-Q points.}
\item{bin.size}{Window size used for marker density plots.}
\item{bin.breaks}{Optional marker density legend break points.}
\item{LOG10}{Logical; transform p-values with \code{-log10}.}
\item{scaled}{Logical; when \code{TRUE}, trait columns are interpreted as
already-scaled \code{-log10(p)} values. Threshold arguments remain on the
p-value scale when \code{LOG10=TRUE}.}
\item{pch}{Point plotting character.}
\item{type}{Base graphics point or line type.}
\item{band}{Spacing between chromosomes.}
\item{H}{Track height for circular Manhattan plots.}
\item{ylim}{Optional y-axis limits, either a vector or per-trait list.}
\item{skip}{Optional lower \code{-log10(p)} display cutoff. Points below this
value are omitted from Manhattan/Q-Q display coordinates.}
\item{cut}{Optional upper \code{-log10(p)} display cutoff. Values above this
point are compressed for display.}
\item{cutfactor}{Compression factor used for values above \code{cut}.}
\item{axis.cex}{Axis label size multiplier.}
\item{axis.lwd}{Axis line width.}
\item{lab.cex}{Axis title and label size multiplier.}
\item{lab.font}{Font used for labels.}
\item{plot.type}{Character vector selecting plot types: \code{"m"},
\code{"c"}, \code{"q"}, \code{"d"}, \code{"mqq"}, and/or \code{"qqm"}.}
\item{multracks}{Logical; draw multi-trait plots as separate tracks.}
\item{multracks.xaxis}{Logical; show x-axis labels on each multi-track panel.}
\item{multraits}{Logical; draw multiple traits in one plot.}
\item{points.alpha}{Integer alpha channel from 0 to 255 for point colors.}
\item{r}{Inner radius for circular Manhattan plots.}
\item{cex}{Point size settings.}
\item{outward}{Logical; draw circular tracks outward when \code{TRUE}.}
\item{xlab}{X-axis title.}
\item{ylab}{Y-axis title.}
\item{ylab.pos}{Axis title line position.}
\item{xticks.pos}{Tick label position.}
\item{mar}{Base graphics margins.}
\item{mar.between}{Margin between multi-track panels.}
\item{threshold}{Optional significance threshold or thresholds.}
\item{threshold.col}{Threshold line color.}
\item{threshold.lwd}{Threshold line width.}
\item{threshold.lty}{Threshold line type.}
\item{amplify}{Logical; redraw points beyond thresholds with signal styling.}
\item{signal.cex}{Point size for amplified signal points.}
\item{signal.pch}{Point character for amplified signal points.}
\item{signal.col}{Color for amplified signal points.}
\item{signal.line}{Line width for circular signal guide lines.}
\item{highlight}{Marker IDs to highlight. For multi-track or multi-trait
plots, use a list with one element per trait when each trait has its own
target markers.}
\item{highlight.cex}{Point size multiplier for highlighted markers.}
\item{highlight.pch}{Point character for highlighted markers.}
\item{highlight.type}{Point or line type for highlighted markers.}
\item{highlight.col}{Color for highlighted markers.}
\item{highlight.text}{Labels for highlighted markers. For multi-track or
multi-trait plots, use a list with one element per trait to match
\code{highlight}.}
\item{highlight.text.col}{Color for highlighted labels.}
\item{highlight.text.cex}{Size multiplier for highlighted labels.}
\item{highlight.text.font}{Font for highlighted labels.}
\item{highlight.text.top.cex}{Additional size multiplier for top labels.}
\item{highlight.text.mode}{Label placement mode: \code{"scatter"},
\code{"top"}, or \code{"nearby"}.}
\item{highlight.text.side}{Preferred side for nearby labels.}
\item{highlight.text.nearby.offset}{Horizontal offset for nearby labels.}
\item{highlight.text.nearby.step}{Vertical step used to avoid nearby-label overlap.}
\item{highlight.text.line.col}{Connector line color for highlighted labels.}
\item{highlight.text.line.lwd}{Connector line width.}
\item{highlight.text.line.lty}{Connector line type.}
\item{highlight.text.line.mode}{Connector style: \code{"auto"},
\code{"straight"}, \code{"elbow"}, or \code{"none"}.}
\item{highlight.text.line.bend}{Bend amount for elbow connectors.}
\item{highlight.text.arrow}{Logical; draw arrowheads on connectors.}
\item{highlight.text.arrow.length}{Arrowhead length.}
\item{highlight.text.top.space}{Space between the plot area and top labels.}
\item{highlight.text.top.inside}{Logical; place top labels inside the current
plot panel instead of reserving external top margin. This is useful for compact
multi-track Manhattan plots.}
\item{highlight.text.top.margin}{Top margin reserved for top labels.}
\item{highlight.text.optimize}{Logical; optimize top-label placement.}
\item{highlight.text.lanes}{Number of lanes available for top labels.}
\item{highlight.text.lane.gap}{Vertical gap between top-label lanes.}
\item{highlight.text.min.gap}{Minimum horizontal gap between top labels.}
\item{annotation.file}{Optional annotation file path or data frame.}
\item{annotation.snp.col}{Column name or index containing marker IDs.}
\item{annotation.label.col}{Column name or index containing annotation labels.}
\item{annotation.trait.col}{Column name or index containing trait IDs.}
\item{annotation.sep}{Field separator for annotation files.}
\item{chr.labels}{Optional chromosome labels.}
\item{chr.border}{Logical; draw chromosome border lines.}
\item{chr.labels.angle}{Angle for chromosome labels.}
\item{chr.den.col}{Color vector for chromosome density tracks.}
\item{chr.pos.max}{Logical; use chromosome maximum position in density plots.}
\item{cir.band}{Spacing between circular tracks.}
\item{cir.chr}{Logical; draw chromosome ring in circular Manhattan plots.}
\item{cir.chr.h}{Circular chromosome ring height.}
\item{cir.axis}{Logical; draw circular y-axis.}
\item{cir.axis.col}{Circular axis color.}
\item{cir.axis.grid}{Logical; draw circular axis grid lines.}
\item{conf.int}{Logical; draw Q-Q confidence interval.}
\item{conf.int.col}{Q-Q confidence interval color.}
\item{file.output}{Logical; write plots to files.}
\item{file.name}{Optional output file name or names.}
\item{file}{Output file type.}
\item{dpi}{Output resolution for raster files.}
\item{height}{Output height.}
\item{width}{Output width.}
\item{main}{Plot title or titles.}
\item{main.cex}{Plot title size multiplier.}
\item{main.font}{Plot title font.}
\item{legend.ncol}{Number of columns in legends.}
\item{legend.cex}{Legend size multiplier.}
\item{legend.pos}{Legend position.}
\item{mqqratio}{Two positive numbers controlling the relative panel widths in
combined Manhattan plus Q-Q plots.}
\item{sig.line}{Logical; draw the genome-wide significance line when
\code{threshold} is not explicitly supplied.}
\item{sig.level}{P-value level for the genome-wide significance line.}
\item{sig.line.col}{Color for the genome-wide significance line.}
\item{sig.line.lwd}{Line width for the genome-wide significance line.}
\item{sig.line.lty}{Line type for the genome-wide significance line.}
\item{suggestive.line}{Logical; draw a suggestive significance line when
\code{threshold} is not explicitly supplied.}
\item{suggestive.level}{P-value level for the suggestive line.}
\item{suggestive.line.col}{Color for the suggestive line.}
\item{suggestive.line.lwd}{Line width for the suggestive line.}
\item{suggestive.line.lty}{Line type for the suggestive line.}
\item{additional.line}{Optional extra p-value levels to draw when
\code{threshold} is not explicitly supplied.}
\item{additional.line.col}{Color or colors for additional lines.}
\item{additional.line.lwd}{Line width or widths for additional lines.}
\item{additional.line.lty}{Line type or types for additional lines.}
\item{box}{Logical; draw plot box.}
\item{verbose}{Logical; print progress messages.}
}
\value{
For marker density plots, invisibly returns the window-level marker density
table. For other plot types, returns \code{NULL} invisibly.
}
\details{
This package is implemented with R's standard graphics stack and does not
require external CRAN or Bioconductor packages.
}
\keyword{hplot}

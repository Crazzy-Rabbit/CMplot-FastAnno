# Rectangular Manhattan plot branches used by CMplot().
# This function draws single-trait, multi-track, and multi-trait Manhattan plots.

draw_rectangular_manhattan_plots <- function(env) {
    # Run the plot.type == "m" branch inside the CMplot() runtime environment.
    evalq({
        if("m" %in% plot.type){
        
            is_visable <- list()
            for(i in 1:R){
                if(file.output){
                    ht=ifelse(is.null(height), 6, height)
                    wh=ifelse(is.null(width), 14, width)
                    is_visable[[i]] <- filter_visible_points(pvalue.posN, logpvalueT[,i], wh, ht, dpi=dpi)
                }else{
                    is_visable[[i]] <- rep(TRUE, nrow(logpvalueT))
                }
            }
        
            if(multracks | multraits){
                if(R < 2)   stop("need more than one trait.")
                if(multracks){
                    if(file.output){
                        ht=ifelse(is.null(height), 6, height)
                        wh=ifelse(is.null(width), 14, width)
                        if(file=="jpg") jpeg(paste("Multi-tracks_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi*R,res=dpi,quality=100)
                        if(file=="pdf") pdf(paste("Multi-tracks_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht*R)
                        if(file=="tiff")    tiff(paste("Multi-tracks_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi*R,res=dpi)
                        if(file=="png") png(paste("Multi-tracks_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi*R,res=dpi,bg=NA)
                        par(mfcol=c(R,1), xaxs="i")
                    }
                    if(!file.output){
                        ht=ifelse(is.null(height), 6, height)
                        wh=ifelse(is.null(width), 14, width)
                        if(is.null(dev.list())) dev.new(width=wh, height=ht)
                        # par(xpd=TRUE)
                    }
                    for(i in 1:R){
                        # Add room for x axis, if there are multiple
                        btwn_adj=if(multracks.xaxis) 2 else 0
                        if(i == 1)  plot.mar <- c(mar.between + btwn_adj, mar[2]+1, mar[3], 0)
                        if(i == R)  plot.mar <- c(mar[1]+1, mar[2]+1, 0, 0)
                        if(i != 1 & i != R) plot.mar <- c(mar.between + btwn_adj, mar[2]+1, 0, 0)
                        par(mar=highlight_mar(plot.mar))
                        if(verbose) cat(paste(" Multi-tracks Manhattan plotting ",trait[i],".\n",sep=""))
                        colx=col[i,]
                        colx=colx[!is.na(colx)]
                        pvalue=pvalueT[,i]
                        logpvalue=logpvalueT[,i]
                        if(is.null(ylim)){
                            if(!is.null(threshold[[i]])){
                                # if(sum(threshold!=0)==length(threshold)){
                                    if(LOG10){
                                        Max=round_y_axis_max(max_without_na(c((-log10(min_without_na(pvalue))),-log10(min_without_na(threshold[[i]])))))
                                        Min <- round_y_axis_min(min_without_na(c((-log10(max_without_na(pvalue))),-log10(max_without_na(threshold[[i]])))))
                                    }else{
                                        Max=round_y_axis_max(max_without_na(c((max_without_na(pvalue)),max_without_na(threshold[[i]]))))
                                        #if(abs(Max)<=1)    Max=max_without_na(c(max_without_na(pvalue),max_without_na(threshold)))
                                        Min<-round_y_axis_min(min_without_na(c((min_without_na(pvalue)),min_without_na(threshold[[i]]))))
                                        #if(abs(Min)<=1)    Min=min_without_na(min_without_na(pvalue),min_without_na(threshold))
                                    }
                            }else{
                                if(LOG10){
                                        Max=round_y_axis_max((-log10(min_without_na(pvalue))))
                                        Min<-round_y_axis_min((-log10(max_without_na(pvalue))))
                                }else{
                                        Max=round_y_axis_max((max_without_na(pvalue)))
                                        #if(abs(Max)<=1)    Max=max_without_na(max_without_na(pvalue))
                                        Min=round_y_axis_min((min_without_na(pvalue)))
                                        #if(abs(Min)<=1)    Min=min_without_na(min_without_na(pvalue))
                                        # }else{
                                            # Max=max_without_na(ceiling(max_without_na(pvalue)))
                                        # }
                                }
                            }
                            if((Max-Min)<=1){
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]*(R/2)+1,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,max_without_na(pvalue.posN)+band),ylim=c(Min,Max),ann=FALSE,
                                    cex.axis=axis.cex*(R/2),font=lab.font,axes=FALSE,yaxs="i")
                            }else{
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]*(R/2)+1,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,max_without_na(pvalue.posN)+band),ylim=c(Min,Max),ann=FALSE,
                                    cex.axis=axis.cex*(R/2),font=lab.font,axes=FALSE,yaxs="i")
                            }
                            mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex*(R/2), font=lab.font, xpd=TRUE)
                        }else{
                            Max <- max_without_na(ylim[[i]])
                            Min <- min_without_na(ylim[[i]])
                            plot(pvalue.posN[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]*(R/2)+1,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,max_without_na(pvalue.posN)+band),ylim=ylim[[i]],ann=FALSE,
                                cex.axis=axis.cex*(R/2),font=lab.font,axes=FALSE,yaxs="i")
                            mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex*(R/2), font=lab.font, xpd=TRUE)
                        }
        
                        if(chr.border){
                            for(b in 1:length(chr.border.pos)){
                                segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=axis.lwd, lty=2)
                            }
                        }
        
                        #add the names of traits on plot 
                        if(legend.pos=="left"){
                            text(min_without_na(pvalue.posN),Max,labels=trait[i],adj=c(-0.2, 1.2),font=4,cex=legend.cex*(R/2),xpd=TRUE) 
                        }else if(legend.pos=="middle"){
                            text((max_without_na(pvalue.posN)+min_without_na(pvalue.posN))/2,Max,labels=trait[i],adj=c(0.5, 1.2),font=4,cex=legend.cex*(R/2),xpd=TRUE) 
                        }else if(legend.pos=="right"){
                            text(max_without_na(pvalue.posN),Max,labels=trait[i],adj=c(1.2, 1.2),font=4,cex=legend.cex*(R/2),xpd=TRUE) 
                        }
                    
                        if(i == R || multracks.xaxis){
                            draw_rect_chr_axis(scale=R/2, Min=Min, Max=Max, title=(i == R))
                        }
                        #if(i==1) mtext("Manhattan plot",side=3,padj=-1,font=lab.font,cex=xn)
                        if(is.null(ylim)){
                            if((Max-Min)>1){
                                axis(2, las=1,lwd=axis.lwd*(R/2),cex.axis=axis.cex*(R/2),font=lab.font)
                                axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=axis.lwd*(R/2))
                            }else{
                                axis(2,las=1,lwd=axis.lwd*(R/2),cex.axis=axis.cex*(R/2),font=lab.font)
                                axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=axis.lwd*(R/2))
                            }
                        }else{
                            axis(2, las=1,lwd=axis.lwd*(R/2),cex.axis=axis.cex*(R/2),font=lab.font)
                            axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=axis.lwd*(R/2))
                        }
                        if(!is.null(threshold[[i]])){
                            for(thr in 1:length(threshold[[i]])){
                                h <- threshold_to_y(threshold[[i]][thr])
                                segments(0, h, max_without_na(pvalue.posN), h, col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                            }
                            if(amplify==TRUE){
                                if(LOG10){
                                    threshold[[i]] <- sort(threshold[[i]])
                                    sgline1=threshold_to_y(max_without_na(threshold[[i]]))
                                }else{
                                    threshold[[i]] <- sort(threshold[[i]], decreasing=TRUE)
                                    sgline1=threshold_to_y(min_without_na(threshold[[i]]))
                                }
                                sgindex=which(logpvalue>=sgline1)
                                HY1=logpvalue[sgindex]
                                HX1=pvalue.posN[sgindex]
                                
                                #cover the points that exceed the threshold with the color "white"
                                points(HX1,HY1,pch=pch,cex=cex[2]*R,col="white")
                                
                                for(ll in 1:length(threshold[[i]])){
                                    if(ll == 1){
                                        if(LOG10){
                                            sgline1=threshold_to_y(threshold[[i]][ll])
                                        }else{
                                            sgline1=threshold_to_y(threshold[[i]][ll])
                                        }
                                        sgindex=which(logpvalue>=sgline1)
                                        HY1=logpvalue[sgindex]
                                        HX1=pvalue.posN[sgindex]
                                    }else{
                                        if(LOG10){
                                            sgline0=threshold_to_y(threshold[[i]][ll-1])
                                            sgline1=threshold_to_y(threshold[[i]][ll])
                                        }else{
                                            sgline0=threshold_to_y(threshold[[i]][ll-1])
                                            sgline1=threshold_to_y(threshold[[i]][ll])
                                        }
                                        sgindex=which(logpvalue>=sgline1 & logpvalue < sgline0)
                                        HY1=logpvalue[sgindex]
                                        HX1=pvalue.posN[sgindex]
                                    }
        
                                    if(is.null(signal.col)){
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*R,col=rep(rep(colx,N[i]),add[[i]])[sgindex])
                                    }else{
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*R,col=signal.col[ll])
                                    }
                                }
                            }
                        }
        
                        if(!is.null(highlight)){
                            # points(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],pch=pch,cex=cex[2]*R,col="white")
                            if(!is.na(highlight_index[[i]][1])){
                                # [FEATURE] Multi-track highlighted points keep the base point size by default.
                                # highlight.cex is a multiplier, so users can still scale annotated points.
                                highlight.point.cex <- highlight.cex * cex[2] * (R/2)
                                if(is.null(highlight.col)){
                                    draw_highlight_labels(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_without_na(pvalue.posN)-band,max_without_na(pvalue.posN)+band),ylim=c(Min,Max),words=highlight.text[[i]],point.cex=highlight.point.cex,text.cex=highlight.text.cex*R/2, pch=highlight.pch,type=highlight.type,point.col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]],text.col=highlight.text.col,text.font=highlight.text.font,top.cex=highlight.text.top.cex,mode=highlight.text.mode,side=highlight.text.side,nearby.offset=highlight.text.nearby.offset,nearby.step=highlight.text.nearby.step,line.col=highlight.text.line.col,line.lwd=highlight.text.line.lwd,line.lty=highlight.text.line.lty,line.mode=highlight.text.line.mode,line.bend=highlight.text.line.bend,arrow=highlight.text.arrow,arrow.length=highlight.text.arrow.length,top.space=highlight.text.top.space,top.inside=highlight.text.top.inside,optimize=highlight.text.optimize,lanes=highlight.text.lanes,lane.gap=highlight.text.lane.gap,min.gap=highlight.text.min.gap)
                                }else{
                                    draw_highlight_labels(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_without_na(pvalue.posN)-band,max_without_na(pvalue.posN)+band),ylim=c(Min,Max),words=highlight.text[[i]],point.cex=highlight.point.cex,text.cex=highlight.text.cex*R/2, pch=highlight.pch,type=highlight.type,point.col=highlight_col[[i]],text.col=highlight.text.col,text.font=highlight.text.font,top.cex=highlight.text.top.cex,mode=highlight.text.mode,side=highlight.text.side,nearby.offset=highlight.text.nearby.offset,nearby.step=highlight.text.nearby.step,line.col=highlight.text.line.col,line.lwd=highlight.text.line.lwd,line.lty=highlight.text.line.lty,line.mode=highlight.text.line.mode,line.bend=highlight.text.line.bend,arrow=highlight.text.arrow,arrow.length=highlight.text.arrow.length,top.space=highlight.text.top.space,top.inside=highlight.text.top.inside,optimize=highlight.text.optimize,lanes=highlight.text.lanes,lane.gap=highlight.text.lane.gap,min.gap=highlight.text.min.gap)
                                }
                            }
                        }
                        if(!is.null(main) & R == 1)  title(main=main[1], cex.main=main.cex, font.main= main.font)
                        if(box) box(lwd=axis.lwd)
                        #if(!is.null(threshold) & !is.null(signal.line))    abline(v=pvalue.posN[which(pvalueT[,i] < min_without_na(threshold))],col="grey",lty=2,lwd=signal.line)
                    }
                    if(file.output) dev.off()
                }
                if(multraits){
                    if(file.output){
                        ht=ifelse(is.null(height), 6, height)
                        wh=ifelse(is.null(width), 14, width)
                        if(file=="jpg") jpeg(paste("Multi-traits_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                        if(file=="pdf") pdf(paste("Multi-traits_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht)
                        if(file=="tiff")    tiff(paste("Multi-traits_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                        if(file=="png") png(paste("Multi-traits_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                        if(!is.null(legend.ncol) && legend.pos=="middle"){
                            mar[3] = mar[3] + ceiling(length(trait) / legend.ncol)
                        }
                        par(mar=highlight_mar(mar),xaxs="i",yaxs="i")
                    }
                    if(!file.output){
                        ht=ifelse(is.null(height), 6, height)
                        wh=ifelse(is.null(width), 14, width)
                        if(is.null(dev.list())) dev.new(width=wh, height=ht)
                        # par(xpd=TRUE)
                    }
                    
                    # [OPTIMIZED] Flatten trait columns through a numeric matrix; as.vector(data.frame) returns a list.
                    pvalue <- as.numeric(as.matrix(Pmap[,3:(R+2)]))
                    if(is.null(ylim)){
                        if(!is.null(threshold)){
                            if(LOG10){
                                Max=round_y_axis_max(max_without_na(c((-log10(min_without_na(pvalue))),-log10(min_without_na(unlist(threshold))))))
                                Min<-round_y_axis_min(min_without_na(c((-log10(max_without_na(pvalue))),-log10(max_without_na(unlist(threshold))))))
                            }else{
                                Max=round_y_axis_max(max_without_na(c((max_without_na(pvalue)),max_without_na(unlist(threshold)))))
                                # if(abs(Max)<=1)   Max=max_without_na(c(max_without_na(pvalue),max_without_na(threshold)))
                                Min <- round_y_axis_min(min_without_na(c((min_without_na(pvalue)),min_without_na(unlist(threshold)))))
                                # if(abs(Min)<=1)   Min=min_without_na(c(min_without_na(pvalue),min_without_na(threshold)))
                            }
                        }else{
                            if(LOG10){
                                    Max=round_y_axis_max((-log10(min_without_na(pvalue))))
                                    Min=round_y_axis_min((-log10(max_without_na(pvalue))))
                            }else{
                                    Max=round_y_axis_max((max_without_na(pvalue)))
                                    # if(abs(Max)<=1)   Max=max_without_na(max_without_na(pvalue))
                                    Min<- round_y_axis_min((min_without_na(pvalue)))
                                    # if(abs(Min)<=1)   Min=min_without_na(min_without_na(pvalue))
                                    # }else{
                                        # Max=max_without_na(ceiling(max_without_na(pvalue)))
                            }
                        }
                        if((Max-Min)<=1){
                            if(cir.density){
                                plot(NULL,xlim=c(min_without_na(pvalue.posN)-band,band+1.05*max_without_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold, Max),ann=FALSE,
                                    cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }else{
                                plot(NULL,xlim=c(min_without_na(pvalue.posN)-band,band+max_without_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                    cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }
                        }else{
                            if(cir.density){
                                plot(NULL,xlim=c(min_without_na(pvalue.posN)-band,band+1.05*max_without_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold,Max),ann=FALSE,
                                    cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }else{
                                plot(NULL,xlim=c(min_without_na(pvalue.posN)-band,band+max_without_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                    cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }
                        }
                        mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                    }else{
                        Max <- max_without_na(unlist(ylim))
                        Min <- min_without_na(unlist(ylim))
                        if(cir.density){
                            plot(NULL,xlim=c(min_without_na(pvalue.posN)-band,band+1.05*max_without_na(pvalue.posN)),ylim=c(Min-Max/den.fold,Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                        }else{
                            plot(NULL,xlim=c(min_without_na(pvalue.posN)-band,band+max_without_na(pvalue.posN)),ylim=c(Min, Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                        }
                        mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                    }
        
                    # Max1 <- Max
                    # Min1 <- Min
                    # if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
                    # if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
                    if(length(unique(col)) == 1 && is.null(signal.col)) stop("'signal.col' is NULL.")
                    if(length(unique(col)) == 1 && amplify == FALSE)    stop("'amplify' is FALSE.")
                    legend_col <- t(col)[1:R]
                    if(length(unique(col)) == 1)    legend_col <- rep(signal.col, R)[1:R]
                    if(legend.pos=="middle"){
                        if(is.null(legend.ncol)){
                            legend((max_without_na(pvalue.posN)+min_without_na(pvalue.posN))*0.5,Max,trait,col=legend_col,pch=pch,text.font=6,cex=legend.cex,box.col=NA,horiz=TRUE,xjust=0.5,yjust=0,xpd=TRUE)
                        }else{
                            legend((max_without_na(pvalue.posN)+min_without_na(pvalue.posN))*0.5,Max,trait,col=legend_col,pch=pch,text.font=6,cex=legend.cex,box.col=NA,horiz=FALSE,ncol=legend.ncol,xjust=0.5,yjust=0,xpd=TRUE)
                        }
                    }else if(legend.pos=="left" || legend.pos=="right"){
                        if(is.null(legend.ncol)){
                            legend(ifelse(legend.pos=="left","topleft","topright"),trait,col=legend_col,pch=pch,text.font=6,cex=legend.cex,box.col=NA,horiz=FALSE,xpd=TRUE)
                        }else{
                            legend(ifelse(legend.pos=="left","topleft","topright"),trait,col=legend_col,pch=pch,text.font=6,cex=legend.cex,box.col=NA,horiz=FALSE,ncol=legend.ncol,xpd=TRUE)
                        }
                    }
        
                    draw_rect_chr_axis(scale=1, Min=Min, Max=Max, title=TRUE)
                    if(is.null(ylim)){
                        if((Max-Min)>1){
                            #print(seq(0,(Max+1),ceiling((Max+1)/10)))
                            axis(2,las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                            axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                            legend.y <- Max
                        }else{
                            axis(2,las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                            axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                            legend.y <- Max
                        }
                    }else{
                        axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                        axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                        legend.y <- Max
                    }
                    if(chr.border){
                        for(b in 1:length(chr.border.pos)){
                            segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=axis.lwd, lty=2)
                        }
                    }
        
                    if(length(unique(col)) != 1){
                        sam.index <- list()
                        trait_max_n <- 0
                        trait_max <- 0
                        for(l in 1:R){
                            sam.index[[l]] <- c(1:nrow(Pmap))[is_visable[[l]] & !is.na(logpvalueT[,l])]
                            if(length(sam.index[[l]]) >= trait_max_n){
                                trait_max_n=length(sam.index[[l]])
                                trait_max=l
                            }
                        }
                        
                        #change the sample number according to Pmap
                        #sam.num <- ceiling(nrow(Pmap)/100)
                        sam.num <- 1000
                        cat_bar <- seq(1, 100, 1)
                        trait_n <- sapply(sam.index, length)
                        trait_sams <- ceiling(trait_n / sam.num)
                        trait_max_sams <- max(trait_sams)
                        trait_1st_sam <- trait_max_sams - trait_sams + 1
                        trait_full_sams <- floor(trait_n / sam.num)
                        trait_1st_full_sam <- trait_max_sams - trait_full_sams + 1
                        for(sam in 1:trait_max_sams) {
                            for(i in 1:R){
                                if(sam < trait_1st_sam[i]){
                                    # nothing
                                }else{
                                    if(sam < trait_1st_full_sam[i]){
                                        plot.index <- sample(sam.index[[i]], trait_n[i] %% sam.num, replace=FALSE)
                                    }else{
                                        plot.index <- sample(sam.index[[i]], sam.num, replace=FALSE)
                                    }
                                    sam.index[[i]] <- sam.index[[i]][-which(sam.index[[i]] %in% plot.index)]
                                    logpvalue=logpvalueT[plot.index,i]
                                    if(!is.null(ylim)){indexx <- logpvalue>=min_without_na(ylim[[i]])}else{indexx <- 1:length(logpvalue)}
                                    points(pvalue.posN[plot.index][indexx],logpvalue[indexx],pch=pch[i],type=type,lwd=cex[2]+1,cex=cex[2],col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                                }
                            }
                            if(verbose){
                                progress <- round((nrow(Pmap) - length(sam.index[[trait_max]])) * 100 / nrow(Pmap))
                                if(progress %in% cat_bar){
                                    cat(" Multi-traits Rectangular plotting ... (finished ", progress, "%)\r", sep="")
                                    cat_bar <- cat_bar[cat_bar != progress]
                                    if(progress == 100) cat("\n")
                                }
                            }
                        }
                    }else{
                        for(i in 1:R){
                            logpvalue=logpvalueT[,i]
                            if(!is.null(ylim)){indexx <- logpvalue>=min_without_na(ylim[[i]])}else{indexx <- 1:length(logpvalue)}
                            points(pvalue.posN[indexx],logpvalue[indexx],pch=pch[i],type=type,lwd=cex[2]+1,cex=cex[2],col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                        }
                    }
                    if(!is.null(threshold)){
                        for(thr in 1:length(threshold[[i]])){
                            h <- threshold_to_y(threshold[[i]][thr])
                            segments(0, h, max_without_na(pvalue.posN), h, col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                        }
                        if(amplify==TRUE){
                            if(length(unique(col)) != 1){
                                for(i in 1:R){
                                    logpvalue=logpvalueT[, i]
                                    for(ll in 1:length(threshold[[i]])){
                                        if(ll == 1){
                                            if(LOG10){
                                                sgline1=threshold_to_y(threshold[[i]][ll])
                                            }else{
                                                sgline1=threshold_to_y(threshold[[i]][ll])
                                            }
                                            sgindex=which(logpvalue>=sgline1)
                                            HY1=logpvalue[sgindex]
                                            HX1=pvalue.posN[sgindex]
                                        }else{
                                            if(LOG10){
                                                sgline0=threshold_to_y(threshold[[i]][ll-1])
                                                sgline1=threshold_to_y(threshold[[i]][ll])
                                            }else{
                                                sgline0=threshold_to_y(threshold[[i]][ll-1])
                                                sgline1=threshold_to_y(threshold[[i]][ll])
                                            }
                                            sgindex=which(logpvalue>=sgline1 & logpvalue < sgline0)
                                            HY1=logpvalue[sgindex]
                                            HX1=pvalue.posN[sgindex]
                                        }
                                        points(HX1,HY1,pch=pch[i],cex=cex[2],col="white")
                                        if(is.null(signal.col)){
                                            points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                                        }else{
                                            points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rgb(t(col2rgb(signal.col[ll])), alpha=points.alpha, maxColorValue=255))
                                        }
                                        
                                    }
                                }
                            }else{
                                for(i in 1:R){
                                    logpvalue=logpvalueT[, i]
                                    if(LOG10){
                                        sgindex = which(logpvalue > -log10(min(unlist(threshold))))
                                    }else{
                                        sgindex = which(logpvalue > max(unlist(threshold)))
                                    }
                                    HY1=logpvalue[sgindex]
                                    HX1=pvalue.posN[sgindex]
                                    points(HX1,HY1,pch=pch[i],cex=cex[2],col="white")
                                    points(HX1,HY1,pch=rep(signal.pch, R)[i],cex=rep(signal.cex, R)[i],col=rgb(t(col2rgb(rep(signal.col, R)[i])), alpha=points.alpha, maxColorValue=255))
                                }
                            }
                        }
                    }
        
                    if(is.null(ylim)){ymin <- Min}else{ymin <- min_without_na(unlist(ylim))}
                    if(cir.density){
                        for(yll in 1:length(pvalue.posN.list)){
                            polygon(c(min_without_na(pvalue.posN.list[[yll]]), min_without_na(pvalue.posN.list[[yll]]), max_without_na(pvalue.posN.list[[yll]]), max_without_na(pvalue.posN.list[[yll]])), 
                                c(ymin-0.5*(Max-Min)/den.fold, ymin-1.5*(Max-Min)/den.fold, 
                                ymin-1.5*(Max-Min)/den.fold, ymin-0.5*(Max-Min)/den.fold), 
                                col="grey", border="grey")
                        }
                        is_visable_den <- filter_visible_points(pvalue.posN, ymin-0.5*(Max-Min)/den.fold, wh, ht, dpi=dpi)
                        segments(
                            pvalue.posN[is_visable_den],
                            ymin-0.5*(Max-Min)/den.fold,
                            pvalue.posN[is_visable_den],
                            ymin-1.5*(Max-Min)/den.fold,
                            col=density.list$den.col[is_visable_den], lwd=0.5
                        )
                        legend(
                            x=max_without_na(pvalue.posN)+band,
                            y=legend.y,
                            title="", legend=density.list$legend.y, pch=15, pt.cex=2.5, col=density.list$legend.col,
                            cex=legend.cex*0.8, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=0.9, xjust=0, xpd=TRUE
                        )          
                    }
                    if(!is.null(main))  title(main=main[1], cex.main=main.cex, font.main= main.font)
                    if(box) box(lwd=axis.lwd)
                    if(file.output) dev.off()
                }
            }else{
                #print("Starting Rectangular-Manhattan plot!",quote=F)
                if(!is.null(file.name) && length(file.name) != R)   stop(paste("please provide a vector containing file names of all", R, "traits."))
                for(i in 1:R){
                    colx=col[i,]
                    colx=colx[!is.na(colx)]
                    if(verbose) cat(paste(" Rectangular Manhattan plotting ",trait[i],".\n",sep=""))
                        if(file.output){
                            ht=ifelse(is.null(height), 6, height)
                            wh=ifelse(is.null(width), 14, width)
                            if(file=="jpg") jpeg(paste("Rect_Manhtn.",ifelse(is.null(file.name),trait[i],file.name[i]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                            if(file=="pdf") pdf(paste("Rect_Manhtn.",ifelse(is.null(file.name),trait[i],file.name[i]),".pdf",sep=""), width=wh,height=ht)
                            if(file=="tiff")    tiff(paste("Rect_Manhtn.",ifelse(is.null(file.name),trait[i],file.name[i]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                            if(file=="png") png(paste("Rect_Manhtn.",ifelse(is.null(file.name),trait[i],file.name[i]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                            par(mar=highlight_mar(mar),xaxs="i",yaxs="i")
                        }
                        if(!file.output){
                            ht=ifelse(is.null(height), 6, height)
                            wh=ifelse(is.null(width), 14, width)
                            if(is.null(dev.list())) dev.new(width=wh, height=ht)
                            # par(xpd=TRUE)
                        }
                        
                        pvalue=pvalueT[,i]
                        logpvalue=logpvalueT[,i]
                        if(is.null(ylim)){
                            y.candidates <- logpvalue[is.finite(logpvalue)]
                            if(!length(y.candidates)) stop("No finite values remain for Manhattan plotting after skip/cut filtering.")
                            if(!is.null(threshold[[i]])) y.candidates <- c(y.candidates, threshold_to_y(threshold[[i]]))
                            Max=round_y_axis_max(max_without_na(y.candidates))
                            Min<-round_y_axis_min(min_without_na(y.candidates))
                            if(!is.null(skip)) Min <- max(Min, skip)
                            if(Max <= Min) Max <- Min + 1
                            if((Max-Min)<=1){
                                if(cir.density){
                                    plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,band+1.05*max_without_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold, Max),ann=FALSE,
                                        cex.axis=axis.cex,font=lab.font,axes=FALSE)
                                }else{
                                    plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,band+max_without_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                    cex.axis=axis.cex,font=lab.font,axes=FALSE)
                                }
                            }else{
                                if(cir.density){
                                    plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,band+1.05*max_without_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold,Max),ann=FALSE,
                                    cex.axis=axis.cex,font=lab.font,axes=FALSE)
                                }else{
                                    plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,band+max_without_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                    cex.axis=axis.cex,font=lab.font,axes=FALSE)
                                }
                            }
                            mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                        }else{
                            Max <- max_without_na(ylim[[i]])
                            Min <- min_without_na(ylim[[i]])
                            if(cir.density){
                                plot(pvalue.posN[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,band+1.05*max_without_na(pvalue.posN)),ylim=c(min_without_na(ylim[[i]])-(Max-Min)/den.fold, max_without_na(ylim[[i]])),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }else{
                                plot(pvalue.posN[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_without_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_without_na(pvalue.posN)-band,band+max_without_na(pvalue.posN)),ylim=ylim[[i]],ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }
                            mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                        }
                        # Max1 <- Max
                        # Min1 <- Min
                        # if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
                        # if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
                        if(chr.border){
                            for(b in 1:length(chr.border.pos)){
                                segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=axis.lwd, lty=2)
                            }
                        }
        
                        draw_rect_chr_axis(scale=1, Min=Min, Max=Max, title=TRUE)
                        if(is.null(ylim)){
                            if((Max-Min)>1){
                                axis(2, las=1, lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                                axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                                legend.y <- Max
                            }else{
                                axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                                axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                                legend.y <- Max
                            }
                        }else{
                            axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                            axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                            legend.y <- tail(ylim[[i]][2], 1)
                        }
                        if(!is.null(threshold[[i]])){
                            for(thr in 1:length(threshold[[i]])){
                                h <- threshold_to_y(threshold[[i]][thr])
                                # print(h)
                                # print(threshold.col[thr])
                                # print(threshold.lty[thr])
                                # print(threshold.lwd[thr])
                                segments(0, h, max_without_na(pvalue.posN), h,col=threshold.col[thr],lty=threshold.lty[thr],lwd=threshold.lwd[thr])
                            }
                            if(amplify == TRUE){
                                if(LOG10){
                                    threshold[[i]] <- sort(threshold[[i]])
                                    sgline1=threshold_to_y(max_without_na(threshold[[i]]))
                                }else{
                                    threshold[[i]] <- sort(threshold[[i]], decreasing=TRUE)
                                    sgline1=threshold_to_y(min_without_na(threshold[[i]]))
                                }
        
                                sgindex=which(logpvalue>=sgline1)
                                HY1=logpvalue[sgindex]
                                HX1=pvalue.posN[sgindex]
                                
                                #cover the points that exceed the threshold with the color "white"
                                points(HX1,HY1,pch=pch,cex=cex[2],col="white")
                                
                                for(ll in 1:length(threshold[[i]])){
                                    if(ll == 1){
                                        if(LOG10){
                                            sgline1=threshold_to_y(threshold[[i]][ll])
                                        }else{
                                            sgline1=threshold_to_y(threshold[[i]][ll])
                                        }
                                        sgindex=which(logpvalue>=sgline1)
                                        HY1=logpvalue[sgindex]
                                        HX1=pvalue.posN[sgindex]
                                    }else{
                                        if(LOG10){
                                            sgline0=threshold_to_y(threshold[[i]][ll-1])
                                            sgline1=threshold_to_y(threshold[[i]][ll])
                                        }else{
                                            sgline0=threshold_to_y(threshold[[i]][ll-1])
                                            sgline1=threshold_to_y(threshold[[i]][ll])
                                        }
                                        sgindex=which(logpvalue>=sgline1 & logpvalue < sgline0)
                                        HY1=logpvalue[sgindex]
                                        HX1=pvalue.posN[sgindex]
                                    }
        
                                    if(is.null(signal.col)){
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rep(rep(colx,N[i]),add[[i]])[sgindex])
                                    }else{
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=signal.col[ll])
                                    }
                                    
                                }
                            }
        
                        }
        
                        if(!is.null(highlight)){
                            # points(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],pch=pch,cex=cex[2],col="white")
                            if(!is.na(highlight_index[[i]][1])){
                                if(is.null(highlight.col)){
                                    draw_highlight_labels(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_without_na(pvalue.posN)-band,max_without_na(pvalue.posN)+band),ylim=c(Min,Max),words=highlight.text[[i]],point.cex=highlight.cex,text.cex=highlight.text.cex, pch=highlight.pch,type=highlight.type,point.col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]],text.col=highlight.text.col,text.font=highlight.text.font,top.cex=highlight.text.top.cex,mode=highlight.text.mode,side=highlight.text.side,nearby.offset=highlight.text.nearby.offset,nearby.step=highlight.text.nearby.step,line.col=highlight.text.line.col,line.lwd=highlight.text.line.lwd,line.lty=highlight.text.line.lty,line.mode=highlight.text.line.mode,line.bend=highlight.text.line.bend,arrow=highlight.text.arrow,arrow.length=highlight.text.arrow.length,top.space=highlight.text.top.space,top.inside=highlight.text.top.inside,optimize=highlight.text.optimize,lanes=highlight.text.lanes,lane.gap=highlight.text.lane.gap,min.gap=highlight.text.min.gap)
                                }else{
                                    draw_highlight_labels(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_without_na(pvalue.posN)-band,max_without_na(pvalue.posN)+band),ylim=c(Min,Max),words=highlight.text[[i]],point.cex=highlight.cex,text.cex=highlight.text.cex, pch=highlight.pch,type=highlight.type,point.col=highlight_col[[i]],text.col=highlight.text.col,text.font=highlight.text.font,top.cex=highlight.text.top.cex,mode=highlight.text.mode,side=highlight.text.side,nearby.offset=highlight.text.nearby.offset,nearby.step=highlight.text.nearby.step,line.col=highlight.text.line.col,line.lwd=highlight.text.line.lwd,line.lty=highlight.text.line.lty,line.mode=highlight.text.line.mode,line.bend=highlight.text.line.bend,arrow=highlight.text.arrow,arrow.length=highlight.text.arrow.length,top.space=highlight.text.top.space,top.inside=highlight.text.top.inside,optimize=highlight.text.optimize,lanes=highlight.text.lanes,lane.gap=highlight.text.lane.gap,min.gap=highlight.text.min.gap)
                                }
                            }
                        }
        
                        #if(!is.null(threshold) & !is.null(signal.line))    abline(v=pvalue.posN[which(pvalueT[,i] < min_without_na(threshold))],col="grey",lty=2,lwd=signal.line)
                
                        if(is.null(ylim)){ymin <- Min}else{ymin <- min_without_na(ylim[[i]])}
                        if(cir.density){
                            for(yll in 1:length(pvalue.posN.list)){
                                polygon(c(min_without_na(pvalue.posN.list[[yll]]), min_without_na(pvalue.posN.list[[yll]]), max_without_na(pvalue.posN.list[[yll]]), max_without_na(pvalue.posN.list[[yll]])), 
                                    c(ymin-0.5*(Max-Min)/den.fold, ymin-1.5*(Max-Min)/den.fold, 
                                    ymin-1.5*(Max-Min)/den.fold, ymin-0.5*(Max-Min)/den.fold), 
                                    col="grey", border="grey", xpd=TRUE)
                            }
                            is_visable_den <- filter_visible_points(pvalue.posN, ymin-0.5*(Max-Min)/den.fold, wh, ht, dpi=dpi)
                            segments(
                                pvalue.posN[is_visable_den],
                                ymin-0.5*(Max-Min)/den.fold,
                                pvalue.posN[is_visable_den],
                                ymin-1.5*(Max-Min)/den.fold,
                                col=density.list$den.col[is_visable_den], lwd=0.5,xpd=TRUE
                            )
                            legend(
                                x=max_without_na(pvalue.posN)+band,
                                y=legend.y,
                                title="", legend=density.list$legend.y, pch=15, pt.cex=2.5, col=density.list$legend.col,
                                cex=legend.cex*0.8, bty="n",
                                y.intersp=1,
                                x.intersp=1,
                                yjust=0.9, xjust=0, xpd=TRUE
                            )
                            
                        }
                    if(!is.null(main))  title(main=main[i], cex.main=main.cex, font.main= main.font)
                    if(box) box(lwd=axis.lwd)
                    if(file.output)  dev.off()
                }
            }
        }
    }, envir=env)
}


# Circular Manhattan plot branch used by CMplot().
# This function draws circular single- and multi-trait Manhattan plots.

draw_circular_manhattan_plot <- function(env) {
    # Run the plot.type == "c" branch inside the CMplot() runtime environment.
    evalq({
        #plot circle Manhattan
        if("c" %in% plot.type){
        
            signal.line.index <- NULL
            if(!is.null(threshold)){
                if(!is.null(signal.line)){
                    for(l in 1:R){
                        if(!is.null(threshold[[l]])){
                            if(LOG10){
                                signal.line.index <- c(signal.line.index,which(pvalueT[,l] < min_without_na(threshold[[l]])))
                            }else{
                                signal.line.index <- c(signal.line.index,which(pvalueT[,l] > max_without_na(threshold[[l]])))
                            }
                        }
                    }
                    signal.line.index <- unique(signal.line.index)
                }
                signal.line.index <- pvalue.posN[signal.line.index]
            }
        
            if(file.output){
                ht=ifelse(is.null(height), 10, height)
                wh=ifelse(is.null(width), 10, width)
                if(file=="jpg") jpeg(paste("Cir_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                if(file=="pdf") pdf(paste("Cir_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht)
                if(file=="tiff")    tiff(paste("Cir_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                if(file=="png") png(paste("Cir_Manhtn.",ifelse(is.null(file.name),taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
            }
            if(!file.output){
                ht=ifelse(is.null(height), 10, height)
                wh=ifelse(is.null(width), 10, width)
                if(is.null(dev.list())) dev.new(width=wh, height=ht)
                par(pty="s", xpd=TRUE)
            }
            RR <- r+H*R+cir.band*R
            # [OPTIMIZED] Pre-compute sin/cos for circular plots (used 50+ times)
            angle_factor <- 2 * base::pi / TotalN
            pvalue_angles <- (pvalue.posN - round(band/2) - circleMin) * angle_factor
            sin_pvalue_angles <- sin(pvalue_angles)
            cos_pvalue_angles <- cos(pvalue_angles)
            if(cir.density){
                plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
            }else{
                plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
            }
            if(!is.null(signal.line)){
                if(!is.null(signal.line.index)){
                    X1chr <- (RR)*sin(2*base::pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                    Y1chr <- (RR)*cos(2*base::pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                    X2chr <- (r)*sin(2*base::pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                    Y2chr <- (r)*cos(2*base::pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                    segments(X1chr,Y1chr,X2chr,Y2chr,lty=2,lwd=signal.line,col="grey")
                }
            }
            for(i in 1:R){
            
                #get the colors for each trait
                colx <- col[i,]
                colx <- colx[!is.na(colx)]
        
                if(verbose) cat(paste(" Circular Manhattan plotting ",trait[i],".\n",sep=""))
                pvalue <- pvalueT[,i]
                logpvalue <- logpvalueT[,i]
                if(is.null(ylim)){
                    if(LOG10){
                        Max <- round_y_axis_max(-log10(min_without_na(pvalue)))
                        Min <- round_y_axis_min(-log10(max_without_na(pvalue)))
                    }else{
                        Max <- round_y_axis_max(max_without_na(pvalue))
                        #if(abs(Max)<=1)    Max <- max_without_na(pvalue)
                        Min <- round_y_axis_min(min_without_na(pvalue))
                        #if(abs(Min)<=1)    Min <- min_without_na(pvalue)
                    }
                }else{
                    Max <- ylim[[i]][2]
                    Min <- ylim[[i]][1]
                }
                Cpvalue <- (H*(logpvalue-Min))/(Max-Min)
                ylimIndx <- logpvalue >= Min & logpvalue <= Max
                if(outward==TRUE){
                    if(cir.chr==TRUE & i == 1){
                        
                        #plot the boundary which represents the chromosomes
                        polygon.num <- 1000
                        for(k in 1:length(chr)){
                            if(k==1){
                                polygon.index <- seq(round(band/2)+1,-round(band/2)-circleMin+max_without_na(pvalue.posN.list[[1]]), length=polygon.num)
                                #change the axis from right angle into circle format
                                X1chr=(RR)*sin(2*base::pi*(polygon.index)/TotalN)
                                Y1chr=(RR)*cos(2*base::pi*(polygon.index)/TotalN)
                                X2chr=(RR+cir.chr.h)*sin(2*base::pi*(polygon.index)/TotalN)
                                Y2chr=(RR+cir.chr.h)*cos(2*base::pi*(polygon.index)/TotalN)
                                if(is.null(chr.den.col)){
                                    polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])   
                                }else{
                                    if(cir.density){
                                            polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                    }else{
                                            polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                    }
                                }
                            }else{
                                polygon.index <- seq(1+round(band/2)+max_without_na(pvalue.posN.list[[k-1]])-circleMin,-round(band/2)-circleMin+max_without_na(pvalue.posN.list[[k]]), length=polygon.num)
                                X1chr=(RR)*sin(2*base::pi*(polygon.index)/TotalN)
                                Y1chr=(RR)*cos(2*base::pi*(polygon.index)/TotalN)
                                X2chr=(RR+cir.chr.h)*sin(2*base::pi*(polygon.index)/TotalN)
                                Y2chr=(RR+cir.chr.h)*cos(2*base::pi*(polygon.index)/TotalN)
                                if(is.null(chr.den.col)){
                                    polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])
                                }else{
                                    if(cir.density){
                                            polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                    }else{
                                            polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                    }
                                }       
                            }
                        }
                        
                        if(cir.density){
        
                            if(file.output){
                                is_visable <- filter_visible_points((RR+cir.chr.h)*sin_pvalue_angles, (RR+cir.chr.h)*cos_pvalue_angles, wh, ht, dpi=dpi)
                            }else{
                                is_visable <- rep(TRUE, length(pvalue.posN))
                            }
                            segments(
                                (RR)*sin_pvalue_angles[is_visable],
                                (RR)*cos_pvalue_angles[is_visable],
                                (RR+cir.chr.h)*sin_pvalue_angles[is_visable],
                                (RR+cir.chr.h)*cos_pvalue_angles[is_visable],
                                col=density.list$den.col[is_visable], lwd=0.5
                            )
                            legend(
                                x=RR+4*cir.chr.h,
                                y=(RR+4*cir.chr.h)/2,
                                title="", legend=density.list$legend.y, pch=15, pt.cex=3, col=density.list$legend.col,
                                cex=legend.cex, bty="n",
                                y.intersp=1,
                                x.intersp=1,
                                yjust=0.3, xjust=0, xpd=TRUE
                            )
                            
                        }
                        
                        # XLine=(RR+cir.chr.h)*sin(2*base::pi*(1:TotalN)/TotalN)
                        # YLine=(RR+cir.chr.h)*cos(2*base::pi*(1:TotalN)/TotalN)
                        # lines(XLine,YLine,lwd=1.5)
                        if(cir.density){
                            draw_circle(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
                            draw_circle(myr=RR,lwd=1.5,add=TRUE,col='grey')
                        }else{
                            draw_circle(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
                            draw_circle(myr=RR,lwd=1.5,add=TRUE)
                        }
        
                    }
                    
                    X=(Cpvalue[ylimIndx]+r+H*(i-1)+cir.band*(i-1))*sin_pvalue_angles[ylimIndx]
                    Y=(Cpvalue[ylimIndx]+r+H*(i-1)+cir.band*(i-1))*cos_pvalue_angles[ylimIndx]
                    if(file.output){
                        is_visable <- filter_visible_points(X, Y, wh, ht, dpi=dpi)
                    }else{
                        is_visable <- rep(TRUE, length(X))
                    }
        
                    if(cir.axis && cir.axis.grid){
                        draw_circle(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                        draw_circle(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                        draw_circle(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                        draw_circle(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                        draw_circle(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    }
        
                    points(X[is_visable],Y[is_visable],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]])[ylimIndx][is_visable])
                    
                    #plot the legend for each trait
                    if(cir.axis==TRUE){
                        #try to get the number after radix point
                        if((Max-Min) > 1) {
                            round.n=2
                        }else{
                            if(Max == 1){
                                round.n=1
                            }else{
                                round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
                            }
                        }
                        segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
        
                        lab=seq(round(Min+(Max-Min)*0,round.n), round(Min+(Max-Min)*1,round.n), length=5)
                        text(-H/20,r+H*(i-0.94)+cir.band*(i-1),lab[1],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                        text(-H/20,r+H*(i-0.75)+cir.band*(i-1),lab[2],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                        text(-H/20,r+H*(i-0.5)+cir.band*(i-1),lab[3],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                        text(-H/20,r+H*(i-0.25)+cir.band*(i-1),lab[4],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                        text(-H/20,r+H*(i-0.06)+cir.band*(i-1),lab[5],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    }
                    
                    if(!is.null(threshold[[i]])){
                        if(sum(threshold[[i]]!=0)==length(threshold[[i]])){
                            for(thr in 1:length(threshold[[i]])){
                                significantline1=ifelse(LOG10, H*(-log10(threshold[[i]][thr])-Min)/(Max-Min), H*(threshold[[i]][thr]-Min)/(Max-Min))
                                #s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*base::pi*(0:TotalN)/TotalN)
                                #s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*base::pi*(0:TotalN)/TotalN)
                                if(significantline1<H){
                                    #lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
                                    draw_circle(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                                }else{
                                    warning(paste("No significant points for ",trait[i]," pass the threshold level using threshold=",threshold[[i]][thr],"!",sep=""))
                                }
                            }
                        }
                    }
                    
                    if(!is.null(threshold[[i]])){
                        if(sum(threshold[[i]]!=0)==length(threshold[[i]])){
                            if(amplify==TRUE){
                                if(LOG10){
                                    threshold[[i]] <- sort(threshold[[i]])
                                    significantline1=H*(-log10(max_without_na(threshold[[i]]))-Min)/(Max-Min)
                                }else{
                                    threshold[[i]] <- sort(threshold[[i]], decreasing=TRUE)
                                    significantline1=H*(min_without_na(threshold[[i]])-Min)/(Max-Min)
                                }
                                
                                p_amp.index <- which(Cpvalue>=significantline1)
                                HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin_pvalue_angles[p_amp.index]
                                HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos_pvalue_angles[p_amp.index]
        
                                #cover the points that exceed the threshold with the color "white"
                                points(HX1,HY1,pch=19,cex=cex[1],col="white")
        
                                for(ll in 1:length(threshold[[i]])){
                                    if(ll == 1){
                                        if(LOG10){
                                            significantline1=H*(-log10(threshold[[i]][ll])-Min)/(Max-Min)
                                        }else{
                                            significantline1=H*(threshold[[i]][ll]-Min)/(Max-Min)
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1)
                                        HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin_pvalue_angles[p_amp.index]
                                        HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos_pvalue_angles[p_amp.index]
                                    }else{
                                        if(LOG10){
                                            significantline0=H*(-log10(threshold[[i]][ll-1])-Min)/(Max-Min)
                                            significantline1=H*(-log10(threshold[[i]][ll])-Min)/(Max-Min)
                                        }else{
                                            significantline0=H*(threshold[[i]][ll-1]-Min)/(Max-Min)
                                            significantline1=H*(threshold[[i]][ll]-Min)/(Max-Min)
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
                                        HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin_pvalue_angles[p_amp.index]
                                        HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos_pvalue_angles[p_amp.index]
                                    }
        
                                    if(is.null(signal.col)){
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
                                    }else{
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=signal.col[ll])
                                    }
                                }
                            }
                        }
                    }
        
                    if(!is.null(highlight)){
                        HX1=(Cpvalue[highlight_index[[i]]]+r+H*(i-1)+cir.band*(i-1))*sin_pvalue_angles[highlight_index[[i]]]
                        HY1=(Cpvalue[highlight_index[[i]]]+r+H*(i-1)+cir.band*(i-1))*cos_pvalue_angles[highlight_index[[i]]]
                        points(HX1,HY1,pch=19,cex=cex[1],col="white")
                        if(is.null(highlight.col)){
                            points(HX1,HY1,pch=highlight.pch,cex=highlight.cex,col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]])
                        }else{
                            points(HX1,HY1,pch=highlight.pch,cex=highlight.cex,col=highlight_col[[i]])
                        }
                    }
        
                    if(cir.chr==TRUE){
                        ticks1=(RR+1.5*cir.chr.h)*sin(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                        ticks2=(RR+1.5*cir.chr.h)*cos(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                        if(is.null(chr.labels)){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=lab.font,cex=lab.cex-0.5, adj=c(0.5, 0))
                            }
                        }else{
                            if(Nchr == 1){
                                for(t in 1:length(ticks)){
                                    angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                    text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab, sep=""),srt=angle, adj=c(0.5, 0),font=lab.font,cex=lab.cex-0.5)
                                }
                            }else{
                                for(t in 1:length(ticks)){
                                    angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                    text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=lab.font,cex=lab.cex-0.5, adj=c(0.5, 0))
                                }
                            }
                        }
                    }else{
                        ticks1=1.01*RR*sin(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                        ticks2=1.01*RR*cos(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                        # ticks1=(0.9*r)*sin(2*base::pi*(ticks-round(band/2))/TotalN)
                        # ticks2=(0.9*r)*cos(2*base::pi*(ticks-round(band/2))/TotalN)
                        if(is.null(chr.labels)){
                            for(t in 1:length(ticks)){
                            angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                            text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }else{
                            if(Nchr == 1){
                                for(t in 1:length(ticks)){
                                    angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                    text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab, sep=""),srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                                }
                            }else{
                                for(t in 1:length(ticks)){
                                    angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                    text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                                }
                            }
                        }
                    }
                }
                if(outward==FALSE){
                    if(cir.chr==TRUE & i == 1){
                        # XLine=(2*cir.band+RR+cir.chr.h)*sin(2*base::pi*(1:TotalN)/TotalN)
                        # YLine=(2*cir.band+RR+cir.chr.h)*cos(2*base::pi*(1:TotalN)/TotalN)
                        # lines(XLine,YLine,lwd=1.5)
        
                        polygon.num <- 1000
                        for(k in 1:length(chr)){
                            if(k==1){
                                polygon.index <- seq(round(band/2)+1,-round(band/2)-circleMin+max_without_na(pvalue.posN.list[[1]]), length=polygon.num)
                                X1chr=(RR)*sin(2*base::pi*(polygon.index)/TotalN)
                                Y1chr=(RR)*cos(2*base::pi*(polygon.index)/TotalN)
                                X2chr=(RR+cir.chr.h)*sin(2*base::pi*(polygon.index)/TotalN)
                                Y2chr=(RR+cir.chr.h)*cos(2*base::pi*(polygon.index)/TotalN)
                                    if(is.null(chr.den.col)){
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])   
                                    }else{
                                        if(cir.density){
                                            polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                        }else{
                                            polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                        }
                                    }
                            }else{
                                polygon.index <- seq(1+round(band/2)+max_without_na(pvalue.posN.list[[k-1]])-circleMin,-round(band/2)-circleMin+max_without_na(pvalue.posN.list[[k]]), length=polygon.num)
                                X1chr=(RR)*sin(2*base::pi*(polygon.index)/TotalN)
                                Y1chr=(RR)*cos(2*base::pi*(polygon.index)/TotalN)
                                X2chr=(RR+cir.chr.h)*sin(2*base::pi*(polygon.index)/TotalN)
                                Y2chr=(RR+cir.chr.h)*cos(2*base::pi*(polygon.index)/TotalN)
                                if(is.null(chr.den.col)){
                                    polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])   
                                }else{
                                        if(cir.density){
                                            polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                        }else{
                                            polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                        }
                                }   
                            }
                        }
                        if(cir.density){
        
                            if(file.output){
                                is_visable <- filter_visible_points((RR+cir.chr.h)*sin_pvalue_angles, (RR+cir.chr.h)*cos_pvalue_angles, wh, ht, dpi=dpi)
                            }else{
                                is_visable <- rep(TRUE, length(pvalue.posN))
                            }
                            segments(
                                (RR)*sin_pvalue_angles[is_visable],
                                (RR)*cos_pvalue_angles[is_visable],
                                (RR+cir.chr.h)*sin_pvalue_angles[is_visable],
                                (RR+cir.chr.h)*cos_pvalue_angles[is_visable],
                                col=density.list$den.col[is_visable], lwd=0.5
                            )
                            legend(
                                x=RR+4*cir.chr.h,
                                y=(RR+4*cir.chr.h)/2,
                                title="", legend=density.list$legend.y, pch=15, pt.cex=3, col=density.list$legend.col,
                                cex=legend.cex, bty="n",
                                y.intersp=1,
                                x.intersp=1,
                                yjust=0.3, xjust=0, xpd=TRUE
                            )
                            
                        }
                        
                        if(cir.density){
                            draw_circle(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
                            draw_circle(myr=RR,lwd=1.5,add=TRUE,col='grey')
                        }else{
                            draw_circle(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
                            draw_circle(myr=RR,lwd=1.5,add=TRUE)
                        }
        
                    }
        
                    X=(-Cpvalue[ylimIndx]+r+H*i+cir.band*(i-1))*sin_pvalue_angles[ylimIndx]
                    Y=(-Cpvalue[ylimIndx]+r+H*i+cir.band*(i-1))*cos_pvalue_angles[ylimIndx]
                    if(file.output){
                        is_visable <- filter_visible_points(X, Y, wh, ht, dpi=dpi)
                    }else{
                        is_visable <- rep(TRUE, length(X))
                    }
        
                    if(cir.axis && cir.axis.grid){
                        draw_circle(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                        draw_circle(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                        draw_circle(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                        draw_circle(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                        draw_circle(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    }
        
                    points(X[is_visable],Y[is_visable],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]])[ylimIndx][is_visable])
                    
                    if(cir.axis==TRUE){
                        
                        #try to get the number after radix point
                        if((Max-Min)<=1) {
                            if(Max == 1){
                                round.n=1
                            }else{
                                round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
                            }
                        }else{
                            round.n=2
                        }
                        segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                        
                        lab=seq(round(Min+(Max-Min)*0,round.n), round(Min+(Max-Min)*1,round.n), length=5)
                        text(-H/20,r+H*(i-0.06)+cir.band*(i-1),lab[1],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                        text(-H/20,r+H*(i-0.25)+cir.band*(i-1),lab[2],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                        text(-H/20,r+H*(i-0.5)+cir.band*(i-1),lab[3],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                        text(-H/20,r+H*(i-0.75)+cir.band*(i-1),lab[4],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                        text(-H/20,r+H*(i-0.94)+cir.band*(i-1),lab[5],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    }
                    
                    if(!is.null(threshold[[i]])){
                        if(sum(threshold[[i]]!=0)==length(threshold[[i]])){
                        
                            for(thr in 1:length(threshold[[i]])){
                                significantline1=ifelse(LOG10, H*(-log10(threshold[[i]][thr])-Min)/(Max-Min), H*(threshold[[i]][thr]-Min)/(Max-Min))
                                #s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
                                #s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
                                if(significantline1<H){
                                    #lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
                                    draw_circle(myr=(-significantline1+r+H*i+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                                }else{
                                    warning(paste("No significant points for ",trait[i]," pass the threshold level using threshold=",threshold[[i]][thr],"!",sep=""))
                                }
                            }
                            if(amplify==TRUE){
                                if(LOG10){
                                    threshold[[i]] <- sort(threshold[[i]])
                                    significantline1=H*(-log10(max_without_na(threshold[[i]]))-Min)/(Max-Min)
                                }else{
                                    threshold[[i]] <- sort(threshold[[i]], decreasing=TRUE)
                                    significantline1=H*(min_without_na(threshold[[i]])-Min)/(Max-Min)
                                }
                                p_amp.index <- which(Cpvalue>=significantline1)
                                HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin_pvalue_angles[p_amp.index]
                                HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos_pvalue_angles[p_amp.index]
        
                                #cover the points that exceed the threshold with the color "white"
                                points(HX1,HY1,pch=19,cex=cex[1],col="white")
        
                                    for(ll in 1:length(threshold[[i]])){
                                        if(ll == 1){
                                            if(LOG10){
                                                significantline1=H*(-log10(threshold[[i]][ll])-Min)/(Max-Min)
                                            }else{
                                                significantline1=H*(threshold[[i]][ll]-Min)/(Max-Min)
                                            }
                                            p_amp.index <- which(Cpvalue>=significantline1)
                                            HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin_pvalue_angles[p_amp.index]
                                            HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos_pvalue_angles[p_amp.index]
                                        }else{
                                            if(LOG10){
                                                significantline0=H*(-log10(threshold[[i]][ll-1])-Min)/(Max-Min)
                                                significantline1=H*(-log10(threshold[[i]][ll])-Min)/(Max-Min)
                                            }else{
                                                significantline0=H*(threshold[[i]][ll-1]-Min)/(Max-Min)
                                                significantline1=H*(threshold[[i]][ll]-Min)/(Max-Min)
                                            }
                                            p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
                                            HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin_pvalue_angles[p_amp.index]
                                            HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos_pvalue_angles[p_amp.index]
        
                                        }
        
                                        if(is.null(signal.col)){
                                            points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
                                        }else{
                                            points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=signal.col[ll])
                                        }
                                    }
                            }
                        }
                    }
        
                    if(!is.null(highlight)){
                        HX1=(-Cpvalue[highlight_index[[i]]]+r+H*i+cir.band*(i-1))*sin_pvalue_angles[highlight_index[[i]]]
                        HY1=(-Cpvalue[highlight_index[[i]]]+r+H*i+cir.band*(i-1))*cos_pvalue_angles[highlight_index[[i]]]
                        points(HX1,HY1,pch=19,cex=cex[1],col="white")
                        if(is.null(highlight.col)){
                            points(HX1,HY1,pch=highlight.pch,cex=highlight.cex,col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]])
                        }else{
                            points(HX1,HY1,pch=highlight.pch,cex=highlight.cex,col=highlight_col[[i]])
                        }
                    }
        
                    if(cir.chr==TRUE){
                        ticks1=(RR+1.5*cir.chr.h)*sin(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                        ticks2=(RR+1.5*cir.chr.h)*cos(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                        if(is.null(chr.labels)){
                            for(t in 1:length(ticks)){
                              angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                              text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }else{
                            if(Nchr == 1){
                                for(t in 1:length(ticks)){
                                    angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                    text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab,sep=""),srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                                }
                            }else{
                                for(t in 1:length(ticks)){
                                    angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                    text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                                }
                            }
                        }
                    }else{
                        ticks1=1.01*RR*sin(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                        ticks2=1.01*RR*cos(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                        # ticks1=RR*sin(2*base::pi*(ticks-round(band/2))/TotalN)
                        # ticks2=RR*cos(2*base::pi*(ticks-round(band/2))/TotalN)
                        if(is.null(chr.labels)){
                            for(t in 1:length(ticks)){
                            
                                #adjust the angle of labels of circle plot
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }else{
                            if(Nchr == 1){
                                for(t in 1:length(ticks)){
                                    angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                    text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab,sep=""),srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                                }
                            }else{
                                for(t in 1:length(ticks)){
                                    angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                    text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                                }
                            }
                        }   
                    }
                }
            }
            if(file.output) dev.off()
            #print("Circular-Manhattan has been finished!",quote=F)
        }
    }, envir=env)
}


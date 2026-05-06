# Q-Q plot branches used by CMplot().
# This function draws single-trait, multi-track, and multi-trait Q-Q plots.

draw_qq_plots <- function(env) {
    # Run the plot.type == "q" branch inside the CMplot() runtime environment.
    evalq({
        if("q" %in% plot.type){
        
            signal.col <- rep(signal.col,R)
            signal.pch <- rep(signal.pch,R)
            signal.cex <- rep(signal.cex*1.1,R)
        
            if(multracks | multraits){
                if(R < 2)   stop("need more than one trait.")
                if(multracks){
                    if(file.output){
                        ht=ifelse(is.null(height), 5.5, height)
                        wh=ifelse(is.null(width), 3.5, width)
                        if(file=="jpg") jpeg(paste("Multi-tracks_QQplot.",ifelse(is.null(file.name),taxa,file.name[1]),".jpg",sep=""), width=R*wh*dpi,height=ht*dpi,res=dpi,quality=100)
                        if(file=="pdf") pdf(paste("Multi-tracks_QQplot.",ifelse(is.null(file.name),taxa,file.name[1]),".pdf",sep=""), width=R*wh,height=ht)
                        if(file=="tiff")    tiff(paste("Multi-tracks_QQplot.",ifelse(is.null(file.name),taxa,file.name[1]),".tiff",sep=""), width=R*wh*dpi,height=ht*dpi,res=dpi)
                        if(file=="png") png(paste("Multi-tracks_QQplot.",ifelse(is.null(file.name),taxa,file.name[1]),".png",sep=""), width=R*wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                        par(mfcol=c(1,R),xpd=TRUE)
                    }else{
                        ht=ifelse(is.null(height), 5.5, height)
                        wh=ifelse(is.null(width), 3.5, width)
                        if(is.null(dev.list())) dev.new(width=wh*R, height=ht)
                        par(xpd=TRUE)
                    }
                    for(i in 1:R){
                        if(i == 1)  par(mar=c(mar[2], mar[2], mar[3], 0))
                        if(i == R)  par(mar=c(mar[2], 1.5, mar[3], mar[4]))
                        if(i != 1 & i != R) par(mar=c(mar[2], 1.5, mar[3], 0))
                        if(verbose) cat(paste(" Multi-tracks Q-Q plotting ",trait[i],".\n",sep=""))        
                        qq.data <- prepare_qq_data(i)
                        N <- qq.data$N
                        log.Quantiles <- qq.data$log.Quantiles
                        log.P.values <- qq.data$log.P.values
                        
                        #calculate the confidence interval of QQ-plot
                        # [OPTIMIZED] Vectorized qbeta() instead of for-loop (~2.6x faster)
                        if(conf.int){
                            N1=length(log.Quantiles)
                            xi <- ceiling((10^-log.Quantiles) * N)
                            xi[xi == 0] <- 1
                            c95 <- qbeta(0.95, xi, N - xi + 1)
                            c05 <- qbeta(0.05, xi, N - xi + 1)
                            index=length(c95):1
                        }else{
                            c05 <- 1
                            c95 <- 1
                        }
                        
                        YlimMax <- max_without_na(c(floor(max_without_na(c(max_without_na(-log10(c05)), max_without_na(-log10(c95))))+1), floor(max_without_na(log.P.values)+1)))
                        if(is.null(ylim)){
                            plot(NULL, xlim=c(0,floor(max_without_na(log.Quantiles)+1)), axes=FALSE, cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0,YlimMax),xlab ="", ylab="")
                        }else{
                            plot(NULL, xlim=c(0,floor(max_without_na(log.Quantiles)+1)), axes=FALSE, cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0,max(ylim[[i]])),xlab ="", ylab="")
                        }
                        axis(1, mgp=c(3,xticks.pos,0), at=seq(0,floor(max_without_na(log.Quantiles)+1),ceiling((max_without_na(log.Quantiles)+1)/10)), lwd=axis.lwd,labels=seq(0,floor(max_without_na(log.Quantiles)+1),ceiling((max_without_na(log.Quantiles)+1)/10)), cex.axis=axis.cex)
                        axis(2, las=1, lwd=axis.lwd,cex.axis=axis.cex)
                        axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(ylim[[i]]))), labels=c("",""), tcl=0, lwd=axis.lwd)
                        
                        #plot the confidence interval of QQ-plot
                        if(conf.int){
                            if(is.null(conf.int.col)){
                                polygon(c(log.Quantiles[index],log.Quantiles),c(plot_y(-log10(c05))[index],plot_y(-log10(c95))),col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                            }else{
                                polygon(c(log.Quantiles[index],log.Quantiles),c(plot_y(-log10(c05))[index],plot_y(-log10(c95))),col=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255))
                            }
                        }
                        if(!is.null(threshold.col)){par(xpd=FALSE); lines(log.Quantiles, plot_y(log.Quantiles), lwd=threshold.lwd[1], lty=threshold.lty[1], col=threshold.col[1]); par(xpd=TRUE)}
                        is_visable <- filter_visible_points(log.Quantiles, log.P.values, wh, ht, dpi=dpi)
                        if(!is.null(threshold[[i]])){
                            # if(sum(threshold!=0)==length(threshold)){
                                thre.line=threshold_to_y(min_without_na(threshold[[i]]))
                                if(amplify==TRUE){
                                    thre.index <- is.finite(log.P.values) & log.P.values<thre.line
                                    if(sum(!thre.index)!=0){
                                        points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col=t(col)[i],pch=19,cex=cex[3])
                                    
                                        #cover the points that exceed the threshold with the color "white"
                                        # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                        if(is.null(signal.col)){
                                            points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=t(col)[i],pch=signal.pch[i],cex=signal.cex[i])
                                        }else{
                                            points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=signal.col[i],pch=signal.pch[i],cex=signal.cex[i])
                                        }
                                    }else{
                                        points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                                    }
                                }else{
                                    points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                                }
                            # }
                        }else{
                            points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                        }
                        mtext(side=1, text=expression(Expected~~-log[10](italic(p))), line=ylab.pos+2, cex=lab.cex, font=lab.font, xpd=TRUE)
                        if(i == 1)  mtext(side=2, text=expression(Observed~~-log[10](italic(p))), line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                        if(!is.null(main)) {
                            title(main=main[i], cex.main=main.cex, font.main= main.font)
                        }else{
                            title(main=trait[i], cex.main=main.cex, font.main= main.font) 
                        }
                        if(box) box(lwd=axis.lwd)
                    }
                    if(file.output) dev.off()
                }
                if(multraits){
                    signal.col <- NULL
                    log.Quantiles.max_without_na <- NULL
                    for(i in 1:R){
                        P.values=as.numeric(Pmap[,i+2])
                        P.values=P.values[!is.na(P.values)]
                        p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                        log.Quantiles <- -log10(p_value_quantiles)
                        log.Quantiles.max_without_na <- c(log.Quantiles.max_without_na, max_without_na(log.Quantiles))
                    }
                    if(file.output){
                        ht=ifelse(is.null(height), 5.5, height)
                        wh=ifelse(is.null(width), 5.5, width)
                        if(file=="jpg") jpeg(paste("Multi-traits_QQplot.",ifelse(is.null(file.name),taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                        if(file=="pdf") pdf(paste("Multi-traits_QQplot.",ifelse(is.null(file.name),taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht)
                        if(file=="tiff")    tiff(paste("Multi-traits_QQplot.",ifelse(is.null(file.name),taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                        if(file=="png") png(paste("Multi-traits_QQplot.",ifelse(is.null(file.name),taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                        par(mar=c(mar[2],mar[2],mar[3],mar[4]),xpd=TRUE)
                    }else{  
                        ht=ifelse(is.null(height), 5.5, height)
                        wh=ifelse(is.null(width), 5.5, width)
                        dev.new(width=wh, height=ht)
                        par(xpd=TRUE)
                    }
                    p_value_quantiles=(1:nrow(Pmap))/(nrow(Pmap)+1)
                    log.Quantiles <- -log10(p_value_quantiles)
                                                
                    # calculate the confidence interval of QQ-plot
                    # [OPTIMIZED] Vectorized qbeta() instead of for-loop (~2.6x faster)
                    if(conf.int){
                        N1=length(log.Quantiles)
                        xi <- ceiling((10^-log.Quantiles) * N1)
                        xi[xi == 0] <- 1
                        c95 <- qbeta(0.95, xi, N1 - xi + 1)
                        c05 <- qbeta(0.05, xi, N1 - xi + 1)
                        index=length(c95):1
                    }
                    
                    if(!conf.int){c05 <- 1; c95 <- 1}
                    
                    if(is.null(ylim)){
                        # [OPTIMIZED] Use a numeric vector for multi-trait QQ ylim calculations.
                        Pmap.min_without_na <- as.numeric(as.matrix(Pmap[,3:(R+2)]))
                        YlimMax <- max_without_na(c(floor(max_without_na(c(max_without_na(-log10(c05)), max_without_na(-log10(c95))))+1), -log10(min_without_na(Pmap.min_without_na[Pmap.min_without_na > 0]))))
                        plot(NULL, xlim=c(0,floor(max_without_na(log.Quantiles.max_without_na)+1)), axes=FALSE, xlab="", ylab="", cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0, floor(YlimMax+1)), main = "QQplot", cex.main=main.cex, font.main=main.font)
                    }else{
                        plot(NULL, xlim=c(0,floor(max_without_na(log.Quantiles.max_without_na)+1)), axes=FALSE, xlab="", ylab="", cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0, max(unlist(ylim))),main = "QQplot", cex.main=main.cex, font.main=main.font)
                    }
                    legend("topleft",trait,col=rgb(t(col2rgb(t(col)[1:R])), alpha=points.alpha, maxColorValue=255),pch=19,cex=legend.cex,text.font=6,box.col=NA, xpd=TRUE)
                    axis(1, mgp=c(3,xticks.pos,0), at=seq(0,floor(max_without_na(log.Quantiles.max_without_na)+1),ceiling((max_without_na(log.Quantiles.max_without_na)+1)/10)), lwd=axis.lwd,labels=seq(0,floor(max_without_na(log.Quantiles.max_without_na)+1),ceiling((max_without_na(log.Quantiles.max_without_na)+1)/10)), cex.axis=axis.cex)
                    axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex)
                    axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(unlist(ylim)))), labels=c("",""), tcl=0, lwd=axis.lwd)
        
                    mtext(side=1, text=expression(Expected~~-log[10](italic(p))), line=ylab.pos+1, cex=lab.cex, font=lab.font, xpd=TRUE)
                    mtext(side=2, text=expression(Observed~~-log[10](italic(p))), line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                    
                    for(i in 1:R){
                        if(verbose) cat(paste(" Multi-traits Q-Q plotting ",trait[i],".\n",sep=""))
                        qq.data <- prepare_qq_data(i)
                        N <- qq.data$N
                        log.Quantiles <- qq.data$log.Quantiles
                        log.P.values <- qq.data$log.P.values
        
                        #calculate the confidence interval of QQ-plot
                        # [OPTIMIZED] Vectorized qbeta() instead of for-loop (~2.6x faster)
                        if(conf.int){
                            N1=length(log.Quantiles)
                            xi <- ceiling((10^-log.Quantiles) * N)
                            xi[xi == 0] <- 1
                            c95 <- qbeta(0.95, xi, N - xi + 1)
                            c05 <- qbeta(0.05, xi, N - xi + 1)
                            index=length(c95):1
                        }else{
                            c05 <- 1
                            c95 <- 1
                        }
        
                        # plot the confidence interval of QQ-plot
                        if(conf.int){
                            if(is.null(conf.int.col)){
                                polygon(c(log.Quantiles[index],log.Quantiles),c(plot_y(-log10(c05))[index],plot_y(-log10(c95))),col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                            }else{
                                polygon(c(log.Quantiles[index],log.Quantiles),c(plot_y(-log10(c05))[index],plot_y(-log10(c95))),col=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255))
                            }
                        }
                           
                        if((i == R) & !is.null(threshold.col)){par(xpd=FALSE); lines(log.Quantiles, plot_y(log.Quantiles), lwd=threshold.lwd[1], lty=threshold.lty[1], col=threshold.col[1]); par(xpd=TRUE)}
                        # points(log.Quantiles, log.P.values, col=t(col)[i],pch=19,cex=cex[3])
                        is_visable <- filter_visible_points(log.Quantiles, log.P.values, wh, ht, dpi=dpi)
                        if(!is.null(threshold[[i]])){
                            # if(sum(threshold!=0)==length(threshold)){
                                thre.line=threshold_to_y(min_without_na(threshold[[i]]))
                                if(amplify==TRUE){
                                    thre.index <- is.finite(log.P.values) & log.P.values<thre.line
                                    if(sum(!thre.index)!=0){
                                        points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=19,cex=cex[3])
                                
                                        # cover the points that exceed the threshold with the color "white"
                                        # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                        if(is.null(signal.col)){
                                            points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=signal.pch[i],cex=signal.cex[i])
                                        }else{
                                            points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=rgb(t(col2rgb(signal.col[i])), alpha=points.alpha, maxColorValue=255),pch=signal.pch[i],cex=signal.cex[i])
                                        }
                                    }else{
                                        points(log.Quantiles[is_visable], log.P.values[is_visable], col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=19,cex=cex[3])
                                    }
                                }else{
                                    points(log.Quantiles[is_visable], log.P.values[is_visable], col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=19,cex=cex[3])
                                }
                            # }
                        }else{
                            points(log.Quantiles[is_visable], log.P.values[is_visable], col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=19,cex=cex[3])
                        }
                    }
                    if(!is.null(main)) {
                        title(main=main[1], cex.main=main.cex, font.main= main.font)
                    }
                    if(box) box(lwd=axis.lwd)
                    if(file.output) dev.off()
                }
            }else{
                if(!is.null(file.name) && length(file.name) != R)   stop(paste("please provide a vector containing file names of all", R, "traits."))
                for(i in 1:R){
                    if(verbose) cat(paste(" Q-Q plotting ",trait[i],".\n",sep=""))
                    if(file.output){
                        ht=ifelse(is.null(height), 5.5, height)
                        wh=ifelse(is.null(width), 5.5, width)
                        if(file=="jpg") jpeg(paste("QQplot.",ifelse(is.null(file.name),trait[i],file.name[i]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                        if(file=="pdf") pdf(paste("QQplot.",ifelse(is.null(file.name),trait[i],file.name[i]),".pdf",sep=""), width=wh,height=ht)
                        if(file=="tiff") tiff(paste("QQplot.",ifelse(is.null(file.name),trait[i],file.name[i]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                        if(file=="png") png(paste("QQplot.",ifelse(is.null(file.name),trait[i],file.name[i]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                        par(mar=c(mar[2],mar[2],mar[3],mar[4]),xpd=TRUE)
                    }else{
                        ht=ifelse(is.null(height), 5.5, height)
                        wh=ifelse(is.null(width), 5.5, width)
                        if(is.null(dev.list())) dev.new(width=wh, height=ht)
                        par(xpd=TRUE)
                    }
                    qq.data <- prepare_qq_data(i)
                    N <- qq.data$N
                    log.Quantiles <- qq.data$log.Quantiles
                    log.P.values <- qq.data$log.P.values
                    
                    #calculate the confidence interval of QQ-plot
                    # [OPTIMIZED] Vectorized qbeta() instead of for-loop (~2.6x faster)
                    if(conf.int){
                        N1=length(log.Quantiles)
                        xi <- ceiling((10^-log.Quantiles) * N)
                        xi[xi == 0] <- 1
                        c95 <- qbeta(0.95, xi, N - xi + 1)
                        c05 <- qbeta(0.05, xi, N - xi + 1)
                        index=length(c95):1
                    }else{
                        c05 <- 1
                        c95 <- 1
                    }
                    if(is.null(ylim)){
                        YlimMax <- max_without_na(c(floor(max_without_na(c(max_without_na(-log10(c05)), max_without_na(-log10(c95))))+1), floor(max_without_na(log.P.values)+1)))
                        plot(NULL, xlim=c(0,floor(max_without_na(log.Quantiles)+1)), axes=FALSE, cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0,YlimMax),xlab="",ylab="")
                    }else{
                        plot(NULL, xlim=c(0,floor(max_without_na(log.Quantiles)+1)), axes=FALSE, cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0,max(ylim[[i]])),xlab="",ylab="")      
                    }
                    axis(1, mgp=c(3,xticks.pos,0),at=seq(0,floor(max_without_na(log.Quantiles)+1),ceiling((max_without_na(log.Quantiles)+1)/10)), lwd=axis.lwd,labels=seq(0,floor(max_without_na(log.Quantiles)+1),ceiling((max_without_na(log.Quantiles)+1)/10)), cex.axis=axis.cex)
                    axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex)
                    axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(ylim[[i]]))), labels=c("",""), tcl=0, lwd=axis.lwd)
        
                    mtext(side=1, text=expression(Expected~~-log[10](italic(p))), line=ylab.pos+1, cex=lab.cex, font=lab.font, xpd=TRUE)
                    mtext(side=2, text=expression(Observed~~-log[10](italic(p))), line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                    
                    #plot the confidence interval of QQ-plot
                    if(conf.int){
                        if(is.null(conf.int.col)){
                            polygon(c(log.Quantiles[index],log.Quantiles),c(plot_y(-log10(c05))[index],plot_y(-log10(c95))),col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                        }else{
                            polygon(c(log.Quantiles[index],log.Quantiles),c(plot_y(-log10(c05))[index],plot_y(-log10(c95))),col=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255))
                        }
                    }
        
                    if(!is.null(threshold.col)){par(xpd=FALSE); lines(log.Quantiles, plot_y(log.Quantiles), lwd=threshold.lwd[1], lty=threshold.lty[1], col=threshold.col[1]); par(xpd=TRUE)}
                    # points(log.Quantiles, log.P.values, col=t(col)[i],pch=19,cex=cex[3])
                    is_visable <- filter_visible_points(log.Quantiles, log.P.values, wh, ht, dpi=dpi)
                    if(!is.null(threshold[[i]])){
                        # if(sum(threshold!=0)==length(threshold)){
                            thre.line=threshold_to_y(min_without_na(threshold[[i]]))
                            if(amplify==TRUE){
                                thre.index <- is.finite(log.P.values) & log.P.values<thre.line
                                if(sum(!thre.index)!=0){
                                    points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col=t(col)[i],pch=19,cex=cex[3])
                                
                                    #cover the points that exceed the threshold with the color "white"
                                    # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                    # print(signal.col)
                                    # print(signal.pch)
                                    # print(signal.cex)
                                    if(is.null(signal.col)){
                                        points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=t(col)[i],pch=signal.pch[i],cex=signal.cex[i])
                                    }else{
                                        points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=signal.col[i],pch=signal.pch[i],cex=signal.cex[i])
                                    }
                                }else{
                                    points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                                }
                            }else{
                                points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                            }
                        # }
                    }else{
                        points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                    }
                    if(!is.null(main)) {
                        title(main=main[i], cex.main=main.cex, font.main= main.font)
                    }else{
                        title(main=trait[i], cex.main=main.cex, font.main= main.font) 
                    }
                    if(box) box(lwd=axis.lwd)
                    if(file.output) dev.off()
                }
            }
        }
    }, envir=env)
}


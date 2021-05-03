#!/usr/bin/env Rscript


# variable
#red0 <- rgb(1,0,0,0.7)
#blue0 <- rgb(0,1,0,0.7)
#purple0 <- rgb(160/255,32/255,240/255,0.7)
#green0 <- rgb(0,0,1,0.7)
#yellow0 <- rgb(1,1,0,0.7)
#gray0 <- rgb(0.7,0.7,0.7,0.7)
#grey0 <- gray0
#black0 <- rgb(1,1,1,0.7)


#png.0 <- function(filename='Rplot.png',width=1200,height=1200,unit='px',pointsize=22,...){
png.0 <- function(filename='Rplot.png',width=5,height=5,unit='in',res=300,pointsize=10,...){
    return(png(filename=filename,width=width,height=height,unit=unit,pointsize=pointsize,res=res,...))
}

blackred <- function(n){
    colorpanel(n, "black", "red")
}


yellowblue <- function(n){
    colorpanel(n, "yellow", "blue")
}

na.fill <- function(x, f=mean){
    y <- x
    y[is.na(y)] <- f(y,na.rm=T)
    y
}


mean.0 <- function(...){
    mean(..., na.rm=T)
}

sd.0 <- function(...){
    sd(..., na.rm=T)
}


se.0 <- function(x,...){
    sd.0(x,...) / sqrt(sum(!is.na(x))-1)
}

se <- function(x, ...){
    sd(x,...) / sqrt(sum(!is.na(x))-1)
}


min.0 <- function(...){
    min(..., na.rm=T)
}

max.0 <- function(...){
    max(..., na.rm=T)
}

cor.0 <- function(...){
    cor(..., use='pairwise.complete.obs')
}

cor.test.mat <- function(M){

    p <- matrix(0, ncol(M), ncol(M))
    colnames(p) <- colnames(M)
    rownames(p) <- colnames(M)

    for (i in 1:ncol(M)){
        for (j in 1:ncol(M)){
          m <- M[,i] + M[,j]
          if (sum(is.na(m)) > length(m)-1) p[i,j] <- NA
          else p[i,j] <- cor.test.0(M[,i],M[,j])$p.value
        }
    }

    return(p)
}

cor.test.0 <- function(...){
    cor.test(..., use='pairwise.complete.obs')
}

plotExprHMDP <- function(inds,cex=0.7,...){
    for (ind in inds){
        plot.0(datExpr[c(inds_ctrl,inds_iso),ind],
               col=c(rep('gray',length(inds_ctrl)),
                     rep('red',length(inds_iso))),
               xlab='Strains',
               ylab='Expression (log2)',
               pch=19,
               cex=cex,
               main=ifelse(is.na(symbols[ind]),probes[ind],symbols[ind]),
               ...
               )
    }
}

plotExprHMDP_hist <- function(inds_genes, inds_strains=inds_ctrl, breaks=15,lwd=3,col='red',...){
    for (ind in inds_genes){
        X <- datExpr[inds_strains,ind]
        hist.0(X,
               breaks=breaks,
               xlab='Expression (log2)',
               main=ifelse(is.na(symbols[ind]),probes[ind],symbols[ind]),
               freq=F,
               ...
               )
        lines(density(X),
              col=col,
              lwd=lwd,
              ...
              )
    }
}

plotExprHumanHF <- function(inds,...){
    for (ind in inds){
        plot.0(datExpr[c(inds_ctrl,inds_cmp,inds_isc),ind],
               col=c(rep('gray',length(inds_ctrl)),
                     rep('green',length(inds_cmp)),
                     rep('red',length(inds_isc))),
               xlab='Patients',
               ylab='Expression (log2)',
               main=ifelse(is.na(symbols[ind]),probes[ind],symbols[ind]),
               ...
               )
    }
}

plotExprHumanHF_hist <- function(inds_genes, inds_patients=inds_ctrl, breaks=15,lwd=3,col='red',...){
    for (ind in inds_genes){
        X <- datExpr[inds_patients,ind]
        hist.0(X,
               breaks=breaks,
               xlab='Expression (log2)',
               main=ifelse(is.na(symbols[ind]),probes[ind],symbols[ind]),
               freq=F,
               ...
               )
        lines(density(X),
              col=col,
              lwd=lwd,
              ...
              )
    }
}


plotExprHumanHF_miRNAs <- function(inds,...){
    for (ind in inds){
        plot.0(datExpr_miRNAs[c(inds_ctrl,inds_cmp,inds_isc),ind],
               col=c(rep('gray',length(inds_ctrl)),
                     rep('green',length(inds_cmp)),
                     rep('red',length(inds_isc))),
               xlab='Patients',
               ylab='Expression (log2)',
               main=miRNAs[ind],
               ...
               )
    }
}

plotExprRA <- function(inds,...){
    for (ind in inds){
        plot.0(datExpr[c(inds_nonresponder,inds_moderate,inds_responder),ind],
               col=c(rep('gray',length(inds_nonresponder)),
                     rep('green',length(inds_moderate)),
                     rep('red',length(inds_responder))),
               xlab='Patients',
               ylab='Expression (log2)',
               main=probes[ind],
               ...
               )
    }
}

boxplotExprRA <- function(inds,...){
    for (ind in inds){
        M <- list('nonresponder'=datExpr[inds_nonresponder,ind],
                   'moderate'=datExpr[inds_moderate,ind],
                   'responder'=datExpr[inds_responder,ind])
        boxplot.0(M,
                  col=c('gray','green','red'),
                  ylab='Expression (log2)',
                  mar=c(12,5,5,5),
                  las=3,
                  main=probes[ind],
                  ...
                  )
    }
}


plotExprRA_baseline <- function(inds,...){
    inds_nonresponder_baseline <- intersect(inds_baseline, inds_nonresponder)
    inds_responder_baseline <- intersect(inds_baseline, inds_responder)
    inds_moderate_baseline <- intersect(inds_baseline, inds_moderate)
    inds_all <- c(inds_nonresponder_baseline,inds_moderate_baseline,inds_responder_baseline)
    inds_abcon <- inds_all
    inds_abcon[!(inds_abcon %in% which(datPheno$Study=='Abcon'))] <- NA
    inds_b13 <- inds_all
    inds_b13[!(inds_b13 %in% which(datPheno$Study=='Batter-up 2013'))] <- NA
    inds_b14 <- inds_all
    inds_b14[!(inds_b14 %in% which(datPheno$Study=='Batter-up 2014'))] <- NA
    inds_jap <- inds_all
    inds_jap[!(inds_jap %in% which(datPheno$Study=='Japanese'))] <- NA
    for (ind in inds){
        plot.0(datExpr[inds_all,ind],
               col=c(rep('gray',length(inds_nonresponder_baseline)),
                     rep('green',length(inds_moderate_baseline)),
                     rep('red',length(inds_responder_baseline))),
               xlab='Patients',
               ylab='Expression (log2)',
               main=ifelse(is.na(symbols[ind]),probes[ind],symbols[ind]),
               ...
               )
        points(datExpr[inds_abcon,ind],pch=2)
        points(datExpr[inds_b13,ind],pch=3)
        points(datExpr[inds_b14,ind],pch=5)
        points(datExpr[inds_jap,ind],pch=7)
    }
}


plotExprRA_after <- function(inds,...){
    for (ind in inds){
        plot.0(datExpr[c(inds_nonresponder_after,inds_moderate_after,inds_responder_after),ind],
               col=c(rep('gray',length(inds_nonresponder_after)),
                     rep('green',length(inds_moderate_after)),
                     rep('red',length(inds_responder_after))),
               xlab='Patients',
               ylab='Expression (log2)',
               main=probes[ind],
               ...
               )
    }
}


plot.0 <- function(...,mar=c(6,6,5,5), cex.lab=1.3, cex.axis=1.3, 
                   cex=1.3, cex.main=1.3, pch=19){
    par(mar=mar, 
        cex.lab=cex.lab, 
        cex.axis=cex.axis, 
        cex.main=cex.main, 
        pch = pch)
    plot(..., cex=cex)
}

plot.cor <- function(x, y,method='p',cut_by=0.1,main=NULL,cex=0.3,cex.main=1.2,col='gray',
                     fitline=TRUE,log='',xlab='x',ylab='y',
                     sfrac=0.002, lwd=3,pch=19,gap=0, ...){

    require(gplots)

    #x[!is.finite(x)] <- NA
    #y[!is.finite(y)] <- NA

    sp_x <- split(x, cut(x, unique(quantile(x, probs = seq(0, 1, by = cut_by), na.rm=T)))) 
    sp_y <- split(y, cut(x, unique(quantile(x, probs = seq(0, 1, by = cut_by), na.rm=T)))) 
    xm <- unlist(lapply(sp_x,mean, na.rm=T))
    ym <- unlist(lapply(sp_y,mean, na.rm=T))
    xm_sd <- unlist(lapply(sp_x,function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x))-1)))
    ym_sd <- unlist(lapply(sp_y,function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x))-1)))

    if (fitline==TRUE){

        if (log=='x'){

            logx <- log(xm)
            inds <- which(is.na(logx) | !is.finite(logx))
            logx[inds] <- NA

            if (length(inds) == length(x)) return()
            
            fit <- lm(ym~logx, na.action = na.exclude)

        } else if (log=='y'){

            logy <- log(ym)
            inds <- which(is.na(logy) | !is.finite(logy))
            logy[inds] <- NA
            
            if (length(inds) == length(x)) return()
            
            fit <- lm(logy~xm, na.action = na.exclude)

        } else if (log=='xy'){

            logx <- log(xm)
            inds <- which(is.na(logx) | !is.finite(logx))
            logy <- log(ym)
            inds <- union(inds, which(is.na(logy) | !is.finite(logy)))
            
            logx[inds] <- NA
            logy[inds] <- NA
            
            if (length(inds) == length(xm)) return()

            fit <- lm(logy~logx, na.action = na.exclude)

        } else {
            fit <- lm(ym~xm)
        }
    }
    cc = cor.test.0(x,y,method=method)
    r = cc$estimate
    p = cc$p.value

    if (is.null(main)){
        if (method == 's'){
            main <- paste0('Spearman\'s rho = ',signif(r,2),'',
                           ' (p=',signif(p,2),')')
        } else if (method == 'p'){
            main <- paste0('Pearson\'s r = ',signif(r,2),'',
                           ' (p=',signif(p,2),')')
        }

        if (fitline == TRUE){
            main <- paste0(main, '\nFit coefficient = ', signif(coef(fit)[2], 2)) 
        }
    }

    plot.0(x,y,
           cex=cex,
           cex.main=cex.main,
           main=main,
           col=col,
           xlab = xlab,
           ylab = ylab,
           log=log,
           ...
           )
    if (fitline==TRUE){
        #abline(fit, untf=T)
        if (log == 'x'){
            lines(xm, predict(fit, newdata=list(log(xm))))
        } else if (log == 'y'){
            lines(xm, exp(predict(fit, newdata=list(xm))))
        } else if (log == 'xy'){
            lines(xm, exp(predict(fit, newdata=list(log(xm)))))
        } else {
            lines(xm, predict(fit, newdata=list(xm)))
        }
    }
    plotCI(xm,
           ym,
           uiw=ym_sd,
           sfrac=sfrac,
           lwd=lwd,
           cex=cex,
           pch=pch,
           add=T,
           gap=gap
           )
    plotCI(xm,
           ym,
           uiw=xm_sd,
           err='x',
           sfrac=sfrac,
           lwd=lwd,
           cex=cex,
           pch=pch,
           add=T,
           gap=gap
           )

    return()
}


plotCI.1 <- function(x,ymat,se=F,ylim=NA,col='black',add=F,lwd=1,pch=19,...){
    require(gplots)
    mus <- apply(ymat,2,mean,na.rm=T)
    if (se==T){
        sds <- apply(ymat,2,function(x) sd(x,na.rm=T)/sqrt(sum(!is.na(x))-1))
    } else{
        sds <- apply(ymat,2,function(x) sd(x,na.rm=T))
    }

    if (length(ylim)!=2){
        ymin <- min(mus-sds)
        ymax <- max(mus+sds)
        ylim <- c(ymin, ymax)
    }

    if (add==F){
        plot.0(x, mus, ylim=ylim, col=col,lwd=lwd,...)
    } else{
        lines(x, mus, col=col, lwd=lwd,...)
        points(x, mus, col=col, pch=pch,...)
    }
    plotCI(x, mus, 
           uiw=sds,
           sfrac=0.005,
           lwd=lwd,
           col=col,
           cex=0,
           gap=0,
           add=T)
    
}

plotCI.0 <- function(x, mus,uiw,ylim=NA,col='black',add=T,cex=0,type='p',lwd=1,pch=19,gap=0,sfrac=0.005,liw=uiw,...){
    require(gplots)

    if (length(ylim)!=2){
        ymin <- min(mus-liw)
        ymax <- max(mus+uiw)
        ylim <- c(ymin, ymax)
    }

    plotCI(x, mus, 
           uiw=uiw,
           liw=liw,
           sfrac=sfrac,
           lwd=lwd,
           type=type,
           pch=pch,
           cex=cex,
           col=col,
           ylim=ylim,
           gap=gap,
           add=add,
           ...)
}

beanplot.0 <- function(l, las=3, mar=c(6,6,5,5),bw="nrd",col=list('gray'),
                       cex.lab=1.3, cex.axis=1.3, cex=1.3, cex.main=1.3, pch = 19, ...){
    require(beanplot)
    par(mar=mar, 
        cex.lab=cex.lab, 
        cex.axis=cex.axis, 
        cex=cex, 
        cex.main=cex.main, 
        pch = pch)
    beanplot(l,
             what = c(0,1,1,0),
             bw=bw,
             las=las,
             cex = 0,
             col=col,
             ...
             ) 
}

barplot.0 <- function(l, las=3, mar=c(6,6,5,5),border='NA',col='gray',
                      cex.lab=1.3, cex.axis=1.3, cex=1.3, cex.main=1.3, ...){
    par(mar=mar, 
        cex.lab=cex.lab, 
        cex.axis=cex.axis, 
        cex=cex, 
        cex.main=cex.main)

    b <- barplot(l,
                 border=border,
                 las=las,
                 col=col,
                 ...
                 ) 

    return(b)
}

barplot.split <- function(x,y,ylim='',case='character',f='mean',...){
    if (is.numeric(y)) 
        y[!is.finite(y)] <- NA
    if (is.numeric(x))   
        x[!is.finite(x)] <- NA
    sp1 <- split(y,x)
    if (f == 'mean'){
      mu1 <- unlist(lapply(sp1,mean.0))
    } else if (f == 'median'){
      mu1 <- unlist(lapply(sp1,median,na.rm=T))
    }
    
    sd1 <- unlist(lapply(sp1,function(x) se.0(x)))
    if (case == 'character'){
    inds <- order(as.numeric(names(mu1)))
    } else {
    inds <- order(names(mu1))
    }
    if (length(ylim)==1){
        if (min.0(mu1)>0){
            ylim <- c(0,max.0(mu1+sd1))
        } else{
            ylim <- c(min.0(mu1-sd1),max.0(mu1+sd1))
        }
    }
    b <- barplot.0(mu1[inds],
                   ylim=ylim,
                   ...)
    plotCI.0(b, mu1[inds], uiw=ifelse(mu1[inds]<0,NA,sd1[inds]), liw=ifelse(mu1[inds]>0,NA,sd1[inds]),...)
}


setpar <- function(mar=c(6,6,5,5),cex.lab=1.3,cex.axis=1.3,
                   cex=1.3,cex.main=1.5,pch=19,...){
    par(mar=mar, 
        cex.lab=cex.lab, 
        cex.axis=cex.axis, 
        cex=cex, 
        cex.main=cex.main, 
        pch = pch,
        ...)
}

## functions
rgb.0 <- function(color_name, alpha=0.7){
    col <- col2rgb(color_name) / 255
    return(rgb(col[1],col[2],col[3],alpha))
}

change_extension <- function(file, newext){
    "Change extension of file. DON'T include dot."

    filenoext <- sub("\\.[^.]*$", "", file)
    return(paste(filenoext,newext,sep='.'))
}

hist2 <- function(dat1, dat2, nbreaks=50, log.breaks=F,main="", xlab="Value", ylab="", legend1="Data 1", legend2="Data 2",legend.pos='topleft',xlim=NULL,ylim=NULL,cex.axis=1.3,cex.lab=1.3,cex.leg=1.3,cex.main=1.3,col1=rgb(0,0,1,1/4),col2 = rgb(1,0,0,1/4), border=F, freq=FALSE,bty='n',...){
    "Produces two superimposed histograms for two sets of data,
    with alpha color blending."

    par(mar=c(6,6,5,5))

    if (log.breaks==TRUE){
        breaks <- 10^seq(from=log10(min(dat1,dat2)),
                         to=log10(max(dat1,dat2)),
                         length.out=nbreaks)
    } else {
        breaks <- seq(from=min(dat1,dat2),
                      to=max(dat1,dat2),
                      length.out=nbreaks
                      )
    }
    h1=hist(dat1,breaks=breaks, plot = FALSE)
    h2=hist(dat2,breaks=breaks, plot = FALSE)
    if (!length(xlim)) xlim <- c(min(h1$mids,h2$mids), max(h1$mids, h2$mids))
    if (!length(ylim)){
        if (freq==FALSE){
            ylim <- c(min(h1$density,h2$density), max(h1$density, h2$density))
        } else {
            ylim <- c(min(h1$counts,h2$counts), max(h1$counts, h2$counts))
        }
    }

    hist(dat1,
         breaks=breaks,
         freq = freq,
         col=col1,
         xlim = xlim,  
         ylim = ylim,  
         main=main,
         ylab=ylab,
         xlab=xlab,
         cex.main = cex.main,
         cex.lab = cex.lab,
         cex.axis = cex.axis, 
         border=border,
         ...
         )
    hist(dat2, 
         breaks=breaks,
         freq = freq,         
         col=col2, 
         add=TRUE,
         border=border,
         ...
         )
    legend(legend.pos, 
           legend = c(legend1, legend2), 
           border=border,
           bty=bty,
           fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)),
           cex = cex.leg) 
}

hist2.arrow <- function(dat1, dat2, x_arrows, col_arrows=rep('red',length(x_arrows)),
                        fac=4, nbreaks=50, log.breaks=F, main="", 
                        xlab="Value", ylab="", legend1="Data 1", 
                        legend2="Data 2", legend.pos='topleft',
                        xlim = c(0.9 * min(dat1,dat2, x_arrows), 
                                 1.1 * max(dat1,dat2, x_arrows)),
                        ylim=NULL, cex.axis=1.3, cex.lab=1.3,
                        cex.leg=1.3,cex.main=1.3, lwd=6,
                        col1=rgb(0,0,1,1/4),mar=c(6,6,5,5),
                        col2 = rgb(1,0,0,1/4), 
                        border=F, freq=FALSE,bty='n',...){
    "Produces two superimposed histograms for two sets of data,
    with alpha color blending."

    par(mar=mar)

    if (log.breaks==TRUE){
        breaks <- 10^seq(from=log10(min(dat1,dat2)),
                         to=log10(max(dat1,dat2)),
                         length.out=nbreaks)
    } else {
        breaks <- seq(from=min(dat1,dat2),
                      to=max(dat1,dat2),
                      length.out=nbreaks
                      )
    }
    h1=hist(dat1,breaks=breaks, plot = FALSE)
    h2=hist(dat2,breaks=breaks, plot = FALSE)
    #if (!length(xlim)) xlim <- c(min(h1$mids,h2$mids), max(h1$mids, h2$mids))
    if (!length(ylim)){
        if (freq==FALSE){
            ylim <- c(min(h1$density,h2$density), max(h1$density, h2$density))
        } else {
            ylim <- c(min(h1$counts,h2$counts), max(h1$counts, h2$counts))
        }
    }

    hist(dat1,
         breaks=breaks,
         freq = freq,
         col=col1,
         xlim = xlim,  
         ylim = ylim,  
         main=main,
         ylab=ylab,
         xlab=xlab,
         cex.main = cex.main,
         cex.lab = cex.lab,
         cex.axis = cex.axis, 
         border=border,
         ...
         )
    hist(dat2, 
         breaks=breaks,
         freq = freq,         
         col=col2, 
         add=TRUE,
         border=border,
         ...
         )
    legend(legend.pos, 
           legend = c(legend1, legend2), 
           border=border,
           bty=bty,
           fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)),
           cex = cex.leg) 

    for (j in 1:length(x_arrows)){
        x0 <- x_arrows[j]
        arrows(x0=x0,
               y0=max(h1$density,h2$density)/fac,
               y1=0,
               lwd=lwd,
               col=col_arrows[j]
               )
        text(x=x0,
             y=max(h1$density,h2$density)/fac,
             labels=names(x_arrows)[j],
             pos=3,
             cex=1.2,
             col=col_arrows[j]
             )
    }
}


hist2.loglog <- function(dat1, dat2, nbreaks=100, log.breaks=TRUE, main="", xlab="Value", ylab="Frequency", legend1="Data 1", legend2="Data 2",xlim=NULL,ylim=NULL,cex.lab=1.5,cex.axis=1.4,cex=1.2,col1=rgb(0,0,1,3/4),col2 = rgb(1,0,0,3/4)){
    "Produces log plots of the histograms for two sets of data,
    with alpha color blending."

    par(mar=c(6,6,5,5))
    dat1 <- dat1[dat1>0]
    dat2 <- dat2[dat2>0]
    if (log.breaks==TRUE){
        breaks <- 10^seq(from=log10(min(dat1,dat2)),
                         to=log10(max(dat1,dat2)),
                         length.out=nbreaks)
    } else {
        breaks <- seq(from=min(dat1,dat2),
                      to=max(dat1,dat2),
                      length.out=nbreaks
                      )
    }
    
    h1 <- hist(dat1,breaks=breaks, plot = FALSE)
    h2 <- hist(dat2,breaks=breaks, plot = FALSE)
    if (!length(xlim)){
        xlim <- log10(c(min(h1$mids,h2$mids), max(h1$mids, h2$mids)))
    }
    if (!length(ylim)){
        ylim <- log10(c(min(h1$density[h1$density>0],h2$density[h2$density>0]), max(h1$density, h2$density)))
    }

    plot(log10(h1$mids),
         log10(h1$density),
         pch=16,
         col=col1,
         xlab=xlab,
         ylab=ylab,
         xlim = xlim,  
         ylim = ylim,
         xaxt="n", 
         yaxt="n",
         main=main,
         cex=cex,
         cex.main = cex.lab,
         cex.lab = cex.lab
         )
    ticks <- seq(from=floor(xlim[1]),to=ceiling(xlim[2]),by=1)
    labels <- 10^ticks
    axis(1, at=ticks, labels=labels,cex.axis=cex.axis)
    ticks <- seq(from=floor(ylim[1]),to=ceiling(ylim[2]),by=1)
    labels <- 10^ticks
    axis(2, at=ticks, labels=labels,cex.axis=cex.axis)


    points(log10(h2$mids),
           log10(h2$density),
           pch=16,
           cex=cex,
           col=col2
           )
    legend('topright', 
           legend = c(legend1, legend2), 
           col = c(col1, col2), 
           pch=16,
           cex = cex.lab
           ) 
}

hist2.log <- function(dat1, dat2, nbreaks=100, log.breaks=T, main="", xlab="Value", ylab="Frequency", legend1="Data 1", legend2="Data 2",xlim=NULL,ylim=NULL,cex.lab=1.5,cex.axis=1.4,cex=1.2,col1=rgb(0,0,1,3/4),col2 = rgb(1,0,0,3/4)){
    "Produces log plots of the histograms for two sets of data,
    with alpha color blending."

    dat1 <- dat1[dat1>0]
    dat2 <- dat2[dat2>0]
    if (log.breaks==TRUE){
        breaks <- 10^seq(from=log10(min(dat1,dat2)),
                         to=log10(max(dat1,dat2)),
                         length.out=nbreaks)
    } else {
        breaks <- seq(from=min(dat1,dat2),
                      to=max(dat1,dat2),
                      length.out=nbreaks
                      )
    }
    
    h1 <- hist(dat1,breaks=breaks, plot = FALSE)
    h2 <- hist(dat2,breaks=breaks, plot = FALSE)
    if (!length(xlim)){
        xlim <- log10(c(min(h1$mids,h2$mids), max(h1$mids, h2$mids)))
    }
    if (!length(ylim)){
        ylim <- c(min(h1$density[h1$density>0],h2$density[h2$density>0]), max(h1$density, h2$density))
    }

    plot(log10(h1$mids),
         h1$density,
         pch=16,
         col=col1,
         xlab=xlab,
         ylab=ylab,
         xlim = xlim,  
         ylim = ylim,
         xaxt="n", 
         yaxt="n",
         main=main,
         cex=cex,
         cex.main = cex.lab,
         cex.lab = cex.lab
         )
    ticks <- seq(from=floor(xlim[1]),to=ceiling(xlim[2]),by=1)
    labels <- 10^ticks
    axis(1, at=ticks, labels=labels,cex.axis=cex.axis)
    ticks <- seq(from=floor(ylim[1]),to=ceiling(ylim[2]),by=1)
    labels <- 10^ticks
    axis(2, at=ticks, labels=labels,cex.axis=cex.axis)


    points(log10(h2$mids),
           h2$density,
           pch=16,
           cex=cex,
           col=col2
           )
    legend('topright', 
           legend = c(legend1, legend2), 
           col = c(col1, col2), 
           pch=16,
           cex = cex.lab
           ) 
}


hist.loglog <- function(dat1, nbreaks=30, main="", xlab="Value", ylab="Frequency", xlim=NULL,ylim=NULL,col1=rgb(0,0,1,3/4),cex.lab=1.5, cex.main=1.5,cex.axis=1.5,mar=c(6,6,5,5),cex=1.2, ...){
    "Produces log plots of the histograms for a set of data,
    with alpha color blending."

    par(mar=mar)
    par(cex.lab=cex.lab)
    par(cex.main=cex.main)
    par(cex.axis=cex.axis)

    dat1 <- dat1[dat1>0]
    breaks <- 10^seq(from=log10(min(dat1)),
                     to=log10(max(dat1)),
                     length.out=nbreaks)
    breaks <- seq(min(dat1),
                  to=max(dat1),
                  length.out=nbreaks)
    
    
    h1 <- hist(dat1,breaks=breaks, plot = FALSE)
    if (!length(xlim)){
        #xlim <- log10(c(min(h1$mids), max(h1$mids)))
        xlim <- c(min(h1$mids), max(h1$mids))
    }
    if (!length(ylim)){
        #ylim <- log10(c(min(h1$density[h1$density>0]), max(h1$density)))
        ylim <- c(min(h1$density[h1$density>0]), max(h1$density))
    }

    plot(h1$mids,
         h1$density,
         log='xy',
         pch=16,
         col=col1,
         xlab=xlab,
         ylab=ylab,
         xlim = xlim,  
         ylim = ylim,
         #xaxt="n", 
         #yaxt="n",
         main=main,
         cex=cex,
         cex.main = cex.lab,
         cex.lab = cex.lab
         )
    #ticks <- seq(from=floor(xlim[1]),to=ceiling(xlim[2]),by=1)
    #labels <- 10^ticks
    #axis(1, at=ticks, labels=labels,cex.axis=cex.axis)
    #ticks <- seq(from=floor(ylim[1]),to=ceiling(ylim[2]),by=1)
    #labels <- 10^ticks
    #axis(2, at=ticks, labels=labels,cex.axis=cex.axis)

    return(h1)

}

hist.cumul <- function(dat, log = 'xy', main="", xlab="Value", ylab="Cumulative distribution", 
                       xlim=NULL, ylim=NULL, col=rgb(0,0,1,3/4),
                       cex.lab=1.5, cex.main=1.5,cex.axis=1.5,
                       mar=c(6,6,5,5),pch=16,cex=1.2, ...){
    "Produces cumulative distribution of data, with alpha color blending."

    par(mar=mar)
    par(cex.lab=cex.lab)
    par(cex.main=cex.main)
    par(cex.axis=cex.axis)

    x <- na.omit(dat)

    plot(x,
         1 - ecdf(x)(x),
         log=log,
         pch=pch,
         col=col,
         xlab=xlab,
         ylab=ylab,
         xlim = xlim,  
         ylim = ylim,
         main=main,
         cex=cex,
         cex.main = cex.lab,
         cex.lab = cex.lab,
         ...
         )
}

hist.0 <- function(x, breaks=20, col='gray', border='NA', main='', cex.lab=1.5, 
                   cex.main=1.5,cex.axis=1.5,mar=c(6,6,5,5),...){
    par(mar=mar)
    par(cex.lab=cex.lab)
    par(cex.main=cex.main)
    par(cex.axis=cex.axis)
    hist(x,
         breaks=breaks,
         col=col,
         border=border,
         main=main,
         ...
         )
}

hist.arrow <- function(x, x_arrows, col_arrows=rep('red',length(x_arrows)),fac=4, 
                       xlim = c(0.9 * min(x, x_arrows), 1.1 * max(x, x_arrows)),
                       lwd=6, cex.arrow=1.2,
                       breaks=20, col='gray', border='NA', main='', cex.lab=1.5, 
                       cex.main=1.5,cex.axis=1.5,mar=c(6,6,5,5),...){
    par(mar=mar)
    par(cex.lab=cex.lab)
    par(cex.main=cex.main)
    par(cex.axis=cex.axis)
    h <- hist(x,
         breaks=breaks,
         col=col,
         xlim=xlim,
         border=border,
         main=main,
         ...
         )
    for (j in 1:length(x_arrows)){
        x0 <- x_arrows[j]
        arrows(x0=x0,
               y0=max(h$counts)/fac,
               y1=0,
               lwd=lwd,
               col=col_arrows[j]
               )
        text(x=x0,
             y=max(h$counts)/fac,
             labels=names(x_arrows)[j],
             pos=3,
             cex=cex.arrow,
             col=col_arrows[j]
             )
    }
}


hist.arrow.mat <- function(x, inds_arrow, col_arrows=rep('red',length(inds_arrow)), 
                       breaks=20, col='gray', border='NA', main='', cex.lab=1.5, 
                       cex.main=1.5,cex.axis=1.5,mar=c(6,6,5,5),...){
    par(mar=mar)
    par(cex.lab=cex.lab)
    par(cex.main=cex.main)
    par(cex.axis=cex.axis)
    h <- hist(x,
         breaks=breaks,
         col=col,
         border=border,
         main=main,
         ...
         )
    for (j in 1:length(inds_arrow)){
        ind <- inds_arrow[j]
        arrows(x0=x[ind],
               y0=max(h$counts)/4,
               y1=0,
               lwd=6,
               col=col_arrows[j]
               )
        text(x=x[ind],
             y=max(h$counts)/4,
             labels=names(x)[ind],
             pos=3,
             cex=1.2,
             col=col_arrows[j]
             )
    }
}

heatmap.0 <- function(dat, main="", col=colorRampPalette(c("white","red")), 
                      symbreaks=FALSE, cellnote = FALSE, notecex=1, key = TRUE, 
                      Rowv = NA, Colv = NA, symm = FALSE,breaks = 16, hlim=FALSE,
                      dendrogram='none',keysize=1,
                      na.color='gray',
                      cexRow = 1.3, cexCol = 1.3, cex = NULL, margin = c(10,10), ...){

    require(gplots)

    if (length(cex)){
        cexRow <- cex
        cexCol <- cex
    }
    if (length(cellnote) == 1){
        if (as.numeric(cellnote) != 0) cellnote <- signif(dat,2)
        else cellnote <- matrix("", ncol = ncol(dat), nrow = nrow(dat))
    } else if (length(cellnote) > 1){
        cellnote <- signif(cellnote, 2)
    }
    
    heatmap.2(
              dat,
              trace = "none",
              cellnote = cellnote,
              notecex = notecex,
              na.color=na.color,
              notecol = "black",
              col = col,
              main = main,
              symbreaks = symbreaks,
              key = key,
              density.info = "none",
              keysize = keysize, 
              dendrogram = dendrogram,
              margin = margin,
              symm = symm,
              Rowv = Rowv,
              Colv = if(symm) Rowv else Colv,
              breaks = breaks,
              cexRow = cexRow,
              cexCol = cexCol,
              ...
              #revC=TRUE
              )
}

heatmap.sym <- function(dat,main="",col=colorRampPalette(c("blue","white","red")),
                        cellnote = FALSE, notecex=.8, key = TRUE, Rowv = T, Colv = T,
                        breaks = 16, cex = 1, symm=TRUE, dendrogram='both',...){

    heatmap.0(dat,
              main=main,
              col=col,
              symbreaks=TRUE,
              symm=symm,
              cellnote=cellnote,
              notecex=notecex,
              key=key,
              Rowv=Rowv,
              Colv=Colv,
              breaks=breaks,
              dendrogram=dendrogram,
              cexRow=cex,
              cexCol=cex,
              ...
              )

}

heatmap.cor <- function(mat, method='p', by=NA, Nmulti=NA, ...){
    # mat has N columns of variables to correlate along M columns

    M <- nrow(mat)
    
    if (is.na(Nmulti)){
      N <- ncol(mat)
    } else {
      N <- Nmulti
    }

    crs <- sapply(1:1000, function(i) cor.0(runif(M),runif(M)))
    x=seq(0,1,0.001);
    cutoff <- x[which(1-ecdf(abs(crs))(x) < 0.05 / (N*(N-1)/2))[1]]
    if (is.na(by)){
      bru <- seq(cutoff,1,cutoff)
    } else {
      bru <- seq(cutoff,1,by)
    }
    breaks <- c(sort(-bru),bru)

    cc <- cor.0(mat, method=method)
    cc[is.na(cc)] <- 0

    heatmap.sym(cc,
                breaks=breaks,
                main=paste('Bonferroni cutoff =',cutoff),
                key.title='',
                key.xlab='Correlation',
                ...
                )

}

whitered <- function(...){
    return(colorRampPalette(c("white","red"))(...))
}

boxplot.0 <- function(data, main="", ylab="", cex.lab=1.5, cex.main=1.5,cex.axis=1.5,cols='lightgray',cex=0.5,mar=c(6,6,5,5),...){
    par(mar=mar)
    par(cex.lab=cex.lab)
    par(cex.main=cex.main)
    par(cex.axis=cex.axis)
    boxplot(data,
            main = main,
            cex = cex,
            col = cols,
            ylab = ylab,
            pch = 19,
            ...
            ) 
}

boxplot.1 <- function(data, main="", ylab="", cex.lab=1.5, cex.main=1.5,cex.axis=1.5,cols='lightgray',cex=0.5,mar=c(6,6,5,5),...){
    par(mar=mar)
    par(cex.lab=cex.lab)
    par(cex.main=cex.main)
    par(cex.axis=cex.axis)
    boxplot(data,
            main = main,
            cex = 0,
            col = cols,
            ylab = ylab,
            ...
            ) 
    stripchart(data, cex=cex,
               vertical = TRUE, method = "jitter", 
               pch = 21, col = "maroon", bg = "bisque", 
               add = TRUE) 
}



paste.csv <- function(...){
    paste(..., sep=',')
}

paste.tab <- function(...){
    paste(..., sep='\t')
}

pvalue_hyper <- function(q, m, n, k, lower.tail=FALSE){
    # q size intersection
    # m size of B
    # n size of universe
    # k size of A
    if (lower.tail == FALSE){
        phyper(q, m, n, k, lower.tail = FALSE)+ dhyper(q, m, n, k) 
    } else{
        phyper(q, m, n, k)
    }
}


pvalue_hyper.0 <- function(A, B, num_tot, lower.tail=FALSE){

    setA <- unique(A)
    setB <- unique(B)
    nA <- length(setA)
    nB <- length(setB)
    AB <- intersect(setA,setB)
    nAB <- length(AB)
    
    return(pvalue_hyper(nAB, nB, num_tot - nB, nA, lower.tail = lower.tail))

}

kmeansAIC <- function(fit){
    m <- ncol(fit$centers)
    n <- length(fit$cluster)
    k <- nrow(fit$centers)
    D <- fit$tot.withinss
    return(D + 2*m*k)
}

kmeansBIC <-function(fit){
    m <- ncol(fit$centers)
    n <- length(fit$cluster)
    k <- nrow(fit$centers)
    D <- fit$tot.withinss
    return(D + log(n)*m*k)
}

len <- function(...){
    length(...)
}


jaccard <- function(set_A, set_B, asym=F){
    # asym=T for enrichment in first set
    A <- set_A
    B <- set_B
    if (length(A) == 0 & length(B) == 0){
        return(0)
    } else {
        if (asym==F)
            return(length(intersect(A,B)) / length(union(A,B)))
        else
            return(length(intersect(A,B)) / length(A))
    }
}

jaccards <- function(l, asym=F){
    # compute jaccard indices between elements of a list
    # asym=T for enrichment in smallest set
    N <- length(l)
    matJac <- matrix(1,N,N)
    for (i in 1:(N-1)){
        for (j in (i+1):N){
            matJac[i,j] <- jaccard(l[[i]],l[[j]],asym=asym)
            matJac[j,i] <- matJac[i,j]
        }
    }
    rownames(matJac) <- names(l)
    colnames(matJac) <- names(l)
    return(matJac)
}

sourceMy <- function(){
    source("~/Dropbox/Software/Library/R/my.R")
}

loadBioconductor <- function(){
    source("http://bioconductor.org/biocLite.R")
}

loadRA <- function(n='abcon'){
    if (n == 'abcon'){
        load('~/Dropbox/Work/projects/scipher/classifier/RA/Data-third-party/data_Abcon.RData',envir=.GlobalEnv)
    } else if (n == 'abcon_0.5'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Abcon_0.5.RData',envir=.GlobalEnv)
    } else if (n == 'abcon_0.8'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_abcon_0.8.RData',envir=.GlobalEnv)
    } else if (n == 'batterup_2013'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Batterup_2013.RData',envir=.GlobalEnv)
    } else if (n == 'japanese'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Japanese.RData',envir=.GlobalEnv)
    } else if (n == 'batterup_2014'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Batterup_2014.RData',envir=.GlobalEnv)
    } else if (n == 'batterup_2014_0.5'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Batterup_2014_0.5.RData',envir=.GlobalEnv)
    } else if (n == 'batterup_2014_0.8'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Batterup_2014_0.8.RData',envir=.GlobalEnv)
    } else if (n == 'batterup_COMBAT'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Batterup_COMBAT.RData',envir=.GlobalEnv)
    } else if (n == 'batterup_GENENORM'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Batterup_GENENORM.RData',envir=.GlobalEnv)
    } else if (n == 'batterup_BMC'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/data_Batterup_BMC.RData',envir=.GlobalEnv)
    } else if (n == 'all'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/merging_no_japanese/data_all.Rdata',envir=.GlobalEnv)
    } else if (n == 'all_0.2'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/merging_no_japanese/data_all_0.2.RData',envir=.GlobalEnv)
    } else if (n == 'all_0.3'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/merging_no_japanese/data_all_0.3.RData',envir=.GlobalEnv)
    } else if (n == 'all_0.5'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/merging_no_japanese/data_all_0.5.RData',envir=.GlobalEnv)
    } else if (n == 'all_0.8'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/merging_no_japanese/data_all_0.8.RData',envir=.GlobalEnv)
    } else if (n == 'all_japanese'){
        load('~/Dropbox/Work/projects/scipher/Data-third-party/merging/data_all.Rdata',envir=.GlobalEnv)
    }
}

loadACOS <- function(){
    load("/Users/marc/Desktop/Dutch grant/amitabh/ACOS/raw/data.RData",envir=.GlobalEnv)
}

loadHMDP <- function(){
    load('~/Dropbox/Work/data/HMDP/data.RData',envir=.GlobalEnv)
}

loadHMDP_SNPs <- function(){
    load('~/Dropbox/Work/data/HMDP/SNPs+eQTLs/Genotypes/data.RData',envir=.GlobalEnv)
}

loadHumanHF <- function(){
    load('~/Dropbox/Work/data/HMDP/Human/data.RData',envir=.GlobalEnv)
}

loadVDAART <- function(){
    load('~/Dropbox/Work/projects/amitabh/VDAART/data/VDAART/data.RData',envir=.GlobalEnv)
}

loadPPI <- function(v='2016'){
    print(paste('Loading PPI version',v))

    if (v=='2013'){
        load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/PPI-2013/PPI_July_2013.RData',envir=.GlobalEnv)
    } else if (v == 'direct'){
        load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/PPI-2013/PPI_July_2013_direct.RData',envir=.GlobalEnv)
    } else if (v == '2015'){
        load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/PPI-2015/PPI_2015.RData',envir=.GlobalEnv)
    } else if (v == '2016'){
        load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/PPI-2016/PPI_2016.RData',envir=.GlobalEnv)
    } else if (tolower(v)=='huri'){
        load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/PPI-huri/PPI_huri.RData',envir=.GlobalEnv)
    } else if (tolower(v)=='humannet'){
        load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/Functional_network/interactome_HumanNet.RData',envir=.GlobalEnv)
    } else {
        print(paste('Could not load PPI',v))
    }
}

loadPPI_huri <- function(){
    load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/PPI-huri/PPI_huri.RData',envir=.GlobalEnv)
}

loadPPI_humannet <- function(){
    load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/Functional_network/interactome_HumanNet.RData',envir=.GlobalEnv)
}

loadPPI_direct <- function(){
    load('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/PPI-2013/PPI_July_2013_direct.RData',envir=.GlobalEnv)
}

loadZscores <- function(){
    read.csv('~/Dropbox/Work/data/PPI/PPI_updated_may_2012/LCC_Zscores/Zscores_1000.csv')
}

loadPathways <- function(){
    #lpaths <- readLines('~/Dropbox/Work/data/pathways/c2.cp.v4.0.entrez.gmt.txt')
    #lpaths <- readLines('~/Dropbox/Work/data/pathways/msigdb.v3.1.entrez.with.wikipathways_and_custom.gmt')
    lpaths <- readLines('~/Dropbox/Work/data/pathways/pathways_biocarta_kegg_reactome_wikipathways.txt')

    pathways <- list() # list of pathway genes (entrez)
    for (l in lpaths){
        words <- strsplit(l[1],'\t')[[1]]
        if (length(grep('MODULE_', words[1]))) next
        pathways[[words[1]]] <- words[-c(1,2)]
    }

    return(pathways)
}

loadDiseases <- function(){
    lpaths <- readLines('~/Dropbox/Work/data/disease_genes/HuGE/all/Disease-Gene_for_jorg.tsv')

    diseases <- list() # list of pathway genes (entrez)
    for (l in lpaths){
        words <- strsplit(l[1],'\t')[[1]]
        diseases[[words[1]]] <- words[-c(1,2)]
    }

    return(diseases)
}

loadDiseases_OMIM <- function(){
    lpaths <- readLines('~/Dropbox/Work/data/disease_genes/diseaseConnect/disease_gene_OMIM.tsv')

    diseases <- list() # list of pathway genes (entrez)
    for (l in lpaths){
        words <- strsplit(l[1],'\t')[[1]]
        diseases[[words[1]]] <- words[-c(1,2)]
    }

    return(diseases)
}

loadDiseases_GWAS <- function(){
    lpaths <- readLines('~/Dropbox/Work/data/disease_genes/GWAS_Jorg/GWAS_all.tsv')

    diseases <- list() # list of pathway genes (entrez)
    for (l in lpaths){
        words <- strsplit(l[1],'\t')[[1]]
        diseases[[words[1]]] <- words[-c(1,2)]
    }

    return(diseases)
}


loadDiseaseConnect <- function(){
    lpaths <- readLines('~/Dropbox/Work/data/disease_genes/diseaseConnect/disease_gene_all.tsv')

    diseases <- list() # list of pathway genes (entrez)
    for (l in lpaths){
        words <- strsplit(l[1],'\t')[[1]]
        diseases[[words[1]]] <- words[-c(1,2)]
    }

    return(diseases)
}



loadPathwaysAndDiseases <- function(){
    #lpaths <- readLines('~/Dropbox/Work/data/pathways/c2.cp.v4.0.entrez.gmt.txt')
    lpaths <- readLines('~/Dropbox/Work/data/pathways/msigdb.v3.1.entrez.with.wikipathways.and.custom.and.huge.gmt')

    pathways <- list() # list of pathway genes (entrez)
    for (l in lpaths){
        words <- strsplit(l[1],'\t')[[1]]
        pathways[[words[1]]] <- words[-c(1,2)]
    }

    return(pathways)
}



loadmirDIP <- function(){
    load('~/Dropbox/Work/data/miRNA/mirDIP/data.RData', envir=.GlobalEnv)
}

convert_hg19_to_entrez <- function(names){
    t <- read.delim('~/Dropbox/Work/data/gene_names_converter/entrez_to_hg19.tsv',header=F)
    t[,2] <- as.character(t[,2])
    t[which(is.na(t[,2])),2] <- NaN
    return(t[match(as.character(names),as.character(t[,2])),1])
}

convert_entrez_to_hg19 <- function(names){
    t <- read.delim('~/Dropbox/Work/data/gene_names_converter/entrez_to_hg19.tsv',header=F)
    return(as.character(t[match(names,t[,1]),2]))
}


head0 <- function(A,n=10,m=10){
    A[1:n,1:m]
}


plot.memory <- function(){
    # create function to return matrix of memory consumption
    object.sizes <- function()
    {
        return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name) 
                               object.size(get(object.name))))))
    }

    # print to console in table format
    object.sizes()

    # draw pie chart
    pie(object.sizes(), main="Memory usage by object")
}

save.append <- function(x, file) {
    old.objects <- load(file, new.env())
    save(list = c(old.objects, deparse(substitute(x))), file = file)
}

rm.all <- function(){
    rm(list=ls())
}

 
write.delim <-function(tab, file='', ...){
    write.table(tab, 
                file=file, 
                sep='\t',
                quote=F,
                row.names=F,
                ...)
}

install.packages.url <- function(url, ...){
    install.packages(url, repos=NULL, method="libcurl", ...)
}


#write.csv <-function(tab, file="", ...){
    #write.table(tab, 
                #file=file, 
                #col.names=T,
                #row.names=T,
                #sep=',',
                #quote=F,
                #...)
#}

compEigengene <- function(datExpr){
    # Compute the eigengene given a NxM expression matrix
    # N are individuals
    # M are genes
    # returns weights of the eigengene and its value
    texpr_scale <- t(scale(datExpr))
    cor_expr <- cor(t(texpr_scale), use="pairwise.complete.obs")
    e <- eigen(cor_expr)
    eigen_weights <- e$values[1] * e$vectors[,1]
    mat <- as.matrix(texpr_scale)
    mat[is.na(mat)] <- 0
    eigen_best <- as.numeric(eigen_weights %*% mat)
    mean_expr <- apply(texpr_scale,2,mean,na.rm=T)
    corr_mean_expr <- cor(eigen_best, mean_expr, use="pairwise.complete.obs")
    if (sign(corr_mean_expr) < 0){
        eigen_best <- -eigen_best
        eigen_weights <- -eigen_weights
    }
    return(list(
                weights=eigen_weights,
                eigengene=eigen_best
                )
    )
}

compBottomUp <- function(datExpr, pheno, maxval=2000, nbins=30, method='pearson', inds_ord=NA){
    # 1/ Ranks genes by correlation with pheno
    # 2/ Computes eigengene correlation with pheno as function of rank
    if (length(inds_ord) == 1 & is.na(inds_ord)){
        cc <- cor(datExpr, pheno, use='pairwise.complete.obs', method=method)
        inds_ord <- order(-abs(cc))
    }

    nseq <- unique(floor(exp(seq(0, log(min(maxval, ncol(datExpr))),length.out=nbins))))

    ccs_eig <- NULL
    for (i in nseq){
        print(paste('Number of genes:',i,'/',max(nseq)))
        eig <- compEigengene(datExpr[,inds_ord[1:i]])
        ccs_eig <- c(ccs_eig, abs(cor(eig$eigengene,pheno,use='pairwise.complete.obs')))
    }
    return(list(
                'inds_ord'=inds_ord,
                'nseq'=nseq,
                'ccs_eig'=ccs_eig
                ))

}



compBottomUpCV <- function(datExpr, pheno, PROP_TRAIN=0.7, N_CROSSVAL=10, maxval=2000, nbins=30, method='pearson'){
    # Divide in training and test set with PROP_TRAIN
    # compute eigengene on training set and test correlation on test set
    
    nseq <- unique(floor(exp(seq(0, log(min(maxval, ncol(datExpr))),length.out=nbins))))

    ##################################
    ## Indices for cross-validation ##
    ##################################
    N <- nrow(datExpr)
    # create lists of indices for cross-validation
    N_train <- floor(PROP_TRAIN * N)
    inds_trains <- list()
    inds_tests <- list()
    for (ncross in 1:N_CROSSVAL){
        inds_t <- sample(1:N, N_train)
        inds_trains[[ncross]] <- inds_t
        inds_tests[[ncross]] <- (1:N)[-inds_t]
    }

    #############################
    ## Ordering by correlation ##
    #############################
    # first we compute correlations
    inds_ords <- list()
    for (ncross in 1:N_CROSSVAL){ #crossvalidation
        inds_train <- inds_trains[[ncross]]
        inds_test <- inds_tests[[ncross]]
        cc <- cor(datExpr[inds_train,], pheno[inds_train], 
                  use='pairwise.complete.obs', method=method)
        inds_ords[[ncross]] <- order(-abs(cc))
    }
    
    #######################################
    ## Computing eigengenes correlations ##
    #######################################
    all_cors <- matrix(ncol=length(nseq), nrow=N_CROSSVAL)
    for (i in 1:length(nseq)){
        print(paste('Number of genes:',nseq[i],'/',max(nseq)))
        
        for (ncross in 1:N_CROSSVAL){ #crossvalidation
            #print(paste('Cross-validation step',ncross,'/',N_CROSSVAL))
            
            inds_train <- inds_trains[[ncross]]
            inds_test <- inds_tests[[ncross]]
            inds_ord <- inds_ords[[ncross]]
            inds_genes <- inds_ord[1:nseq[i]]
            
            eig <- compEigengene(datExpr[inds_train,inds_genes])
            eigen_weights <- eig$weights
            mat <- as.matrix(t(datExpr[inds_test,inds_genes]))
            mat[is.na(mat)] <- 0
            eigen_test <- as.numeric(eigen_weights %*% mat)
            all_cors[ncross,i] <- abs(cor(eigen_test, pheno[inds_test], use="pairwise.complete.obs", method=method))
            
        }
        print(paste0('Mean correlation: ', mean(all_cors[,i])))
    }
    

    return(list(
                'inds_trains'=inds_trains,
                'inds_tests'=inds_tests,
                'inds_ords',inds_ords,
                'nseq'=nseq,
                'all_cors'=all_cors
                ))

}


list2matrix <- function(l){
    M <- matrix(unlist(l),ncol=length(l))
    colnames(M) <- names(l)
    return(M)
}

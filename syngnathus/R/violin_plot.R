

#' Create violin plots
#' @author Sarah P. Flanagan
#' @param x Data to plot, either a vector or a dataframe with multiple columns.
#' @param ... Additional datasets (each should be a vector)
#' @param range Range
#' @param h h
#' @param ylim y-axis limits
#' @param names Axis names
#' @param horizontal Plot horizontally or vertically
#' @param col Color for the violin backgrounds
#' @param border Color for the violin borders
#' @param lty Line type
#' @param lwd Line width
#' @param rectCol Rectangle color
#' @param colMed color med
#' @param pchMed shape for med
#' @param at plot param at
#' @param add Add the plot or not
#' @param wex wex
#' @param drawRect Add a rectangle to the figure
#' @param summaryLines Add lines at values given by data.frame with each column corresponds to a factor/data (default is NULL)
#' @param sumLineAdj An adjustment to add a a little bit beyond the violin borders, only when summaryLines is provided. (default is 0)
#' @param plot.axes Whether or not to plot axes
#' @param axis.box Add the axis box
#' @param plot.ann Add the annotations to the plot
#' @return Optionally returns the summary statistics of the data.
#' @note Modified from vioplot() in library(vioplot)
#' @examples
#' mu<-2
#' si<-0.6
#' bimodal<-c(rnorm(1000,-mu,si),rnorm(1000,mu,si))
#' uniform<-runif(2000,-4,4)
#' normal<-rnorm(2000,0,3)
#' violin_plot(bimodal,uniform,normal,col=c("red","blue","green"))
#' @import sm
#' @export
violin_plot <- function(x,...,range=1.5,h=NULL,ylim=NULL,names=NULL, horizontal=FALSE,
                        col="magenta", border="black", lty=1, lwd=1, rectCol="black", colMed="white", pchMed=19, at, add=FALSE, wex=1,
                        drawRect=TRUE,plot.axes=TRUE,axis.box=FALSE,plot.ann=TRUE, summaryLines=NULL,sumLineAdj=0)
{
  require(sm)
  # process multiple datas
  
  if(ncol(x)>1){
    datas<-as.list(x)
    others<-list(...)
    if(length(others)>0){
      nd<-length(datas)
      for(d in 1:length(others)){
        datas[[nd+d]]<-others[[d]]
      }
    }
  }else{
    datas <- list(x,...)
  }
  
  n <- length(datas)
  
  if(missing(at)) at <- 1:n

  # pass 1
  #
  # - calculate base range
  # - estimate density
  #

  # setup parameters for density estimation
  upper  <- vector(mode="numeric",length=n)
  lower  <- vector(mode="numeric",length=n)
  q1     <- vector(mode="numeric",length=n)
  q3     <- vector(mode="numeric",length=n)
  med    <- vector(mode="numeric",length=n)
  base   <- vector(mode="list",length=n)
  height <- vector(mode="list",length=n)
  baserange <- c(Inf,-Inf)

  # global args for sm.density function-call
  args <- list(display="none")

  if (!(is.null(h)))
    args <- c(args, h=h)

  for(i in 1:n) {
    data<-datas[[i]]

    # calculate plot parameters
    #   1- and 3-quantile, median, IQR, upper- and lower-adjacent
    data.min <- min(data)
    data.max <- max(data)
    q1[i]<-quantile(data,0.25)
    q3[i]<-quantile(data,0.75)
    med[i]<-median(data)
    iqd <- q3[i]-q1[i]
    upper[i] <- min( q3[i] + range*iqd, data.max )
    lower[i] <- max( q1[i] - range*iqd, data.min )

    #   strategy:
    #       xmin = min(lower, data.min))
    #       ymax = max(upper, data.max))
    #

    est.xlim <- c( min(lower[i], data.min), max(upper[i], data.max) )

    # estimate density curve
    smout <- do.call("sm.density", c( list(data, xlim=est.xlim), args ) )

    # calculate stretch factor
    #
    #  the plots density heights is defined in range 0.0 ... 0.5
    #  we scale maximum estimated point to 0.4 per data
    #
    hscale <- 0.4/max(smout$estimate) * wex

    # add density curve x,y pair to lists
    base[[i]]   <- smout$eval.points
    height[[i]] <- smout$estimate * hscale

    # calculate min,max base ranges
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1],t[1])
    baserange[2] <- max(baserange[2],t[2])

  }

  # pass 2
  #
  # - plot graphics

  # setup parameters for plot
  if(!add){
    xlim <- if(n==1)
      at + c(-.5, .5)
    else
      range(at) + min(diff(at))/2 * c(-1,1)

    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  } else {
    label <- names
  }
  # setup colors and borders
  if(length(col) < n){
    col<-c(col,rep(col,n-length(col)))
  }
  if(length(border)<n){
    border<-c(border,rep(border,n-length(border)))
  }
  if(length(colMed)<n){
    colMed<-c(colMed,rep(colMed,n-length(colMed)))
  }
  boxwidth <- 0.05 * wex

  # setup plot
  if(!add)
    plot.new()
  if(!horizontal) {
    if(!add){
      plot.window(xlim = xlim, ylim = ylim)
      if(plot.axes){
        if(plot.ann){
          axis(2)
          axis(1,at = at, label=label )
        }else{
          axis(2,labels=F)
          axis(1,at = at, labels=F )
        }
      }
    }

    if(axis.box){ box() }
    for(i in 1:n) {
      # plot left/right density curve
      polygon( c(at[i]-height[[i]], rev(at[i]+height[[i]])),
               c(base[[i]], rev(base[[i]])),
               col = col[i], border=border[i], lty=lty, lwd=lwd)

      if(drawRect){
        # plot IQR
        lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty)

        # plot 50% KI box
        rect( at[i]-boxwidth/2, q1[i], at[i]+boxwidth/2, q3[i], col=rectCol)

        # plot median point
        points( at[i], med[i], pch=pchMed, col=colMed )
      }

      if(!is.null(summaryLines)){

        for(sl in summaryLines[,i]){
          diffs<- abs(base[[i]]-sl)
          ht<-height[[i]][which.min(diffs)]+sumLineAdj
          lines( c(at[i]-ht,at[i]+ht),
                 c(sl, sl) ,lwd=lwd, lty=lty,col=border[i])
        }
      }
    }

  }
  else {
    if(!add){
      plot.window(xlim = ylim, ylim = xlim,bty=axis.bty)
      if(plot.axes){
        if(plot.ann){
          axis(2)
          axis(1,at = at, label=label )
        }else{
          axis(2,labels=F)
          axis(1,at = at, labels=F )
        }
      }
    }

    if(axis.box){ box() }
    for(i in 1:n) {
      # plot left/right density curve
      polygon( c(base[[i]], rev(base[[i]])),
               c(at[i]-height[[i]], rev(at[i]+height[[i]])),
               col = col[i], border=border[i], lty=lty, lwd=lwd)

      if(drawRect){
        # plot IQR
        lines( c(lower[i], upper[i]), at[c(i,i)] ,lwd=lwd, lty=lty)

        # plot 50% KI box
        rect( q1[i], at[i]-boxwidth/2, q3[i], at[i]+boxwidth/2,  col=rectCol)

        # plot median point
        points( med[i], at[i], pch=pchMed, col=colMed )
      }

      if(!is.null(summaryLines)){

        for(sl in summaryLines[,i]){
          diffs<- abs(base[[i]]-sl)
          ht<-height[[i]][which.min(diffs)]+sumLineAdj
          lines( c(at[i]-ht,at[i]+ht),
                 c(sl, sl) ,lwd=lwd, lty=lty,col=border[i])
        }
      }

    }

  }
  invisible (list( upper=upper, lower=lower, median=med, q1=q1, q3=q3))
}

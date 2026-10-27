# function to plot two variables on different axes
myplot <- function(x, y1, y2, xlab, ylab1, ylab2, add.line, name.out, path.out){
  
  # set png device
  png(paste0(path.out, name.out), width=1600, height=900) 
  
  # set larger right margin
  par(oma = c(0, 0, 0, 2))
  
  # plot y1 on x
  plot(x, y1, xlab = xlab, ylab = ylab1, type = "p", pch = 20, col = "gray", axes = FALSE)
  
  # draw axes
  axis(1)
  axis(2)
  
  # plot y2 on x
  par(new = TRUE)
  plot(x, y2, xlab = "", ylab = "", type = "l", col = "cyan", axes = FALSE)
  
  # add 1% point of N(0,1) if necessary
  if(add.line){
    lines(x, rep(qnorm(0.01), length(x)), col = "orange", lty = 3)    
    lines(x, rep(qnorm(0.0001), length(x)), col = "red", lty = 3)    
  }
  
  
  # draw additional axis
  mtext(ylab2, side = 4, line = 3)
  axis(4)
  
  
  
  # add box
  box()
  
  # set graphic parameters
  my.legend <- c(ylab1, ylab2)
  my.col <- c("gray", "cyan")
  my.lty <- c(0, 1)
  my.pch <- c(20, NA_integer_)
  my.ncol <- 2
  if (add.line) {
    # my.legend <- c(my.legend, "1%", "0.01%")
    my.legend <- c(my.legend, NA, NA)
    my.col <- c(my.col, "orange", "red")
    my.lty <- c(my.lty, 3, 3)
    my.pch <- c(my.pch, rep(NA_integer_, 2))
    my.ncol <- my.ncol + 2
  }
  
  # add legend outside the box
  par(xpd = TRUE)
  cor.x <- par()$usr[1]
  cor.y <- par()$usr[4] + 0.1 * abs(par()$usr[4] - par()$usr[3])
  legend(cor.x, cor.y, legend = my.legend, col = my.col, lty = my.lty, pch = my.pch, ncol = my.ncol)
  
  # device off
  dev.off()
}

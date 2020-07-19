
EPSR_figure <- function(epsr, MCMAX){

		mycolsi <- rainbow(ncol(epsr), s = 1, v = 1, start = 0,
				end = max(1, ncol(epsr) - 1)/ncol(epsr), alpha = 1)

		repxlim <- c(1:(MCMAX-1))
		png('EPSR.png')
		plot(x = repxlim , y = epsr[,1], type="l",
			xlab = "iterations", ylab = "EPSR", ylim = c(0,3), col = mycolsi[1])
		for (i in 2:ncol(epsr))
		{
			lines(x = repxlim, y = epsr[,i], col = mycolsi[i])
		}
		dev.off()
		#try_plot = try(savePlot(filename = "EPSR",
		#		type = "png",
		#		device = dev.cur(),
		#		restoreConsole = TRUE))
		#if("try-error" %in% class(try_plot))
		#{
		#  cat('If you are working with RStudio, the plot can be exported from menu in plot panel (lower right-pannel).   \n')
		#}
}
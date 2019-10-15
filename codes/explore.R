library(ggplot2)
library(reshape2)


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


data = read.table("chr22_60smp_train.txt",header=T)

size = nrow(data)


#### H3K4ME3

data_H3K4ME3 = rbind(cbind(rep("H3K4ME3_1",size),data[,1],data[,10]),cbind(rep("H3K4ME3_2",size),data[,2],data[,10]),cbind(rep("H3K4ME3_3",size),data[,3],data[,10]))
data_H3K4ME3 = data.frame("group" = factor(data_H3K4ME3[,1]),"read" = as.numeric(data_H3K4ME3[,2]),"HiC" = factor(data_H3K4ME3[,3]))

data_H3K4ME3$read = log2(data_H3K4ME3$read)

g1 = ggplot(data_H3K4ME3,aes(x=read)) + 
       geom_density(aes(col=HiC,linetype=group))+ 
	   xlab("log2(Read)") + 
	   ggtitle("Fig 1a. Distribution of H3K4ME3 reads") + 
	   theme(plot.title=element_text(size=10))


##### H3K4ME1

data_H3K4ME1 = rbind(cbind(rep("H3K4ME1_1",size),data[,4],data[,10]),cbind(rep("H3K4ME1_2",size),data[,5],data[,10]),cbind(rep("H3K4ME1_3",size),data[,6],data[,10]))
data_H3K4ME1 = data.frame("group" = factor(data_H3K4ME1[,1]),"read" = as.numeric(data_H3K4ME1[,2]),"HiC" = factor(data_H3K4ME1[,3]))

data_H3K4ME1$read = log2(data_H3K4ME1$read)

g2 = ggplot(data_H3K4ME1,aes(x=read)) + geom_density(aes(col=HiC,linetype=group))+
  xlab("log2(Read)") + ggtitle("Fig 1b. Distribution of H3K4ME1 reads") +
  theme(plot.title=element_text(size=10))


#### H3K27AC
data_H3K27AC = rbind(cbind(rep("H3K27AC_1",size),data[,7],data[,10]),cbind(rep("H3K27AC_2",size),data[,8],data[,10]),cbind(rep("H3K27AC_3",size),data[,9],data[,10]))
data_H3K27AC = data.frame("group" = factor(data_H3K27AC[,1]),"read" = as.numeric(data_H3K27AC[,2]),"HiC" = factor(data_H3K27AC[,3]))


data_H3K27AC$read = log2(data_H3K27AC$read)

g3 = ggplot(data_H3K27AC,aes(x=read)) + geom_density(aes(col=HiC,linetype=group))+
  xlab("log2(Read)") + ggtitle("Fig 1c. Distribution of H3K27AC reads")+
  theme(plot.title=element_text(size=10))






##### correlation plot

data_cor = melt(cor(data[,c(1:9)],use="p"))

pdf("Correlation between nine bins.pdf")

g4 = ggplot(data_cor,aes(x=Var1,y=Var2)) + 
   geom_tile(aes(fill=value)) +
   scale_fill_gradientn(colours=c("green", "red"),limits=c(0,1))+
   labs(title="Fig 1d. Correlation between features") + 
   xlab("") + 
   ylab("") +
   theme(axis.text.x=element_text(angle=45, hjust=1, size=7),axis.text.y=element_text(vjust=0, size=7),plot.title=element_text(size=10))



### multi plot
pdf("Fig1.pdf")
multiplot(g1,g2,g3,g4,cols=2)
dev.off()



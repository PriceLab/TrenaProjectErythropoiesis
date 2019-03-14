load("KLF1.FLI1.atac.srm.RData")
tbl <- as.data.frame(t(mtx))
tbl$time <- seq_len(nrow(tbl))
tbl$tp <- as.factor(rownames(tbl))
colnames(tbl) <- c("KLF1.srm", "FLI1.srm", "KLF1.total.hits", "KLF1.top.hits", "KLF1.total.density", "KLF1.top.density",
                   "FLI1.total.hits", "FLI1.top.hits", "FLI1.total.density", "FLI1.top.density")
ggplot(data=tbl, mapping = aes(x=KLF1.srm, y = KLF1.top.density)) + geom_point()
ggplot(data=tbl, mapping = aes(x = rownames(tbl), y=KLF1.srm)) + geom_line()
ggplot(data=tbl) + geom_smooth(mapping = aes(x = rownames(tbl), y=KLF1.srm))

ggplot(data=mtcars) + geom_smooth(mapping = aes(x=disp, y=mpg))
ggplot(data=tbl) + geom_point(mapping = aes(x=time, y=KLF1.srm)) + geom_smooth(mapping = aes(x=time, y=KLF1.srm))



x.axis.labels <- c("4r1", "4r2", "8", "10r1", "10r2", "11r1", "11r2", "12r1", "12r2")
m <- aes(x=time, y=KLF1.srm)
m <- aes(x=time, y=KLF1.total.density)
m <- aes(x=time, y=FLI1.srm)
m <- aes(x=time, y=FLI1.total.density)
m <- aes(x=time, y=FLI1.total.hits)
p <- ggplot(data=tbl) + geom_point(mapping=m) + geom_smooth(mapping=m) + xlab("Day & replicate") +
      scale_x_continuous(breaks=1:9, labels=x.axis.labels)
p

m <- aes(x=FLI1.total.density, y=FLI1.srm)
p <- ggplot(data=tbl) + geom_point(mapping=m) + geom_smooth(mapping=m) +
#      scale_x_continuous(breaks=1:9, labels=x.axis.labels)
p

library(grid)
grob1 <- grobTree(textGrob(paste("Pearson Correlation : ", round(cor(tbl$FLI1.total.density, tbl$FLI1.srm), 4) ),
                           x = 0.63, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))

ggplot(tbl, aes(x=FLI1.total.density, y=FLI1.srm)) + geom_point() +
   ggtitle("FLI1.total.density vs FLI1.srm") + geom_smooth(method=lm, se=FALSE) +
   #scale_x_continuous(name = "FLI1.total.density", limits = c(5, 15), breaks = seq(5, 15, 2)) +
   #scale_y_continuous(name = "FLI1.srm", limits = c(4000,8000), breaks = seq(5, 15, 2)) +
   annotation_custom(grob1) +
   theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(color="black"), axis.line.x = element_line(color="black")) +
   theme_bw()






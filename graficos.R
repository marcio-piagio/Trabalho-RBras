save(resultados.2.2, file = "resultados.2.2.RData")


14535.41/3600
# Parar os Clusters

parallel::stopCluster(parallel::makeCluster(4))
future::plan("sequential")

################################################################################

# Gráficos 
resultados.2.2
resultados <- resultados.2.2

grafic.result <- function(k){
  metricas <- c("AIC","BIC","AICc","HQIC","TRV-95")
  library(RColorBrewer)
  cor <- brewer.pal(7, "Dark2")
  
  par(mar=c(2, .1, 0.1, 0.1))
  layout(matrix(c(rep(1, 9),1,2,3,6,7,10,11,14,1,1,4,5,8,9,12,13,14,1,rep(1, 9)),9,4),
         heights=c(0.1,3,1,3,1,3,1,0.6), widths=c(0.1, 10, 10, 0.1))
  
  plot.new()
  box(col="grey40")
  
  par(mar= c(.1, 5, 2, 1))
  
  for (i in 1:6){
    
    
    plot(c.true, resultados[[i]]$AIC[,k], ylab = expression(1-beta), xlab = "c",
         main = paste0('K = ', K.true[i]),
         axes=F, type="b", pch=1, col=cor[1], lty=1, ylim=c(-0.05,1))
    points(c.true, resultados[[i]]$BIC[,k],    type="b", col=cor[2], pch=2, lty=2)
    points(c.true, resultados[[i]]$AICc[,k],   type="b", col=cor[3], pch=3, lty=3)
    points(c.true, resultados[[i]]$HQIC[,k], type="b", col=cor[4], pch=4, lty=4)
    points(c.true, resultados[[i]]$TRV.95[,k], type="b", col=cor[5], pch=5, lty=5)
    # points(c.true, resultados[[k]]$TRV.99[,i], type="b", col=cor[6], pch=6, lty=6)
    axis(1, at= seq(-2,2,0.4), tck = 0.015, mgp = c(1, 0.5, 0),
         col="grey40", col.axis="grey20")
    axis(2,at= seq(0,1,0.1), las=1, tck = 0.015, mgp=c(1, 0.5, 0),
         col="grey40", col.axis="grey20")
    abline(h=0.05, col="red", lty=2)
    box(bty = "l", col="grey40")
    plot.new()
    text(0.5, 0.5, "c", cex=1, font=2, col="grey20")
    
  } 
  # 
  # par(mar=c(.5,18,.5,18))
  # plot.new()
  # text(0.5, 0.5, paste0('K = ', K.true[4]), cex=1.5, font=2, col="grey20")
  # 
  # plot(c.true, resultados[[4]]$AIC[,k], ylab = "", xlab = "c",
  #      axes=F, type="b", pch=1, col=cor[1], lty=1, ylim=c(-0.05,1))
  # points(c.true, resultados[[4]]$BIC[,k],    type="b", col=cor[2], pch=2, lty=2)
  # points(c.true, resultados[[4]]$AICc[,k],   type="b", col=cor[3], pch=3, lty=3)
  # points(c.true, resultados[[4]]$HQIC[,k], type="b", col=cor[4], pch=4, lty=4)
  # points(c.true, resultados[[4]]$TRV.95[,k], type="b", col=cor[5], pch=5, lty=5)
  # # points(c.true, resultados[[k]]$TRV.99[,3], type="b", col=cor[6], pch=6, lty=6)
  # axis(1, seq(-2,2,0.4), tck = 0.015, mgp = c(1.5, 0.5, 0),
  #      col="grey40", col.axis="grey20")
  # axis(2,at= seq(0,1,0.1), las=1, tck = 0.015, mgp=c(1.5, 0.5, 0),
  #      col="grey40", col.axis="grey20")
  # abline(h=0.05, col="red", lty=2)
  # box(bty = "l", col="grey40")
  
  par(mar=c(0.01,2,.01,.01))
  plot(1, type="n", axes=F)
  legend(x="bottom", legend=metricas, col=cor, box.col="grey40",
         lwd=1, lty=1:7,pch=1:7, horiz=T, cex=1, text.col= cor)
}

grafic.result(1)


################################################################################

grafic.result.2 <- function(N,k1,k2,k3,
                            ymim1,ymax1,yby1,
                            ymim2,ymax2,yby2,
                            ymim3,ymax3,yby3){
  metricas <- c("k = 5","k = 10","k = 15")
  library(RColorBrewer)
  cor <- brewer.pal(7, "Dark2")
  
  par(mar=c(0.5,0.1, 0.5, 1.5))
  layout(matrix(c(rep(1, 7),1,2,3,6,7,8,1,1,4,5,6,7,8,1,rep(1, 7)),7,4),
         heights=c(1, 2, 9, 2, 9, 3), widths=c(0.2, 10, 10, 0.2))
  
  plot.new()
  box(col="grey40")
  
  par(mar= c(1.5, 3, .3, 2))
  
  for (i in c(k1,k2)){
    plot.new()
    text(0.5, 0.5, paste0('N = ', N.true[i]), cex=1.5, font=2, col="grey20")
    
    plot(c.true, resultados[[i]]$Vies.Mtb.N[,N], ylab = "", xlab = "c",
         axes=F, type="b", pch=1, col=cor[1], lty=1, ylim=c(ifelse(i==1,ymim1,ymim2),
                                                            ifelse(i==1,ymax1,ymax2)))
    points(c.true, resultados[[2]][[k]][,i],    type="b", col=cor[2], pch=2, lty=2)
    points(c.true, resultados[[3]][[k]][,i],   type="b", col=cor[3], pch=3, lty=3)
    axis(1, at= seq(-2,2,0.4), tck = 0.015, mgp = c(1.5, 0.5, 0),
         col="grey40", col.axis="grey20")
    axis(2,at= seq(ifelse(i==1,ymim1,ymim2),
                   ifelse(i==1,ymax1,ymax2),
                   ifelse(i==1,yby1,yby2)),
         las=1, tck = 0.015, mgp=c(1.5, 0.5, 0),
         col="grey40", col.axis="grey20")
    abline(h=0.05, col="red", lty=2)
    box(bty = "l", col="grey40")
    
  } 
  
  par(mar=c(.5,18,.5,18))
  plot.new()
  text(0.5, 0.5, paste0('N = ', N.true[3]), cex=1.5, font=2, col="grey20")
  
  plot(c.true, resultados[[1]][[k]][,3], ylab = "", xlab = "c",
       axes=F, type="b", pch=1, col=cor[1], lty=1, ylim=c(ymim3,ymax3))
  points(c.true, resultados[[2]][[k]][,3],    type="b", col=cor[2], pch=2, lty=2)
  points(c.true, resultados[[3]][[k]][,3],   type="b", col=cor[3], pch=3, lty=3)
  axis(1, seq(-2,2,0.4), tck = 0.015, mgp = c(1.5, 0.5, 0),
       col="grey40", col.axis="grey20")
  axis(2,at= seq(ymim3,ymax3,yby3), las=1, tck = 0.015, mgp=c(1.5, 0.5, 0),
       col="grey40", col.axis="grey20")
  abline(h=0.05, col="red", lty=2)
  box(bty = "l", col="grey40")
  
  par(mar=c(0.2,.5,.5,.5))
  plot(1, type="n", axes=F)
  legend(x="bottom", legend=metricas, col=cor, box.col="grey40",
         lwd=1, lty=1:7,pch=1:7, horiz=T, cex=1, text.col= cor)
}

# Vies.Mt.N
grafic.result.2(1,1,2,50,-150,125,20,-30,240,15)
# Vies.Mtb.N
grafic.result.2(9,-20,5,2,-5,400,20,-10,600,20)
# Vies.Mtb.c
grafic.result.2(10,-1,0.1,0.1,-0.1,1,0.1,-0.2,1,0.1)
# EQM.Mt.N
grafic.result.2(11,0,320000,20000,0,1600000,100000,0,4800000,300000)
# EQM.Mtb.N
grafic.result.2(12,0,18000,2000,0,400000,20000,0,800000,50000)
# EQM.Mtb.c
grafic.result.2(13,0,1,0.1,0,1.4,0.15,0,1,0.1)


################################################################################


# create a data frame
C <- factor(round(as.numeric(rep(c.true, each = no.simulation*2)), 2))
M <- rep(c("Mt","Mtb"),41,each = no.simulation)

k5.N50 <- data.frame("N" = matrix.N[[1]][,1],"Modelo"=M,"c"=C)
k6.N50 <- data.frame("N" = matrix.N[[2]][,1],"Modelo"=M,"c"=C)
k7.N50 <- data.frame("N" = matrix.N[[3]][,1],"Modelo"=M,"c"=C)
k8.N50 <- data.frame("N" = matrix.N[[4]][,1],"Modelo"=M,"c"=C)
k9.N50 <- data.frame("N" = matrix.N[[5]][,1],"Modelo"=M,"c"=C)
k10.N50 <- data.frame("N" = matrix.N[[6]][,1],"Modelo"=M,"c"=C)
k11.N50 <- data.frame("N" = matrix.N[[7]][,1],"Modelo"=M,"c"=C)
k12.N50 <- data.frame("N" = matrix.N[[8]][,1],"Modelo"=M,"c"=C)
k13.N50 <- data.frame("N" = matrix.N[[9]][,1],"Modelo"=M,"c"=C)
k14.N50 <- data.frame("N" = matrix.N[[10]][,1],"Modelo"=M,"c"=C)
k15.N50 <- data.frame("N" = matrix.N[[11]][,1],"Modelo"=M,"c"=C)

library(ggplot2)
library(RColorBrewer)
cor <- brewer.pal(7, "Dark2")
# grouped boxplot
ggplot(k11.N50, aes(x=c, y=N, fill=Modelo)) + 
  geom_boxplot(show.legend=F)+
  geom_hline(yintercept = 50, color="red", linetype = "dashed")+
  # facet_wrap(~Modelo)+
  ggtitle("N = 50")+
  ylab("")+xlab("")+
  # theme_black()+
  theme(legend.position = "bottom", legend.background = element_rect(colour = "gray50", size = .6))+
  scale_fill_manual("", values = c("Mt" = cor[1], "Mtb" = cor[4]))+
  ylim(-0, 200)


k5.N500 <- data.frame("N" = matrix.N[[1]][,10],"Modelo"=M,"c"=C)
k6.N500 <- data.frame("N" = matrix.N[[2]][,10],"Modelo"=M,"c"=C)
k7.N500 <- data.frame("N" = matrix.N[[3]][,10],"Modelo"=M,"c"=C)
k8.N500 <- data.frame("N" = matrix.N[[4]][,10],"Modelo"=M,"c"=C)
k9.N500 <- data.frame("N" = matrix.N[[5]][,10],"Modelo"=M,"c"=C)
k10.N500 <- data.frame("N" = matrix.N[[6]][,10],"Modelo"=M,"c"=C)
k11.N500 <- data.frame("N" = matrix.N[[7]][,10],"Modelo"=M,"c"=C)
k12.N500 <- data.frame("N" = matrix.N[[8]][,10],"Modelo"=M,"c"=C)
k13.N500 <- data.frame("N" = matrix.N[[9]][,10],"Modelo"=M,"c"=C)
k14.N500 <- data.frame("N" = matrix.N[[10]][,10],"Modelo"=M,"c"=C)
k15.N500 <- data.frame("N" = matrix.N[[11]][,10],"Modelo"=M,"c"=C)

library(ggplot2)
library(RColorBrewer)
cor <- brewer.pal(7, "Dark2")
# grouped boxplot
ggplot(k11.N500, aes(x=c, y=N, fill=Modelo)) + 
  geom_boxplot(show.legend=F)+
  geom_hline(yintercept = 500, color="red", linetype = "dashed")+
  # facet_wrap(~Modelo)+
  ggtitle("N = 500")+
  ylab("")+xlab("")+
  # theme_black()+
  theme(legend.position = "bottom", legend.background = element_rect(colour = "gray50", size = .6))+
  scale_fill_manual("", values = c("Mt" = cor[1], "Mtb" = cor[4]))
  # ylim(-0, 200)


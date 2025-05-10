grafic.result <- function(k, V){
  metricas <- c("AIC","BIC","AICc","HQIC", "CAIC","TRV-90","TRV-95")
  library(RColorBrewer)
  cor <- brewer.pal(7, "Dark2")
  
  par(mar=c(2, .1, 0.1, 0.1))
  layout(matrix(c(rep(1, 9),1,2,3,6,7,10,11,14,1,1,4,5,8,9,12,13,14,1,rep(1, 9)),9,4),
         heights=c(0.1,3,1,3,1,3,1,0.6), widths=c(0.1, 10, 10, 0.1))
  
  plot.new()
  box(col="grey40")
  
  par(mar= c(.1, 5, 2, 1))
  
  for (i in V){
    
    
    plot(c.true, AIC[[i]][,k], ylab = expression(1-beta), xlab = "c",
         main = paste0('K = ', K.true[i]),
         axes=F, type="b", pch=1, col=cor[1], lty=1, ylim=c(-0.05,1))
    points(c.true, BIC[[i]][,k],    type="b", col=cor[2], pch=2, lty=2)
    points(c.true, AICc[[i]][,k],   type="b", col=cor[3], pch=3, lty=3)
    points(c.true,  HQIC[[i]][,k], type="b", col=cor[4], pch=4, lty=4)
    points(c.true,  CAIC[[i]][,k], type="b", col=cor[5], pch=5, lty=5)
    points(c.true, TRV.90[[i]][,k], type="b", col=cor[6], pch=6, lty=6)
    points(c.true, TRV.95[[i]][,k], type="b", col=cor[7], pch=7, lty=7)
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

grafic.result(10,c(1,3,5,7,9,11))


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


tapply(k5.N50$N,k5.N50$Modelo,summary)
tapply(k6.N50$N,k6.N50$Modelo,summary)


library(ggplot2)
library(RColorBrewer)
cor <- brewer.pal(8, "Dark2")

library(ggplot2)
library(ggtext) 

ggplot(k6.N50, aes(x = c, y = N, fill = Modelo)) + 
  geom_boxplot(show.legend = TRUE) +
  geom_hline(yintercept = 50, color = "red", linetype = "dashed") +
  ggtitle("N = 50") +
  ylab("Estimativa do tamanho populacional") +
  xlab("Estimativa do parâmetro de recaptura") +
  # theme_black() +
  theme(
    legend.position = "bottom", 
    legend.background = element_rect(colour = "gray50", size = 0.6),
    axis.text.y = ggtext::element_markdown() # Permite HTML nos rótulos
  ) +
  scale_fill_manual("", values = c("Mt" = cor[1], "Mtb" = cor[8])) +
  scale_y_continuous(
    breaks = seq(-50, max(k6.N50$N, na.rm = T) + 100, by = 50),
    labels = function(x) {
      ifelse(x == 50, "<span style='color:red;'>50</span>", x)
    }
  )


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


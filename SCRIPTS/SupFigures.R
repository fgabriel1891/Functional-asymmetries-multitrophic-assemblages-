##### Appendix figures 

# Figure S1: Diet PCA


biplot(dietRDA,
       scaling = 2, 
       display = "species" , 
       type = "text", 
       frame = F)

summary(dietRDA)

# Figure S2: Climate PCA

png("FigsLast/ClimaticGradient.png", 
    500, 500, pointsize = 20)
biplot(rdaTEMP, 
       type = "text",
       display = "species",
       frame = F)
dev.off()


#### Figure S3: Correlations found with RLQ

par(oma = c(2,2,2,2))

corrplot::corrplot(as.matrix(RLQ$tab),
                   method = "square",
                   is.corr = F, 
                   mar = c(0,0,0,0))

##### Figure S4: Loadings of traits in trait matching axes 


par(mfrow = c(1,2))
corrplot::corrplot(as.matrix(RLQ$c1), 
                   method = "square",
                   is.corr = F,
                   cl.pos	 = "b",
                   cl.cex = 0.4,
                   mar = c(0,0,0,0))


corrplot::corrplot(as.matrix(RLQ$l1), 
                   method = "square",
                   is.corr = F, 
                   cl.pos	 = "b",
                   cl.cex = 0.4,
                   mar = c(0,0,0,0))

##### Figure S5: Trait matching space 

png("FigsLast/FigureS4_RLQ.png", 
    500, 500, pointsize = 20)
plot(RLQ$mQ[,1:2], 
     frame = F, 
     pch = 21, 
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "Trait-matching axis 1",
     ylab = "Trait-matching axis 2",
     col = "orange", 
     bg = "orange")

points(RLQ$mR[,1:2],
       col = "purple")
legend("topleft", 
       c("palms", "mammals"), 
       bty = "n", 
       cex = 0.7,
       col = c("orange", "purple"), 
       pch = 16)
dev.off()



##### Figure S6: Imputed and observed data in trait matching space (mammals)
 
png("FigsLast/FigureS6_RLQ.png", 
    500, 500, pointsize = 20)
plot(cord2$RQL1, cord2$RQL2,
     frame = F, 
     pch = 1, 
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "Trait-matching axis 1",
     ylab = "Trait-matching axis 2",
     col = scales::alpha("purple", 0.4))# mammals

points(RLQ$mR[,1:2],
       col = "black", 
       pch = "x")
legend("topleft", 
       c("imputed", "observed"), 
       bty = "n", 
       cex = 0.7,
       col = c( "purple"), 
       pch =c("o", "x"))
dev.off()

##### Figure S7: Imputed and observed data in trait matching space (palms)

png("FigsLast/FigureS7_RLQ.png", 
    500, 500, pointsize = 20)
plot(cord$RQL1, cord$RQL2,
     frame = F, 
     pch = 1, 
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "Trait-matching axis 1",
     ylab = "Trait-matching axis 2",
     col = scales::alpha("orange", 0.4))# palms

points(RLQ$mQ[,1:2],
       col = "black", 
       pch = "x")
legend("topleft", 
       c("imputed", "observed"), 
       bty = "n", 
       cex = 0.7,
       col = c( "orange"), 
       pch =c("o", "x"))
dev.off()



plot(cord$RQL1, cord$RQL2) # palms




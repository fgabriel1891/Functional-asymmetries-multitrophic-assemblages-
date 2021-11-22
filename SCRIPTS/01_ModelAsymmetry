


png("FigsLast/AsymmetrySpace.png",  
    500, 500, pointsize = 20)
plot(MyMetric3$ClusAs~MyMetric3$sizeAs,
     pch = 16, 
     col = "grey70", 
     frame = F, 
     xlim = c(0,1),
     xlab = "Functional size asymmetry", 
     ylab = "Functional spacing asymmetry")
abline(v=0.5, h=0, lty = 2)
dev.off()






MyMetric3$ENV1 <- decostand(MyMetric3$ENV1, "range")
MyMetric3$ENV2 <- decostand(MyMetric3$ENV2, "range")
MyMetric3$sizeAs <- decostand(MyMetric3$sizeAs, "range")
MyMetric3$ClusAs <- decostand(MyMetric3$ClusAs, "range")
MyMetric3$ZscorePd <- decostand(MyMetric3$ZscorePd, "range")
MyMetric3$ZscoreMd <- decostand(MyMetric3$ZscoreMd, "range")
MyMetric3$ZscoreP <- decostand(MyMetric3$ZscoreM, "range")
MyMetric3$ZscoreM <- decostand(MyMetric3$ZscoreP, "range")
MyMetric3$Q <- decostand(MyMetric3$Q, "range")
MyMetric3$NODF <- decostand(MyMetric3$NODF, "range")

MyMetric3$INBz <- INB$INBsc[match(paste0(MyMetric3$X1, "_",MyMetric3$X2),
                                  paste0(INB$x, "_",INB$y))]


MyMetric3$INBdz <- INB$INBdZsc[match(paste0(MyMetric3$X1, "_",MyMetric3$X2),
                                     paste0(INB$x, "_",INB$y))]

MyMetric3$INBz <- vegan::decostand(MyMetric3$INBz, "range", na.rm = T)
MyMetric3$INBdz <- vegan::decostand(MyMetric3$INBdz, "range", na.rm = T)


plot(MyMetric3$INBz,MyMetric3$INBdz)
abline(h = 0, v = 0, lty = 2)

plot(MyMetric3$INBz~MyMetric3$ENV1)

boxplot(MyMetric3$INBz~MyMetric3$Dominion)

MyMetric3$Dominion <- neotropic$Dominions[match(MyMetric3$Prov, neotropic$Province_1)]

theme_set(theme_sjplot())
boxplot(MyMetric3$ZscoreM,MyMetric3$ZscoreP)

plot_model(MixModelZm, "pred", show.data = T, grid = T)
plot_model(MixModelZp, "pred", show.data = T, grid = T)
plot_model(MixModelZpd, "pred", show.data = T, grid = T)
plot_model(MixModelZmd, "pred", show.data = T, grid = T)

tab_model(MixModelZm,MixModelZp,MixModelZpd,MixModelZmd)



plot_model(MixModelNODF, "pred", show.data = T, grid = T)
plot_model(MixModelQ, "pred", show.data = T, grid = T)

plot_model(MixModelSA, "pred", show.data = T, grid = T)
plot_model(MixModelCA, "pred", show.data = T, grid = T)

tab_model(MixModelSA,MixModelCA,MixModelNODF,MixModelQ)


length(unique(mamPA_trait2$spName))
length(unique(mamPA_trait$spName))


unique(mamPA_trait$Family[!mamPA_trait$spName %in% mamPA_trait2$spName])




axis(1, labels = F)
par(mar = c(0,0,0,0))
plot(neotropic, col =f(as.numeric(
    as.factor(neotropic$Dominions)),7, "Set2"))

data.frame(neotropic$Dominions,
f(as.numeric(
    as.factor(neotropic$Dominions)),7, "Set2"))



MixModelSA <- lme4::lmer(INBz~ ENV1 + ENV2  + (1|Dominion),REML =F, MyMetric3)
MixModelCA <- lme4::lmer(INBdz~ ENV1 + ENV2  + (1|Dominion),REML =F, MyMetric3)
MixModelNODF <- lme4::lmer(NODF~ ENV1 + ENV2  + (1|Dominion),REML =F, MyMetric3)
MixModelQ <- lme4::lmer(Q~ ENV1 + ENV2  + (1|Dominion),REML =F, MyMetric3)


tab_model(MixModelCA,MixModelSA)


sd(c(lme4::ranef(MixModelSA)$Prov)$`(Intercept)`)
sd(c(lme4::ranef(MixModelCA)$Prov)$`(Intercept)`)


sd(c(lme4::ranef(MixModelNODF)$Prov)$`(Intercept)`)
sd(c(lme4::ranef(MixModelQ)$Prov)$`(Intercept)`)




MixModelZp <- lme4::lmer(ZscoreP~ ENV1 + ENV2  + (1|Dominion),REML =F, MyMetric3)
MixModelZm <- lme4::lmer(ZscoreM~ ENV1 + ENV2  + (1|Dominion),REML =F, MyMetric3)
MixModelZpd <- lme4::lmer(ZscorePd~ ENV1 + ENV2  + (1|Dominion),REML =F, MyMetric3)
MixModelZmd <- lme4::lmer(ZscoreMd~ ENV1 + ENV2  + (1|Dominion),REML =F, MyMetric3)




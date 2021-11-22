
INB$palmRich
colnames(mammalStTrait)[3] <- "taxa"
colnames(palmStTrait)[3] <- "taxa"



AssZscore<-function(INB,
                    INBd,
                    palmStTrait, 
                    mammalStTrait,
                    mammRich,
                    palmRich){
  
  
  AsDis <- replicate(999,log(NicheBreadth(palmStTrait[sample(1:nrow(palmStTrait), 
                                                             palmRich),])/
                               NicheBreadth(mammalStTrait[sample(1:nrow(mammalStTrait), 
                                                                 mammRich),])))
  
  ClDis <- replicate(999,log(NicheDispersion(palmStTrait[sample(1:nrow(palmStTrait), 
                                                                palmRich),])/
                               NicheDispersion(mammalStTrait[sample(1:nrow(mammalStTrait), 
                                                                    mammRich),])))
  
  
  
  INBz <- (INB-mean(AsDis, na.rm = T))/sd(AsDis, na.rm = T)
  INBdz <- (INBd-mean(ClDis, na.rm = T))/sd(ClDis, na.rm = T)
  
  return(c(INBz,INBdz))
  
}

zsc <-c()
 
for(i in 1294:nrow(INB)){
  x<-i
  print(paste("doing", i))
  zsc[[i]] <-AssZscore(INB$INB[x], INB$INBd[x],
            palmStTrait,mammalStTrait,
            INB$mammRich[x],INB$palmRich[x])
  cat("\014")
  }
#zsc[[1293]] <- c("NA", "NA")

zsc1 <- as.numeric(unlist(zsc))
evens <- function(x) subset(x, x %% 2 == 0)

INB$INBdZsc <- zsc1[evens(1:length(zsc1))]
INB$INBsc <-zsc1[-evens(1:length(zsc1))]






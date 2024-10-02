# 1 October 2024 - the "#" allows you to note to yourself and you can run a whole script while
## not having the machine get interrupted by your notes
# always date your script, and put your initials and project near the heading

# Premer, M.I. 
## SICOMP exploratory analysis and summaries - data collected by T. Locke and S. Curry in Summer 2024


dev.off() # clear the plotting
rm(list=ls()) # clear the workspace

require(tidyr) # if you don't have this and the following packages, need to use "install.packages("tidyr")" here, and so forth
require(dplyr)
require(ggplot2)
require(lattice)
require(nlme)


si <- read.csv("~/Google drive/My drive/Research/SICOMP/sicomp.csv") # read in data as a csv
head(si) # look at it
#ws <- dplyr::filter(si,SPECIES=="WS")
xyplot(HT_24~DBH_24|SPECIES, # with lattice package, create a quick scatter plot of HT~DBH across species
       type=c("p","smooth"),
       data=si)

ht.mod2 <- nlme(HT_24~4.5+exp(a+b/(DBH_24+1)), # nonlinear mixed model, fitting a different HT~DBH based on species
                data=si,
                fixed=a+b~1,
                random=a+b~1|SPECIES,
                na.action=na.pass,
                start=c(a=4.5,b=-6),
                control=nlmeControl(returnObject = TRUE,
                                    msMaxIter = 10000,
                                    maxIter = 5000))
summary(ht.mod2) # estimate global coefficients
ranef(ht.mod2) # estimate species specific modifiers to the global average through the random effect

trt <- read.csv("~/Google drive/My drive/Research/SICOMP/trt.csv") # bring in the treatment data frame
si <- left_join(si,trt) # append the treatment list to the treelist for ease of summaries
si$treatment <- as.factor(si$treatment)
si$SPECIES <- as.factor(si$SPECIES)
si$fit <- predict(ht.mod2,si,na.action="na.pass") # impute heights for all trees that don't have a mesaured height
                                                    ## for simplicity, just use modeled heigths for now, ignore mesaured height 
                                                     ### just a RMSE of ~5 ft, not too bad? 

xyplot(fit~DBH_24|SPECIES,data=si,
       type=c("p","smooth"))
si$vol <- (si$DBH_24/12)^1/3*si$fit # crude way to estimate volume
xyplot(vol~DBH_24|SPECIES,data=si)#,
       #type=c("p","smooth"))

si$DBH_24[is.na(si$DBH_24)] <- 0 # get rid of NAs by turning them to 0s
si <- filter(si,DBH_24>0) # get rid of trees that are dead or missing DBH 

dollars <- read.csv("~/Google Drive/My Drive/Research/SICOMP/SIComp_Values.csv") # this is for a quick exercise to estimate the
                                                                                ## financial efficicacy of WS plantation here, 
                                                                              ### stumpage values and cord/tons conversions come from Maine Forest Service 
si <- left_join(si,dollars)
head(si)
require(MEForLab) # this is a package I built locally, you can build locally with the following
#install.packages("devtools")
#devtools::install_github("piketremor/MEForLab"") # just run lines 61 and 62 without the "#" 
require(MEForLab) # load the package
packageDescription('MEForLab')
si.sum <- si%>% # summarize the data frame by calculating standard structural metrics, and apply the cord/ton conversions, apply 2022 stumnpage rates
  mutate(tf=10,
         ba=DBH_24^2*0.005454,na.action="na.omit",
         cords = vol/79,
         tree.tons = cords*Tons,
         tree.value = tree.tons*Ton.Value)%>%
  group_by(plot.id)%>%
  summarize(bapa = sum(tf*ba),
            tpa = sum(tf),
            qmd = qmd(bapa,tpa),
            rd = bapa/sqrt(qmd),
            vol.pa = sum(vol*tf),
            cords.pa = sum(cords*tf),
            tons.pa = sum(tree.tons*tf),
            value.pa = sum(tree.value*tf))
print(si.sum)
si.value <- si.sum%>%
  left_join(.,trt)%>%
  group_by(treatment)%>%
  summarize(npv = mean(value.pa))

h40.frame <- si%>% # need to calculate H40 to estimate si index, to brute force it, select the 4 largest trees by DBH per species per plot 
  group_by(plot.id,SPECIES)%>% # H40 is the average height of the 40 largest DBH trees per acre, so in a 0.1 acre plot, you can just take the 4 largest, but need to determine which species is dominant first
  top_n(.,4,wt=DBH_24)

si.props <- si%>%
  mutate(tf=10,
         ba = DBH_24^2*0.005454)%>%
  group_by(plot.id,SPECIES)%>%
  summarize(sp.bapa = sum(ba*tf),
            sp.tpa = sum(tf))
si.alls <- si.sum%>%
  left_join(.,si.props)%>%
  mutate(props = sp.bapa/bapa)%>%
  group_by(plot.id)%>%
  top_n(.,n=1,wt=props)

sia <- si.alls[c(1,13,10)]
h40.frame2 <- h40.frame%>%
  left_join(.,sia)
h40.frame2$props[is.na(h40.frame2$props)] <- 0
h40.3 <- dplyr::filter(h40.frame2,props>0)

h40.sum <- h40.3%>%
  group_by(plot.id,SPECIES)%>%
  summarize(h40=mean(fit))

h40.sum <- rename(h40.sum,Dominant.Species = SPECIES)
sicomp.summary <- left_join(si.sum,h40.sum)

# just for fun, SI50 is supposed to be a good predictor of volume.. how well does it hold here? 
plot(sicomp.summary$h40,sicomp.summary$cords.pa,ylim=c(0,100),
     xlim=c(0,100))
lm1 <- lm(cords.pa~h40,data=sicomp.summary)
summary(lm1) # not great.... 
abline(lm1)
sicomp.summary$value <- sicomp.summary$cords.pa*50
lm2 <- lm(value~h40,data=sicomp.summary)
plot(sicomp.summary$h40,sicomp.summary$value,ylim=c(0,5000),
     xlim=c(0,100))
abline(lm2)
summary(lm2)

si20 <- left_join(sicomp.summary,trt)
boxplot(cords~treatment,data=si20)
boxplot(rd~treatment,data=si20)


treat.level <- si20%>%
  group_by(treatment)%>%
  summarize(bapa = mean(bapa),
            tpa = mean(tpa),
            qmd = mean(qmd),
            rd = mean(rd),
            vol.pa = mean(vol.pa),
            h40 = mean(h40),
            cords = mean(cords))


tour <- dplyr::filter(si,plot.id=="C3"|plot.id=="HC2"|plot.id=="LC3")
tour <- tour%>%
  mutate(dbh.class = 2*as.integer((tour$DBH_24+(2/2))/2),
         ba = DBH_24^2*0.005454,
         tf = 10,
         bapa = ba*tf)
head(tour)         

#write.csv(si20,"~/Desktop/SI2020.csv")

#png("~/Desktop/SIComp_DDTour.png",units='in',height=6,width=15,res=1000)

ggplot(tour,aes(x=dbh.class,y=bapa,fill=SPECIES))+
  geom_bar(stat="identity")+
  facet_wrap(~plot.id)+
  labs(x="DBH class (in.)",y="Basal area per acre")+
  scale_x_continuous(breaks = c(0,2,4,6,8))+
  theme_set(theme_bw(18))

dev.off()


# equivalence test here


spruce <- dplyr::filter(si,SPECIES=="WS")

spruce$HT_24[is.na(spruce$HT_24)] <- 0
spruce.check <- dplyr::filter(spruce,HT_24>0)
plot(spruce.check$HT_24,spruce.check$fit)
spruce.check$resid <- spruce.check$HT_24-spruce.check$fit
mean(spruce.check$resid)
length(spruce.check$HT_24)
sd(spruce.check$resid)/sqrt(92)
plot(spruce.check$HT_24,spruce.check$fit)
ts <- mean(spruce.check$resid)/(sd(spruce.check$resid)/sqrt(92))
spruce.check$scaled.obs <- spruce.check$HT_24-mean(spruce.check$HT_24)
spruce.check$scaled.res <- spruce.check$resid-mean(spruce.check$resid)
plot(spruce.check$scaled.obs,spruce.check$fit,
     xlab="Measured Height (ft.)",
     ylab="Estimated Height (ft.)")
abline(v=0,col="gray80")
mod1 <- lm(fit~scaled.obs,data=spruce.check)
abline(mod1)
abline(mod1$coefficients[1],mod1$coefficients[2],col="darkred",lwd=3)
abline(confint(mod1)[1],mod1$coefficients[2],lty=2)
abline(confint(mod1)[3],mod1$coefficients[2],lty=2)
abline(h=mean(spruce.check$HT_24)+mean((spruce.check$HT_24*0.25)),col="gray80")
abline(h=mean(spruce.check$HT_24)-mean((spruce.check$HT_24*0.25)),col="gray80")
abline(mod1$coefficients[1],confint(mod1)[2],lty=1,col="darkgreen")
abline(mod1$coefficients[1],confint(mod1)[4],lty=1,col="darkgreen")

tost(x=spruce.check$fit,y=spruce.check$HT_24,
     epsilon=1,
     alpha=0.95,
     var.equal=TRUE)

text(-2,33,"Null hypothesis of difference is rejected at p = 0.17")


#abline(22.18,0.6859,col="darkred",lwd=3)
mean(spruce$HCB_24)*0.25
#abline(h=8.75+2.18)
#abline(h=8.75-2.18)
#confint(mod1)
#abline(20.18,(0.68+0.068))
#abline(20.18,(0.68+0.068))


xyplot(fit~DBH_24|SPECIES,data=spruce)
spruce$fvs.ht <- mapply(wykoff.ht,
                        SPP=spruce$SPECIES,
                        DBH=spruce$DBH_24)
head(spruce)
spruce$resid <- spruce$fit-spruce$fvs.ht
xyplot(resid~DBH_24|SPECIES,data=spruce)

xyplot(fit~DBH_24|SPECIES,data=spruce)
plot(spruce$HT_24,spruce$fvs.ht,
     ylim=c(0,50),
     xlim=c(0,50),
     abline(0,1),
     col="darkgreen",
     pch=15,
     ylab="Modeled Height (ft.)",
     xlab="Observed Height (ft.)")
points(spruce$HT_24,spruce$fit,col="red",
       pch=19)

fvs.mod <- lm(fvs.ht~HT_24,data=spruce)
fit.mod <- lm(fit~HT_24,data=spruce)
summary(fvs.mod)
summary(fit.mod)

equivalence::tost(x=spruce$HT_24,y=spruce$fit,epsilon=1,paired=TRUE,conf.level=0.95,var.equal=TRUE)
equivalence::tost(x=spruce$HT_24,y=spruce$fvs.ht,epsilon=1,paired=TRUE,conf.level=0.95,var.equal=TRUE)



#### done









head(sicomp.summary)
steinman.site(WS,30,24)

si$plot.id <- as.factor(si$plot.id)

str(si)
si$SPECIES <- as.factor(si$SPECIES)


require(dplyr)
si.s <- si%>%
  group_by(plot.id,SPECIES)%>%
  arrange(desc(DBH_24),by.group=TRUE)

si.s


plot(density(si.sum$rd))
trt <- read.csv("~/Desktop/trt.csv")
si.sums <- left_join(si.sum,trt)
si.sums
boxplot(vol.pa~treatment,data=si.sums)

si.sum.sp <- si%>%
  mutate(tf=10,
         ba=DBH_24^2*0.005454,na.action="na.omit")%>%
  group_by(plot.id,SPECIES)%>%
  summarize(sp.bapa = sum(tf*ba),
            sp.tpa = sum(tf),
            sp.vol = sum(tf*vol))
boxplot(sp.vol~SPECIES,data=si.sum.sp)


si.sum.sp

si.sum.frame <- si.sum.sp%>%
  left_join(.,si.sum)%>%
  mutate(iv = ((sp.bapa/bapa)+(sp.tpa/tpa))/2)%>%
  ungroup()%>%
  filter(.,SPECIES=="WS")


keep <- si.sum.frame[c(1,11)]

df <- left_join(si.sum,keep)
df$iv[is.na(df$iv)] <- 0

df

df$cords <- df$vol.pa/79
plot(df$iv,df$cords)
ks <- lm(cords~iv,data=df)
summary(ks)
plot(df$iv,df$cords)

df

plot(si.sum$bapa,si.sum$tpa)
s

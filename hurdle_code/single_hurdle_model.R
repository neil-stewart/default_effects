#-----------------------------------------------------------------------------------
# single hurdle model
#-----------------------------------------------------------------------------------
rm(list = ls(all=TRUE))
library(data.table)
library(mhurdle)

# matched data following Simmon and Hopkins (2005) method
data <- fread("simmons_hoplins_matched_data.csv") 

# exclude card-months with zero balance
dat <- data[repayment>=0 & !is.na(score) & non_bal==0]

# note. "dummy" has a value of 1 for treatment group and 0 for control group
formula <- repayment~ tenure + switcher*dummies | 
                 bal + utilisation + amt_spend + mer_APR + cash_APR + score + switcher*dummies

# assuming independence of two distrubances as in Cragg (1971)
model <- mhurdle(formula, data=dat, dist="ln", h2=FALSE, corr=FALSE, method="bhhh") 
summary(model)

#---------------------
# predictions
#---------------------
new.data <- expand.grid(repayment=100, bal=mean(data$bal), utilisation=mean(data$utilisation), amt_spend=mean(data$amt_spend), 
                      mer_APR=mean(data$mer_APR), cash_APR=mean(data$cash_APR), score=mean(data$score, na.rm=TRUE), 
                      tenure=mean(data$tenure), dummies=as.factor(seq(-6,5)), switcher=c(0,1))
new.data <- data.frame(new.data)

pred <- predict(model, newdata=new.data) 
zero <- pred[,"zero"]
pos  <- pred[,"pos"]

#--------------------------------------
# bootstrap for having prediction CIs
#--------------------------------------
library(plyr)
pbar <- create_progress_bar("text")
pbar$init(1000)


# pairs of treatment and control accounts
res <- fread("simmons_hoplins_matched.csv")
pairs <- res[AccountNumber1 %in% dat[DD_type=="min", AccountNumber],.(AccountNumber1, AccountNumber2)]

# resample matched pairs
N <- 1
Nall <- 1

while(N<1000)
{
  skip_to_next <- FALSE
  print(c(N, Nall))
  Nall <- Nall + 1
  rows <- sample(1:nrow(pairs), nrow(pairs), replace=TRUE)
  sampled <- pairs[rows]
  ac1 <- data.table(AccountNumber=sampled$AccountNumber1)
  ac2 <- data.table(AccountNumber=sampled$AccountNumber2)
  setkey(ac1, AccountNumber)
  setkey(ac2, AccountNumber)
  setkey(dat, AccountNumber)
  b1 <- merge(dat, ac1)
  b2 <- merge(dat, ac2)
  b.dat <- rbind(b1, b2)
  model <- tryCatch(mhurdle(formula, data=b.dat, 
                     dist="ln", h2=FALSE, corr=FALSE, method="bhhh"),
                     error=function(e){skip_to_next <<- TRUE})
  if(skip_to_next==TRUE){next}else{
  pred <- predict(model, newdata=new.data)
  zero <- rbind(zero, pred[,"zero"])
  pos  <- rbind(pos,  pred[,"pos"])
  N <- nrow(pos)
  pbar$step()
  }
}

D1 <- as.numeric(zero[1,])
L1 <- unlist(apply(zero, 2, function(x){quantile(x, .025)}))
U1 <- unlist(apply(zero, 2, function(x){quantile(x, .975)}))

D2 <- as.numeric(pos[1,])
L2 <- unlist(apply(pos, 2, function(x){quantile(x, .025)}))
U2 <- unlist(apply(pos, 2, function(x){quantile(x, .975)}))

library(ggplot2)

pdata <- data.table(D=D1, L=L1, U=U1, x=as.numeric(as.character(new.data$dummies)), 
                    Type=c(rep("Synthetic Control",12), rep("Switcher",12)))

ggplot(pdata, aes(x=x, y=D, ymin=L, ymax=U, color=Type, fill=Type))+geom_point()+geom_line()+
  geom_ribbon(aes(col=NULL), alpha=0.1) + 
  scale_x_continuous(name="Months from (Synthetic) Switching", breaks=seq(-6,5))+
  theme_bw(base_size=23)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), strip.text.x=element_text(size=23, angle=0),
       strip.background=element_blank(), axis.line=element_line(colour="black", size=0.5),
       panel.border=element_blank(), legend.position=c(0.9,0.9), strip.placement = "outside")+
  labs(x="Months from First Autopay", y="P(Repayment = 0)", col="", size=2)+ guides(fill=FALSE,outline=FALSE)
ggsave("zero_single_hurdle_mean_covariates_simmons.pdf", width=15, height=7)

pdata <- data.table(D=D2, L=L2, U=U2, x=as.numeric(as.character(new.data$dummies)), 
                    Type=c(rep("Synthetic Control",12), rep("Switcher",12)))

ggplot(pdata, aes(x=x, y=D, ymin=L, ymax=U, color=Type, fill=Type))+geom_point()+geom_line()+
  geom_ribbon(aes(col=NULL), alpha=0.1) + 
  scale_x_continuous(name="Months from (Synthetic) Switching", breaks=seq(-6,5))+
  theme_bw(base_size=23)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), strip.text.x=element_text(size=23, angle=0),
       strip.background=element_blank(), axis.line=element_line(colour="black", size=0.5),
       panel.border=element_blank(), legend.position=c(0.1,0.2), strip.placement = "outside")+
  labs(x="Months from First Autopay", y="E(Repayment | Repayment >0)", col="", size=2)+ guides(fill=FALSE,outline=FALSE)
ggsave("payment_single_hurdle_mean_covariates_simmons.pdf", width=15, height=7)




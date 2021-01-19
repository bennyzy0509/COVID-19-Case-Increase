library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(forecast)
library(glmnet)
library(leaps)
library(caret)
library(stargazer)
library(gridExtra)
library(gam)
library(mgcv)
library(randomForest)
RNGkind(sample.kind = "Rounding")
################Part I: Data Cleaning################
covid <- read.csv("case.csv",header=T)
demographic <- read.csv("Social_Vulnerability_Index_2018_-__United_States__county.csv", header = T)
mobility <- read_excel("mobility_us.xlsx")
unemp_race <- read_excel("unemployment_racegroup.xlsx")

covid[is.na(covid)] <- 0 # replace NA with 0
covid<-subset(covid, state!="AS" & state!="VI" & state!="MP"& state!="GU" & state!="PR") # we drop the rows because insular area are out of our research scope
# aggregate by state
demo <- demographic %>%
  group_by(ST_ABBR) %>%
  summarise(ppov = mean(EP_POV), pincome = mean(EP_PCI), pnoveh = mean(EP_NOVEH), puninsur = mean(EP_UNINSUR), page65 = mean(EP_AGE65))

covid= merge(covid, demo, by.x = "state", by.y = "ST_ABBR",all.x = TRUE, all.y = FALSE) # merge covid and demographic datasets

# aggregate by date
covidbydate <- covid %>%
  group_by(date) %>%
  summarise(case = sum(positive),case_incrs = sum(positiveIncrease), death = sum(death),
            test = sum(totalTestResults),hos = sum(hospitalizedCurrently),ppov = mean(ppov), pincome = mean(pincome), 
            pnoveh = mean(pnoveh), puninsur = mean(puninsur), page65 = mean(page65))

mobility <- mobility %>%
  group_by(date) %>%
  summarise(residential = mean(residential), grocery = mean(grocery),park = mean(park), retail = mean(retail), 
            transit = mean(transit),workplace = mean(workplace))

covidbydate$date = as.Date(as.character(covidbydate$date),"%Y%m%d")
mobility$date = as.Date(mobility$date)
unemp_race$date = as.Date(unemp_race$date)

covidbydate = merge(covidbydate, mobility, by.x = c("date"), by.y = c("date"),all.x = TRUE, all.y = FALSE) # merge covid and mobility datasets
covidbydate = merge(covidbydate, unemp_race, by.x = c("date"), by.y = c("date"),all.x = TRUE, all.y = FALSE)

# create new variables: outbrk_days and test_pos
covidbydate$outbrk_days = round(as.numeric(difftime(covidbydate$date, "2020-01-20", units = "days")),0) # days since the first day of outbreak
covidbydate$test_pos = (covidbydate$case/covidbydate$test)*100

# clean NA's
summary(covidbydate)
covidbydate <- subset(covidbydate, date >= "2020-02-15" & date <= "2020-07-21") # limit the time range from 2/15 to 7/21
covidbydate$park[is.na(covidbydate$park)] <- 0 # replace 2 NAs with 0
# multiple imputation to fill in NA for race group
attach(covidbydate)
lm <- lm(black_deaths ~ outbrk_days,covidbydate)
covidbydate$black_deaths[is.na(covidbydate$black_deaths)] <- predict(lm, list(covidbydate$outbrk_days[is.na(covidbydate$black_deaths)]))
lm <- lm(white_deaths ~ outbrk_days+case+death,covidbydate)
covidbydate$white_deaths[is.na(covidbydate$white_deaths)] <- predict(lm, list(covidbydate$outbrk_days[is.na(covidbydate$white_deaths)]))
covidbydate$black_deaths[covidbydate$black_deaths <0] <- 0 # since death won't be negative, replace negative values to 0
covidbydate$white_deaths[covidbydate$white_deaths <0] <- 0

covidbydate <- covidbydate[, c(1,21,2:5,22,6,12:18,7:11,19:20)] # reorder the columns
covidbydate = data.frame(covidbydate)

write.csv(covidbydate,file = "finaldata_coronavirus.csv")

################Part II: Data Visualization################
## summary statistics table
stargazer(covidbydate, type = "text", title="Linear Regression Results", align=TRUE,single.row=TRUE)

## daily new case 
avgnewcases <- covidbydate %>%
  dplyr::mutate(
    new_conf_03da = zoo::rollmean(case_incrs, k = 3, fill = NA),
    new_conf_07da = zoo::rollmean(case_incrs, k = 7, fill = NA)) %>% 
  dplyr::ungroup()

df <- data.frame(grp=covidbydate$date,val=covidbydate$case_incrs)
df2 <- data.frame(grp=avgnewcases$date,val=avgnewcases$new_conf_07da)

names(avgnewcases)
ggplot(df, aes(x=grp, y=val)) + geom_bar(stat="identity",fill="steelblue", width = 1, colour="white") + 
  geom_line(data=df2, aes(x=grp, y=val), colour="red",size = 1)+theme_minimal()+
  labs(title = 'Daily New Covid-19 Cases from Feb to Jul in the United States',
       x = 'Month',
       y = 'Daily New Cases') + theme(
         plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
       )
## scatter plots
attach(covidbydate)
plot(covidbydate[2:14], main = "Relationship for Each Pair of Variables")

## visualize relationship between daily new cases and selected variables
par(mfrow=c(2,3))
scatter.smooth(outbrk_days, case_incrs, ylab = "daily new cases", xlab = "outbreak days", col="gray")
scatter.smooth(case, case_incrs, ylab = "daily new cases",xlab = "cumulative cases",col="gray")
scatter.smooth(death, case_incrs, ylab = "daily new cases", col="gray")
scatter.smooth(test_pos, case_incrs, ylab = "daily new cases",xlab = "test positive rate",col="gray")
scatter.smooth(hos, case_incrs, ylab = "daily new cases",xlab = "current hospitalization",col="gray")
scatter.smooth(black_deaths, case_incrs, ylab = "daily new cases", xlab = "deaths of American Black", col="gray")
mtext("Smoothing Scatter Plot for Selected Variables",side = 3, line = -2.5, font = 4, outer = TRUE, cex = 1)

################Part III: Modeling################
formula = case_incrs~outbrk_days+case+death+test+test_pos+hos+residential+grocery+park+retail+transit+workplace+ppov+punemp+pincome+pnoveh+puninsur+page65+black_deaths+white_deaths
## Randomly Split into trainning and test datasets
set.seed(1)
row.number = sample(1:nrow(covidbydate), 0.8 * nrow(covidbydate))
train = covidbydate[row.number, ]
test = covidbydate[-row.number, ]
covid.train=covidbydate[row.number,"case_incrs"]
covid.test=covidbydate[-row.number,"case_incrs"]

## Linear Regression Models
mod1 = regsubsets(formula,train,nvmax =21)
res.sum <- summary(mod1)
data.frame(
  Adj.R2 = which.max(res.sum$adjr2),
  CP = which.min(res.sum$cp),
  BIC = which.min(res.sum$bic)
)
coef(mod1,13)
mod2 = lm(case_incrs~outbrk_days+death+test+test_pos+hos+grocery+retail+transit+workplace+ppov+pnoveh+black_deaths+white_deaths,train)
stargazer(mod2, type = "text", title="Linear Regression Results", align=TRUE,single.row=TRUE)

par(mfrow=c(2,2))
plot(mod2)

pred2.train = predict(mod2,newdata=train)
data.frame(
  MSE.Train = (RMSE(pred2.train, covid.train)^2),
  R2.Train = R2(pred2.train, covid.train)) # 12713787

pred2.test = predict(mod2,newdata=test)
data.frame(
  MSE.Test = (RMSE(pred2.test, covid.test)^2),
  R2.Test = R2(pred2.test, covid.test)) #21133054

## non-linear models
x = death
y = case_incrs
a = ggplot(covidbydate, aes(x, y)) + geom_point() + geom_smooth(method = "lm", formula = y ~x) + labs(title = "Linear Regression Model", x = 'Death', y = 'Daily New Cases')
b = ggplot(covidbydate, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x)) + labs(title = "Generalized Additive Model", x = 'Death', y = 'Daily New Cases')
grid.arrange(a,b)

x = hos
y = case_incrs
a = ggplot(covidbydate, aes(x, y)) + geom_point() + geom_smooth(method = "lm", formula = y ~x) + labs(title = "Linear Regression Model", x = 'Current Hospitalization', y = 'Daily New Cases')
b = ggplot(covidbydate, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x)) + labs(title = "Generalized Additive Model", x = 'Current Hospitalization', y = 'Daily New Cases')
grid.arrange(a,b)


## GAM
## fitting natural cubic spline
mod3 <- gam(case_incrs ~ s(outbrk_days)+s(case)+s(death)+s(test)+s(test_pos) +
              s(retail)+s(transit)+s(workplace)+s(black_deaths)+s(white_deaths)+hos+grocery+
              ppov+pnoveh,data=train)
summary(mod3)  # check the significance of smooth term

mod3 <- update(mod3,.~.-s(test_pos)-s(retail)-s(transit)-s(white_deaths)-s(black_deaths)+test_pos+retail+transit+white_deaths+black_deaths)
summary(mod3) 

pred3.train = predict(mod3,newdata=train)
data.frame(
  MSE.Train = (RMSE(pred3.train, covid.train)^2),
  R2.Train = R2(pred3.train, covid.train))# MSE 2081354

pred3.test = predict(mod3,newdata=test)
data.frame(
  MSE = (RMSE(pred3.test, covid.test)^2),
  R2 = R2(pred3.test, covid.test))

anova(mod2, mod3, test="Chisq") # significance for the second one

## ridge
set.seed(1)
x = model.matrix(formula, covidbydate)[, -1]
y = covidbydate$case_incrs
train = sample(1:nrow(covidbydate), 0.8 * nrow(covidbydate))
test = (-train)
y.train = y[train]
y.test = y[test]
mod4 = glmnet(x[train, ], y[train], alpha = 0, 
              lambda = grid)
pred4 = predict(mod4, s = 4, newx = x[test, ])
mean((pred4 - y.test) ^ 2)
cv.out = cv.glmnet(x[train, ], y[train], alpha = 0)
par(mfrow=c(1,1))
plot(cv.out)
bestlam = cv.out$lambda.min
bestlam
pred4.train = predict(mod4, s = bestlam, newx = x[train, ])
pred4 = predict(mod4, s = bestlam, newx = x[test, ])
mean((pred4.train - y.train) ^ 2)
mean((pred4 - y.test) ^ 2)
out = glmnet(x, y, alpha = 0)
predict(out, type = "coefficients", s = bestlam)[1:20, ]
R2(pred4, covid.test)

## lasso
set.seed (1)
mod5 = glmnet(x[train, ], y[train], alpha = 1, lambda = grid)
plot(mod5)
cv.out = cv.glmnet(x[train, ], y[train], alpha = 1)
plot(cv.out)
bestlam = cv.out$lambda.min
pred5.train = predict(mod5, s = bestlam, newx = x[train, ])
pred5 = predict(mod5, s = bestlam, newx = x[test, ])
mean((pred5.train - y.train) ^ 2)
mean((pred5 - y.test) ^ 2)
out = glmnet(x, y, alpha = 1, lambda = grid)
lasso.coef = predict(out, type = "coefficients", s = bestlam)[1:20, ]
lasso.coef
R2(pred5, covid.test)

## random forest
set.seed (1)
install.packages("devtools")
library(devtools)
devtools::install_github('skinner927/reprtree')
library(reprtree)
covid.test = covidbydate[-train, "case_incrs"]
mod6 = randomForest(formula, data = covidbydate, subset = train, 
                    mtry = 6, importance = TRUE)
mod6
plot(mod6)
reprtree:::plot.getTree(mod6) #tree plot
pred6.train = predict(mod6, newdata = covidbydate[train, ])
pred6.test = predict(mod6, newdata = covidbydate[-train, ])
mean((pred6.train - covid.train) ^ 2)
mean((pred6.test - covid.test) ^ 2)
plot(pred6.test, covid.test)
abline(0, 1)
varImpPlot(mod6)
R2(pred6.test, covid.test)

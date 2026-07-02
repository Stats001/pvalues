# R Code extracted from: How to Get P-Values Less Than 0.05 in Biological Research
# (Second Edition) by Michael Walker
# Extracted automatically from the PDF; chapter/section headers preserved as comments.


#==============================================================================
# How to Get P-Values Less Than 0.05 in Biological Research
#==============================================================================

#==============================================================================
# Preface
#==============================================================================

# ---- A note on what this book is and isn't ----

# ---- Website and companion book ----

# ---- R packages used in the book ----

install.packages(c("ggplot2", "Rfit", "RobStatTM", "car", "vcdExtra", "tibble",
                    "dplyr", "tidyverse", "survival", "survminer", "coin", "rankFD"))

# ---- Dedication ----

# ---- Acknowledgements ----

# ---- About the author ----

# ---- Copyright ----

#==============================================================================
# 1 Avoid t-tests: Control for variables using multiple regression or ANOVA
#==============================================================================

# ---- 1.1 Short story: Of mice and marriage ----

library(tidyverse) # For graphs and data manipulation
mutation.weight.data = tibble::tribble(
~Treatment, ~Sex, ~Grams,
"WildType", "Female", 24.2,
"WildType", "Female", 27.9,
"WildType", "Female", 28.8,
"WildType", "Male", 32.7,
"WildType", "Male", 34.5,
"WildType", "Male", 38.8,
"Mutant", "Female", 20.5,
"Mutant", "Female", 24.6,
"Mutant", "Female", 25.1,
"Mutant", "Male", 29.8,
"Mutant", "Male", 30.9,
"Mutant", "Male", 33.5
)
# Create the plot "Weight of wild type versus mutant mice".
ggplot(mutation.weight.data, aes(y=Grams, x=Treatment)) +
geom_point() +
ggtitle("Weight of wild type versus mutant mice")
t.test(Grams ~ Treatment, var.equal=TRUE, data =
mutation.weight.data)

summary(lm(Grams ~ Treatment, data = mutation.weight.data))

# Create the plot "Weight of wild type versus mutant mice broken out by sex".
ggplot(mutation.weight.data, aes(y=Grams, x=Treatment)) +
geom_point() +
facet_wrap(~Sex) +
ggtitle("Weight of wild type versus mutant mice \nbroken out by sex")
summary(lm(Grams ~ Treatment,
data=subset(mutation.weight.data, Sex=="Male")))

summary(lm(Grams ~ Treatment,
data=subset(mutation.weight.data, Sex=="Female")))

# Do the two-way ANOVA.
# Only change: add Grams ~ Sex + Treatment
summary(lm(Grams ~ Sex + Treatment, data=
mutation.weight.data))

# ---- 1.2 Example: Get a better p-value by controlling for sources of variability in an experiment. ----

pcr.data = tibble::tribble(
~Treatment,
~Column, ~DeltaCT,
"Control",
"A",
11.5,
"LowDose",
"A",
11.6,
"MediumDose",
"A",
11.9,
"Control",
"B",
11.6,
"LowDose",
"B",
11.9,
"MediumDose",
"B",
12.0,
"Control",
"C",
11.9,
"LowDose",
"C",
12.3,
"MediumDose",
"C",
12.8,
"Control",
"D",
11.9,
"LowDose",
"D",
12.5,
"MediumDose",
"D",
13.5,
"Control",
"E",
12.7,
"LowDose",
"E",
11.9,
"MediumDose",
"E",
12.9,
"Control",
"F",
11.9,
"LowDose",
"F",
12.8,
"MediumDose",
"F",
13.0,
"Control",
"G",
12.9,
"LowDose",
"G",
13.1,
"MediumDose",
"G",
13.5,
"Control",
"H",
12.9,
"LowDose",
"H",
13.0,
"MediumDose",
"H",
13.4
)
with(pcr.data, boxplot(DeltaCT ~ Treatment, xlab="Treatment",
ylab = "Delta CT", ylim=c(11,14)))
# Do one-way ANOVA for DeltaCT ~ Treatment
lm.pcr.1factor = lm(DeltaCT ~ Treatment, data= pcr.data)
anova(lm.pcr.1factor)

with(pcr.data,
boxplot(DeltaCT ~ Column, xlab="Column", ylab = "Delta CT",
ylim=c(11,14)))
# Do two-way ANOVA, controlling for Column effect
lm.pcr.2factor = lm(DeltaCT ~ Treatment + Column , data=
pcr.data)
anova(lm.pcr.2factor)

# ---- 1.3 Identify important variables quickly using fractional factorial designs: ----

# ---- 1.4 P-value, null hypothesis, and alternative hypothesis ----

# ---- 1.5 R code and data sets ----

#==============================================================================
# 2 Calculate power and sample size if you want p < 0.05
#==============================================================================

# ---- 2.1 What determines your power and sample size? ----

power.t.test(delta = 4, sd = 2, sig.level = 0.05, power = 0.8)

power.t.test(delta = 2, sd = 2, sig.level = 0.05, power = 0.8)

power.t.test(delta = 2, sd = 3, sig.level = 0.05, power = 0.8)

power.t.test(delta = 6, sd = 4, sig.level = 0.05, n=25)

power.t.test(delta = 6, sd = 4, sig.level = 0.05, power =
0.95)

# ---- 2.2 Software for power and sample size ----

# ---- 2.3 Examples of power and sample size calculations using the base R functions ----

power.t.test(delta = 0.5, sd = 0.5, sig.level = 0.01, power =
0.9)

power.t.test(delta = 0.5, sd = 0.5, sig.level = 0.01, n = 31)

# ---- 2.4 Estimate sample size using data from a pilot study ----

group1=c(1.1, 2.3, 4.8, 5.4, 7.9)
group2=c(3.3, 4.2, 4.7, 7.6, 9.2)
mean(group2)-mean(group1)

sd(group1)

sd(group2)

power.t.test(delta = 1.5, sd = 2.7, sig.level = 0.05,
power=0.8)

# ---- 2.5 Assay variability reduces power, increases sample size, and gives worse p-values ----

power.t.test(delta=5, sd=10, sig.level = 0.05, power=0.8)

power.t.test(delta=5, sd=5, sig.level = 0.05, power=0.8)

#==============================================================================
# 3 Increase sample size last: reduce unexplained variability first
#==============================================================================

# ---- 3.1 Simulations to illustrate the effect of sample size and within-group variability on power ----

set.seed(1234)
placebo= rnorm(10000, 150, 20)
mean(placebo)

sd(placebo)

set.seed(3456)
drug = rnorm(10000, 140, 20)
mean(drug)

sd(drug)

plot(density(placebo), xlim= c(50, 250),ylim=c(0,.025))
lines(density(drug), lty=2)
placebo.sample = sample(placebo, size=30)
drug.sample = sample(drug, size=30)
plot(density(placebo.sample), xlim= c(50, 250),ylim=c(0,
.025))
lines(density(drug.sample), lty=2)
ttest.result = t.test(placebo.sample, drug.sample,
var.equal=TRUE)
ttest.result$p.value

set.seed(1234)
for (i in 1:4) {
placebo.sample = sample(placebo, size=30)
drug.sample = sample(drug, size=30)
ttest.result = t.test(placebo.sample, drug.sample,
var.equal=TRUE)
print(ttest.result$p.value)
}

set.seed(12345)
n = 30
pvalue.list = c()
for (i in 1:1000)
{
placebo.sample = sample(placebo, size=n)
drug.sample = sample(drug, size=n)
pvalue.list[i] = t.test(placebo.sample, drug.sample,
var.equal=TRUE)$p.value
}
hist(pvalue.list, xlim= c(0, 1), breaks=seq(0,1,.05),
ylim=c(0,1000))
set.seed(12345)
n = 50
pvalue.list = c()
for (i in 1:1000)
{
placebo.sample = sample(placebo, size=n)
drug.sample = sample(drug, size=n)
pvalue.list[i] = t.test(placebo.sample, drug.sample,
var.equal=TRUE)$p.value
}
hist(pvalue.list, xlim= c(0, 1), breaks=seq(0,1,.05),
ylim=c(0,1000))
set.seed(1234)
placebo= rnorm(10000, 150, 10)
mean(placebo)

sd(placebo)

set.seed(3456)
drug = rnorm(10000, 140, 10)
mean(drug)

sd(drug)

plot(density(placebo), xlim= c(50, 250), ylim=c(0, .05))
lines(density(drug), lty=2)
set.seed(4321)
placebo.sample = sample(placebo, size=30)
drug.sample = sample(drug, size=30)
plot(density(placebo.sample), xlim= c(50, 250),ylim=c(0, .05))
lines(density(drug.sample), lty=2)
t.test(placebo.sample, drug.sample, var.equal=TRUE)

set.seed(4321)
n = 30
pvalue.list = c()
for (i in 1:1000)
{
placebo.sample = sample(placebo, size=n)
drug.sample = sample(drug, size=n)
pvalue.list[i] = t.test(placebo.sample, drug.sample,
var.equal=TRUE)$p.value
}
hist(pvalue.list, xlim= c(0, 1), breaks=seq(0,1,.05),
ylim=c(0,1000))
sum(pvalue.list < 0.05)/1000

# ---- 3.2 Which of the four factors that affect power should you change to best increase power? ----

# ---- 3.3 Definition of effect size and standardized effect size ----

#==============================================================================
# 4 If you have normally distributed data with outliers: Robust tests
#==============================================================================

# install.packages("RobStatTM")
library(RobStatTM)
library(ggplot2)

# ---- 4.1 Robust methods: Reduce the influence of extreme observations ----

cancer=c(83, 89, 90, 93, 98)
controls=c(99, 100, 103, 104,
141)
Mucin=c(cancer,controls)
Diagnosis=c(rep("Cancer",5), rep("Controls",5))
colon.data=data.frame(Diagnosis, Mucin)
colon.data

with(colon.data, boxplot(Mucin ~ Diagnosis, ylab="Mucin",
ylim=c(50,150)))
aggregate(Mucin ~ Diagnosis, colon.data, mean)

aggregate(Mucin ~ Diagnosis, colon.data, median)

with(colon.data, t.test(Mucin ~ Diagnosis, var.equal=TRUE))

summary(with(colon.data, lmrobM(Mucin ~ Diagnosis)))

anova.2way.data = tibble::tribble(
~Treatment, ~Litter, ~Response,
"Drug", "Litter.1", 4,
"Drug", "Litter.1", 8,
"Drug", "Litter.2", 11,
"Drug", "Litter.2", 12,
"Drug", "Litter.3", 13,
"Drug", "Litter.3", 15,
"Drug", "Litter.4", 20,
"Drug", "Litter.4", 21,
"Placebo", "Litter.1", 10,
"Placebo", "Litter.1", 15,
"Placebo", "Litter.2", 16,
"Placebo", "Litter.2", 17,
"Placebo", "Litter.3", 19,
"Placebo", "Litter.3", 20,
"Placebo", "Litter.4", 21,
"Placebo", "Litter.4", 23)
ggplot(anova.2way.data,
aes(y=Response, x=Treatment)) +
geom_point() + ggtitle("Response by treatment group")
t.test(Response ~ Treatment, var.equal=TRUE, data =
anova.2way.data)

summary(lm (Response ~ Treatment, data = anova.2way.data))

ggplot(anova.2way.data, aes(y=Response, x=Treatment)) +
geom_point() +
facet_wrap(~Litter) +
ggtitle("Response to treatment \nbroken out by Litter") + ylim(0,25)
lm.2way = lm(Response ~ Litter + Treatment, data=
anova.2way.data)
anova(lm.2way)

anova.2way.data.outlier = tibble::tribble(
~Treatment, ~Litter, ~Response,
"Drug", "Litter.1", 4,
"Drug", "Litter.1", 8,
"Drug", "Litter.2", 11,
"Drug", "Litter.2", 12,
"Drug", "Litter.3", 13,
"Drug", "Litter.3", 15,
"Drug", "Litter.4", 20,
"Drug", "Litter.4", 32,
"Placebo", "Litter.1", 10,
"Placebo", "Litter.1", 15,
"Placebo", "Litter.2", 16,
"Placebo", "Litter.2", 17,
"Placebo", "Litter.3", 19,
"Placebo", "Litter.3", 20,
"Placebo", "Litter.4", 21,
"Placebo", "Litter.4", 31)
ggplot(anova.2way.data.outlier, aes(y=Response, x=Treatment))
+ geom_point() +
ggtitle("Response by treatment group \nwith outliers")
t.test(Response ~ Treatment, var.equal=TRUE, data =
anova.2way.data.outlier)

summary(lm (Response ~ Treatment, data =
anova.2way.data.outlier))

ggplot(anova.2way.data.outlier, aes(y=Response, x=Treatment))
+ geom_point() +
facet_wrap(~Litter) +
ggtitle("Response to treatment \nbroken out by Litter, with outliers")
lm.2way.outlier = lm(Response ~ Litter + Treatment,
data= anova.2way.data.outlier)
anova(lm.2way.outlier)

anova(lm(Response ~ Litter * Treatment, data=
anova.2way.data.outlier))

anova.2way.data.outlier.lmrobM = lmrobM(Response ~ Litter +
Treatment,
data = anova.2way.data.outlier)
summary(anova.2way.data.outlier.lmrobM)

anova.2way.data.lmrobM = lmrobM(Response ~ Litter + Treatment,
data = anova.2way.data) # control=cont)
summary(anova.2way.data.lmrobM)

# ---- 4.2 Which method should you use? ----

#==============================================================================
# 5 If you have non-normal data and one explanatory variable: Wilcoxon and Kruskal-Wallis
#==============================================================================

ranked.data
= data.frame(extra = sort(sleep$extra),
ranked.extra = rank(sort(sleep$extra)))
ranked.data

# ---- 5.1 Alternative to two-sample t-test: Wilcoxon rank sum for two independent samples ----

cancer=c(83, 89, 90, 93, 98)
controls=c(99, 100, 103, 104,
141)
Mucin=c(cancer,controls)
Diagnosis=c(rep("Cancer",5), rep("Controls",5))
colon.data=data.frame(Diagnosis, Mucin)
colon.data

with(colon.data, boxplot(Mucin ~ Diagnosis, ylab="Mucin",
ylim=c(50,150)))
aggregate(Mucin ~ Diagnosis, colon.data, mean)

aggregate(Mucin ~ Diagnosis, colon.data, median)

with(colon.data, t.test(Mucin ~ Diagnosis, var.equal=TRUE))

with(colon.data, wilcox.test(Mucin ~ Diagnosis,
conf.int=TRUE))

# ---- 5.2 How do you quantify the difference between two treatment groups using ranks? The Hodges-Lehmann estimator ----

treatment = c(5, 7, 8, 10, 12, 15, 18)
control = c(9, 11, 13, 14, 16, 20, 25)
# Calculate all pairwise treatment - control differences
cross_diffs = outer(treatment, control, "-")
cross_diffs

sort(cross_diffs)

median(cross_diffs)

median(treatment)

median(control)

median(treatment) - median(control)

# ---- 5.3 Sample size for a t-test to give power of 0.8 for the colon data ----

power.t.test(delta = 18.8, sd = 17.8, sig.level = 0.05, power
= 0.8)

# ---- 5.4 Alternative to paired t-test: Wilcoxon signed rank test ----

bears.winter=c(300,470,550,760,800,955,1260,1290)
bears.spring=c(280,445,520,715,750,890,990,900)
Weight=c(bears.winter, bears.spring)
Season=c(rep("Winter",8), rep("Spring",8))
bear.data=data.frame(Season,Weight)
bear.data

plot(bears.winter-bears.spring)
t.test(bears.winter, bears.spring, paired=TRUE)

wilcox.test(bears.winter, bears.spring, paired=TRUE)

# ---- 5.5 Rank alternatives to one-way ANOVA: Kruskal-Wallis and Rfit ----

# install.packages("Rfit") # if not already installed
library(Rfit)
head(quail, 10) # display the first 10 rows of the quail data
set

with(quail, boxplot(ldl ~ treat))
aggregate(ldl~treat, quail, median)

quail.lm.fit = with(quail, lm(ldl ~ treat))
anova(quail.lm.fit)

with(quail, kruskal.test(ldl, treat))

quail.oneway.rfit = with(quail, oneway.rfit(ldl, treat))
quail.oneway.rfit

summary(quail.oneway.rfit, method="tukey")

# ---- 5.6 Do the Wilcoxon tests and Kruskal-Wallis test compare the medians of groups? ----

Group.A=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
Group.B=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,
2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
median(Group.A) # 0

median(Group.B) # 0

sd(Group.A) # 0.5

sd(Group.B) # 1.0

wilcox.test(Group.A, Group.B)

# ---- 5.7 When should I use t-tests and ANOVA versus rank methods? ----

#==============================================================================
# 6 If you have non-normal data and more than one explanatory variable: Rank methods for multiple explanatory variables
#==============================================================================

library(tidyverse)
library(coin)
library(Rfit)
library(rankFD)

# ---- 6.1 Overview of rank methods for multiple explanatory variables ----

# ---- 6.2 Example data set with outliers ----

anova.2way.data.outlier = tibble::tribble(
~Treatment, ~Litter, ~Response,
"Drug", "Litter.1", 4,
"Drug", "Litter.1", 8,
"Drug", "Litter.2", 11,
"Drug", "Litter.2", 12,
"Drug", "Litter.3", 13,
"Drug", "Litter.3", 15,
"Drug", "Litter.4", 20,
"Drug", "Litter.4", 32,
"Placebo", "Litter.1", 10,
"Placebo", "Litter.1", 15,
"Placebo", "Litter.2", 16,
"Placebo", "Litter.2", 17,
"Placebo", "Litter.3", 19,
"Placebo", "Litter.3", 20,
"Placebo", "Litter.4", 21,
"Placebo", "Litter.4", 31
)
# van Elteren stratified Wilcoxon test from the coin R package
# requires that covariates are defined to be factors
anova.2way.data.outlier$Litter =
as.factor(anova.2way.data.outlier$Litter)
anova.2way.data.outlier$Treatment =
as.factor(anova.2way.data.outlier$Treatment)
str(anova.2way.data.outlier)

ggplot(anova.2way.data.outlier, aes(y=Response, x=Treatment))
+ geom_point() +
facet_wrap(~Litter) +
ggtitle("Response to treatment \nbroken out by Litter, with outliers")

# ---- 6.3 Analysis using standard anova ----

lm.2way.outlier = lm(Response ~ Litter + Treatment,
data= anova.2way.data.outlier)
anova(lm.2way.outlier)

# ---- 6.4 Rank-based analysis using the van Elteren stratified Wilcoxon test ----

# van Elteren stratified Wilcoxon test
wilcox_test(Response ~ Treatment | Litter,
data
= anova.2way.data.outlier,
distribution = "asymptotic")

# ---- 6.5 Rank-based analysis using the Rfit package ----

# install.packages("Rfit") # if not already installed
library(Rfit)
anova.2way.data.outlier.rfit = with(anova.2way.data.outlier,
rfit(Response ~ Litter + Treatment))
summary(anova.2way.data.outlier.rfit)

with(anova.2way.data.outlier, raov(Response ~ Litter *
Treatment))

anova.2way.data = tibble::tribble(
~Treatment, ~Litter, ~Response,
"Drug", "Litter.1", 4,
"Drug", "Litter.1", 8,
"Drug", "Litter.2", 11,
"Drug", "Litter.2", 12,
"Drug", "Litter.3", 13,
"Drug", "Litter.3", 15,
"Drug", "Litter.4", 20,
"Drug", "Litter.4", 21,
"Placebo", "Litter.1", 10,
"Placebo", "Litter.1", 15,
"Placebo", "Litter.2", 16,
"Placebo", "Litter.2", 17,
"Placebo", "Litter.3", 19,
"Placebo", "Litter.3", 20,
"Placebo", "Litter.4", 21,
"Placebo", "Litter.4", 23
)
summary(rfit(Response ~ Litter + Treatment, data =
anova.2way.data))

# ---- 6.6 Compare power for Wilcoxon, stratified Wilcoxon and rfit: data with outliers ----

# ---- 6.7 Compare power for Wilcoxon, stratified Wilcoxon and rfit: no outliers ----

# ---- 6.8 Limitations of the rank methods ----

# ---- 6.9 Appendix: How does Rfit rank regression work? ----

#==============================================================================
# 7 Survival Analysis
#==============================================================================

library(survival)
library(ggplot2)
library(tidyverse)
library(survminer)

# ---- 7.1 Avoid log rank tests. Use Cox proportional hazards regression to include additional variables that affect survival, which will increase power. ----

# mouse.lifespan.data
mouse.lifespan.data = data.frame(
treatment = c("CR", "CR", "CR", "CR", "CR", "CR", "CR",
"CR", "CR", "CR",
"Control", "Control", "Control", "Control",
"Control",
"Control", "Control", "Control", "Control",
"Control"),
sex = c("M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
"M", "M", "M", "M", "M", "F", "F", "F", "F", "F"),
months = c(10, 13, 18, 20, 22, 18, 19, 20, 24, 25,
9, 10, 14, 17, 18, 12, 15, 16, 21, 22),
event = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
mouse.lifespan.data

# Create life table estimates broken out by treatment
mouse.lifespan.survfit.by.treatment = survfit(Surv(months,
event == 1) ~ treatment, data = mouse.lifespan.data)
# Plot using ggsurvplot
ggsurvplot(
mouse.lifespan.survfit.by.treatment,
data = mouse.lifespan.data,
title = "Survival by treatment group",
legend.title = "Treatment",
legend = "bottom",
conf.int = FALSE,
pval = FALSE,
risk.table = FALSE,
legend.labs = c("Control", "Caloric restriction"),
palette = c("black", "red"),
linetype = c("solid", "dashed"),
size = 0.5,
ggtheme = theme_minimal()
)
# Create life table estimates broken out by sex
mouse.lifespan.survfit.by.sex = survfit(Surv(months, event ==
1) ~ sex,
data =
mouse.lifespan.data)
# Plot using ggsurvplot
ggsurvplot(
mouse.lifespan.survfit.by.sex,
data = mouse.lifespan.data,
title = "Survival by sex",
legend.title = "Sex",
legend = "bottom",
conf.int = FALSE,
pval = FALSE,
risk.table = FALSE,
legend.labs = c("Female", "Male"),
palette = c("black", "red"),
linetype = c("solid", "dashed"),
size = 0.5,
ggtheme = theme_minimal()
)
# Log rank test using survdiff
surv.diff.mouse.rx = survdiff(Surv(months, event == 1) ~
treatment,
data = mouse.lifespan.data)
surv.diff.mouse.rx

surv.diff.mouse.sex = survdiff(Surv(months, event == 1) ~ sex,
data = mouse.lifespan.data)
surv.diff.mouse.sex

# Cox PH for treatment + sex
surv.diff.mouse.ph = coxph(Surv(months, event == 1) ~
treatment + sex,
data = mouse.lifespan.data)
summary(surv.diff.mouse.ph)

# mouse.lifespan.litter.data
mouse.lifespan.litter.data = data.frame(
treatment = c("CR", "CR", "CR", "CR", "CR", "CR", "CR",
"CR", "CR", "CR",
"Control", "Control", "Control", "Control",
"Control",
"Control", "Control", "Control", "Control",
"Control"),
sex = c("M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
"M", "M", "M", "M", "M", "F", "F", "F", "F", "F"),
litter = c("A", "B", "C", "D", "E", "A", "B", "C", "D", "E",
"A", "B", "C", "D", "E", "A", "B", "C", "D",
"E"),
months = c(13, 16, 17, 23, 20, 21, 19, 22, 22, 24,
9, 13, 11, 14, 16, 19, 16, 17, 21, 23),
event = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
mouse.lifespan.litter.data

## KM plots
# Create life table estimates broken out by treatment
mouse.lifespan.litter.survfit.by.treatment = survfit(
Surv(months, event == 1) ~ treatment, data =
mouse.lifespan.litter.data)
# Plot using ggsurvplot
ggsurvplot(
mouse.lifespan.litter.survfit.by.treatment,
data = mouse.lifespan.litter.data,
title = "Survival by treatment group",
legend.title = "Treatment",
legend = "bottom",
conf.int = FALSE,
pval = FALSE,
risk.table = FALSE,
legend.labs = c("Control", "Caloric restriction"),
palette = c("black", "red"),
linetype = c("solid", "dashed"),
size = 0.5,
ggtheme = theme_minimal()
)
# Create life table estimates broken out by sex
mouse.lifespan.litter.survfit.by.sex = survfit(
Surv(months, event == 1) ~ sex, data =
mouse.lifespan.litter.data)
# Plot using ggsurvplot
ggsurvplot(
mouse.lifespan.litter.survfit.by.sex,
data = mouse.lifespan.litter.data,
title = "Survival by sex",
legend.title = "Sex",
legend = "bottom",
conf.int = FALSE,
pval = FALSE,
risk.table = FALSE,
legend.labs = c("Female", "Male"),
palette = c("black", "red"),
linetype = c("solid", "dashed"),
size = 0.5,
ggtheme = theme_minimal()
)
# Create life table estimates broken out by litter
mouse.lifespan.litter.survfit.by.litter = survfit(
Surv(months, event == 1) ~ litter, data =
mouse.lifespan.litter.data)
# Plot using ggsurvplot
ggsurvplot(
mouse.lifespan.litter.survfit.by.litter,
data = mouse.lifespan.litter.data,
title = "Survival by litter",
legend.title = "Litter",
legend = "bottom",
conf.int = FALSE,
pval = FALSE,
risk.table = FALSE,
legend.labs = c("A", "B", "C", "D", "E"),
palette = c("black", "red", "blue", "purple", "brown"),
linetype = "strata",
size = 0.5,
ggtheme = theme_minimal()
)
# Log rank test using survdiff
surv.diff.mouse.rx = survdiff(Surv(months, event == 1) ~
treatment,
data =
mouse.lifespan.litter.data)
surv.diff.mouse.rx # p = 0.08

surv.diff.mouse.sex = survdiff(Surv(months, event == 1) ~ sex,
data =
mouse.lifespan.litter.data)
surv.diff.mouse.sex # p = 0.02

surv.diff.mouse.litter = survdiff(Surv(months, event == 1) ~
litter,
data =
mouse.lifespan.litter.data)
surv.diff.mouse.litter # p = 0.1

# Cox PH for treatment + sex
surv.diff.mouse.ph = coxph(Surv(months, event == 1) ~
treatment + sex,
data = mouse.lifespan.litter.data)
summary(surv.diff.mouse.ph)

# treatment + sex + litter
surv.diff.mouse.ph = coxph(Surv(months, event == 1) ~
treatment + sex + litter,
data = mouse.lifespan.litter.data)
summary(surv.diff.mouse.ph)

# ---- 7.2 Omitting important predictors shrinks hazard ratios towards the null hypothesis value and reduces power. ----

# Generate simulated data
set.seed(1)
n = 10000
treatment = c(rep(0, n/2), rep(1, n/2))
age = runif(n, min = 3, max = 18)
# Important continuous
predictor
beta1 = 1.0
beta2 = 0.2
# Hazard: h(t) = h0(t) * exp(beta1*treatment + beta2*age)
hazard = exp(beta1 * treatment + beta2 * age)
time = rexp(n, rate = hazard)
event = rep(1, n)
# No censoring
omit.predictor.data = data.frame(time, event, treatment, age)
head(omit.predictor.data)

### KM plots
# Create life table estimates broken out by treatment
omit.predictor.survfit.by.treatment = survfit(Surv(time, event
== 1) ~ treatment, data = omit.predictor.data)
# Plot using ggsurvplot
ggsurvplot(
omit.predictor.survfit.by.treatment,
data = omit.predictor.data,
title = "Survival by treatment",
legend.title = "treatment",
legend = "bottom",
legend.labs = c(0, 1),
palette = c("black", "red"),
linetype = c("solid", "dashed"),
linewidth = 0.5,
ggtheme = theme_minimal()
)
# Create age categories to use in the KM plot
omit.predictor.data = omit.predictor.data %>%
mutate(age.categorical = case_when(
age >= 3
& age < 8
~ "3 to 8",
age >= 8
& age < 13 ~ "8 to 13",
age >= 13 & age < 18 ~ "GT13"))
# Create life table estimates broken out by age category
omit.predictor.survfit.by.age = survfit(Surv(time, event == 1)
~ age.categorical, data = omit.predictor.data)
# Plot using ggsurvplot
ggsurvplot(
omit.predictor.survfit.by.age,
data = omit.predictor.data,
title = "Survival by age",
legend.title = "age",
legend = "bottom",
legend.labs = c("3 to 8", "8 to 13", "GT13"),
palette = c("red", "blue", "black"),
linetype = "strata",
size = 0.5,
ggtheme = theme_minimal()
)
# Fit Full Model (Correct)
model_full = coxph(Surv(time, event) ~ treatment + age, data =
omit.predictor.data)
summary(model_full)$coefficients["treatment", ]

model_reduced = coxph(Surv(time, event) ~ treatment, data =
omit.predictor.data)
round(summary(model_reduced)$coefficients["treatment", ],
digits = 4)

# Create sample of 60 observations from omit.predictor.data
set.seed(123)
omit.predictor.data.60 =
omit.predictor.data[sample(nrow(omit.predictor.data), size =
60),]
# Fit Full Model (Correct)
model_full_60 = coxph(Surv(time, event) ~ treatment + age,
data = omit.predictor.data.60)
summary(model_full_60)$coefficients["treatment", ]

model_reduced_60 = coxph(Surv(time, event) ~ treatment, data =
omit.predictor.data.60)
round(summary(model_reduced_60)$coefficients["treatment", ],
digits = 4)

# ---- 7.3 Balancing covariates across treatment arms is not sufficient. Include covariates in the analysis. ----

# balanced.data
balanced.data = data.frame(
treatment = c("Drug", "Drug", "Drug", "Drug", "Drug",
"Drug", "Drug", "Drug", "Drug", "Drug",
"Control", "Control", "Control", "Control",
"Control",
"Control", "Control", "Control", "Control",
"Control"),
age = c("40", "50", "60", "70", "80", "40", "50", "60",
"70", "80",
"40", "50", "60", "70", "80", "40", "50", "60",
"70", "80"),
months = c(10, 13, 18, 20, 22, 18, 19, 20, 24, 25,
9, 10, 14, 17, 18, 12, 15, 16, 21, 22),
event = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
balanced.data

with(balanced.data, table(age, treatment))

## KM plot by treatment
# Create life table estimates broken out by treatment
balanced.survfit.by.treatment = survfit(Surv(months, event ==
1) ~ treatment, data = balanced.data)
# KM plot using ggsurvplot
ggsurvplot(
balanced.survfit.by.treatment,
data = balanced.data,
title = "Survival by treatment group",
legend.title = "Treatment",
legend = "bottom",
legend.labs = c("Control", "Drug"),
palette = c("black", "red"),
linetype = c("solid", "dashed"),
size = 0.5,
ggtheme = theme_minimal()
)
# Cox PH with only treatment
surv.diff.balanced.ph.treatment = coxph(Surv(months, event ==
1) ~ treatment,
data = balanced.data)
summary(surv.diff.balanced.ph.treatment)

# Cox PH with treatment + age
surv.diff.balanced.ph = coxph(Surv(months, event == 1) ~
treatment + age, data = balanced.data)
summary(surv.diff.balanced.ph)

# ---- 7.4 Avoid testing for a difference in percent surviving at a specific timepoint. Keep the time information and use survival analysis to get better power. ----

n_sim
= 100
# Number of Monte Carlo replications
n_per_arm
= 300
# patients per treatment arm
median_A
= 30
# median survival in treatment group A
(months)
median_B
= 40
# median survival in treatment group B
(months)
max_follow = 60
# administrative censoring (months)
# logistic regression cutpoints:
cutoffs
= c(12, 24, 36, 48, 60)
alpha
= 0.05

#==============================================================================
# 8 For count outcomes consider negative binomial regression and Wilcoxon tests
#==============================================================================

library(MASS) # use the glm.nb function from the MASS library
library(tidyverse)

# ---- 8.1 Example: Do mice with a genetic variant have smaller litters? ----

## Create data set
n.0 = 10
n.1 = 10
litter.sizes.data = data.frame(Genotype = c(rep("Wild type",
n.0), rep("Variant", n.1)),
litter.size = c(1,2,4,4,6,7,7,8,11,12, 0,2,2,2,3,4,4,5,6,8))
# Graph the litter size by genotype
ggplot(litter.sizes.data, aes(y=litter.size,
x=as.factor(Genotype), fill=Genotype)) +
geom_dotplot(binaxis='y', stackdir='center',
stackratio=1.5, dotsize=1.5) + theme_minimal()
+
theme(legend.position="none") + ylab("Litter size") +
xlab("Genotype") + ylim(0,15) +
ggtitle("Litter size by Genotype")
# Mean by Genotype
with(litter.sizes.data, aggregate(litter.size ~ Genotype,
FUN="mean"))

t.test(litter.size ~ Genotype, data=litter.sizes.data)

wilcox.test(litter.size ~ Genotype, data=litter.sizes.data)

# negative binomial model
litter.sizes.nb = glm.nb(litter.size ~ Genotype,
data=litter.sizes.data)
summary(litter.sizes.nb)

# ---- 8.2 Why not use t-tests or ordinary linear regression for count data? ----

# ---- 8.3 Compare power for data from a negative binomial distribution ----

# ---- 8.4 Compare power when data are from a mixture of distributions ----

#==============================================================================
# 9 Avoid turning quantitative variables into categorical variables
#==============================================================================

# ---- 9.1 Example. Creating a categorical response with 2 values ----

# SBP data
set.seed(11111)
drug
= round(rnorm(n=20, mean= -10, sd=20), 0)
placebo= round(rnorm(n=20, mean=
10, sd=20), 0)
treatment
= c(rep("Drug",
length(drug)),
rep("Placebo", length(placebo)))
SBP.change = c(drug, placebo)
SBP.data
= data.frame(treatment, SBP.change)
SBP.data = within(SBP.data, {
SBP.change.category = NA
SBP.change.category[SBP.change <
0] = "Improved"
SBP.change.category[SBP.change >= 0] = "Not Improved"
})
SBP.change.category.table = with(SBP.data,
table(SBP.change.category,
treatment))
head(SBP.data)

with(SBP.data, boxplot(SBP.change ~ treatment))
with(SBP.data, t.test(SBP.change ~ treatment))

SBP.change.category.table

fisher.test(SBP.change.category.table)

power.t.test(delta = 20, sd=20, power=0.8)

# proportion of drug group below 0
pnorm(q=0, mean=-10, sd=20) # 0.69

# proportion of placebo group below 0
pnorm(q=0, mean=10, sd=20) # 0.31

power.prop.test(p1=0.69, p2=0.31, power = 0.8)

# ---- 9.2 If you must create a categorical variable, create more than two categories. ----

# Cholesterol data
drug.mean
= 190
placebo.mean
= 210
sd.cholesterol
= 50
n.per.group
= 100
set.seed(3)
drug
= round(rnorm(n=n.per.group, mean=drug.mean,
sd=sd.cholesterol), 0)
placebo = round(rnorm(n=n.per.group, mean=placebo.mean,
sd=sd.cholesterol), 0)
treatment
= c(rep("Drug",
length(drug)),
rep("Placebo", length(placebo)))
cholesterol = c(drug, placebo)
cholesterol.data = data.frame(treatment, cholesterol)
cholesterol.data = within(cholesterol.data, {
chol.LT200 = NA
chol.LT200[cholesterol <
200] = "LT200"
chol.LT200[cholesterol >= 200] = "GTE200"
})
cholesterol.data = within(cholesterol.data, {
chol.category = NA
chol.category[cholesterol <
200]
=
"LT200"
chol.category[cholesterol >= 200 & cholesterol < 240]=
"GTE200LT240"
chol.category[cholesterol >= 240]
=
"GTE240"
})
LT.200.table = with(cholesterol.data,
table(chol.LT200, treatment))
cholesterol.ordered.table = with(cholesterol.data,
table(treatment,
chol.category))
head(cholesterol.data, 10)

with(cholesterol.data, boxplot(cholesterol ~ treatment))
with(cholesterol.data, t.test(cholesterol ~ treatment))

LT.200.table

fisher.test(LT.200.table)

cholesterol.ordered.table

# install.packages("vcdExtra")
library(vcdExtra)
CMHtest(cholesterol.ordered.table)

# ---- 9.3 Logistic regression and ordinal logistic regression for categorical response variables. ----

# ---- 9.4 When can it be helpful to use categories instead of quantitative measures? ----

#==============================================================================
# 10 Avoid Chi-Squared tests when you have ordered categories: use CMH tests instead
#==============================================================================

library(vcdExtra)

# ---- 10.1 Example: CMH test for table with ordered columns ----

ordered.columns.table = array(
data = c(1, 5, 3, 3, 5, 1),
dim = c(2, 3),
dimnames = list(
Treatment = c("drug", "placebo"),
Result = c("no remission", "partial remission", "complete remission")
)
)
ordered.columns.table

chisq.test(ordered.columns.table)

fisher.test(ordered.columns.table)

CMHtest(ordered.columns.table)

# ---- 10.2 Example: CMH test for table with ordered rows and columns ----

ordered.rows.and.columns.table = array(
data = c(5, 2, 1, 3, 3, 3, 2, 4, 6),
dim = c(3, 3),
dimnames = list(
Treatment = c("placebo", "low dose", "high dose"),
Result = c("no remission", "partial remission", "complete remission")
)
)
ordered.rows.and.columns.table

chisq.test(ordered.rows.and.columns.table)

fisher.test(ordered.rows.and.columns.table)

CMHtest(ordered.rows.and.columns.table)

# ---- 10.3 Watch out for which CMH test your software is performing. ----

#==============================================================================
# 11 Balancing variables across treatment groups is not sufficient to get best p-values
#==============================================================================

# ---- 11.1 Short example 1: Testing a treatment when the response is also affected by patient's age. ----

# Load required libraries
library(ggplot2)
library(car)
# Generate a data set "protein.data"
# Serum protein vs Age and Treatment
protein.data = tibble::tribble(
~Subject, ~Age, ~Treatment,
~Serum.protein,
1,
40, "Placebo",
63,
2,
40, "Drug",
57,
3,
50, "Placebo",
75,
4,
50, "Drug",
65,
5,
60, "Placebo",
82,
6,
60, "Drug",
78,
7,
70, "Placebo",
94,
8,
70, "Drug",
86
)
graph.1 = ggplot(protein.data, aes(y=Serum.protein, x=Age,
color=Treatment)) +
geom_point(size=3) +
ggtitle("Serum protein vs Age and Treatment") +
ylim(50, 100)
plot(graph.1)
t.test(Serum.protein ~ Treatment, data=protein.data,
var.equal=TRUE)

protein.lm = lm(Serum.protein ~ Treatment + Age,
data=protein.data)
summary(protein.lm)

# ---- 11.2 Short example 2: Adjusting for 2 or more variables that affect the response. ----

set.seed(100) # set.seed makes rnorm reproducible.
Treatment = c(rep(0, 8), rep(1, 8))
Instrument = c(rep(c(0, 1),8))
Operator = c(rep(c(0, 0, 1, 1),4))
Response =
20 + Treatment + Instrument + Operator +
rnorm(n=16)
experiment.2.data = data.frame(Treatment = Treatment,
Instrument =Instrument,
Operator = Operator , Response=
Response)
experiment.2.data

with(experiment.2.data, t.test(Response ~ Treatment,
var.equal=TRUE))

with(experiment.2.data,
summary(lm(Response ~ Treatment + Instrument + Operator)))

# ---- 11.3 Detailed presentation for example 1 ----

protein.data

library(car)
ggplot(protein.data, aes(y=Serum.protein, x=Treatment)) +
geom_point() +
ggtitle("Serum protein vs Treatment") + ylim(50, 100) +
geom_point(aes(color=Treatment))
t.test(Serum.protein ~ Treatment, data=protein.data,
var.equal=TRUE)

graph.1
Anova(lm(Serum.protein ~ Treatment + Age, data =
protein.data))

# ---- 11.4 Detailed presentation for example 2: How to adjust for two or more variables. ----

experiment.2.data
Treatment Instrument Operator Response
1
0
0
0 19.49781
2
0
1
0 21.13153
3
0
0
1 20.92108
4
0
1
1 22.88678
5
0
0
0 20.11697
6
0
1
0 21.31863
7
0
0
1 20.41821
8
0
1
1 22.71453
9
1
0
0 20.17474
10
1
1
0 21.64014
11
1
0
1 22.08989
12
1
1
1 23.09627
13
1
0
0 20.79837
14
1
1
0 22.73984
15
1
0
1 22.12338
16
1
1
1 22.97068
with(experiment.2.data, plot(Response ~ Treatment))
with(experiment.2.data, plot(Response ~ Instrument))
with(experiment.2.data, plot(Response ~ Operator))
with(experiment.2.data, t.test(Response ~ Treatment,
var.equal=TRUE))

with(experiment.2.data,
summary(lm(Response ~ Treatment + Instrument + Operator)))

#==============================================================================
# 12 Avoid correlated explanatory variables in multiple regression and ANOVA
#==============================================================================

library(tidyverse)

# ---- 12.1 Example: How correlated explanatory variables affect linear regression ----

response =
c(25,22,18,16,19,18,22,17,18,15,20,23,19,26,22,24,25,21,27,24,27,
19,22,26,27,29,25,22,18,19,27,25,23,18,19,17)
weight =
c(9,10,12,19,20,24,7,11,12,14,17,22,9,11,12,15,18,25,9,10,13,13,17,
20,8,10,11,16,19,25,10,10,12,15,21,24)
dose=c(rep(0,6), rep(1,6), rep(2,6), rep(3,6), rep(4,6),
rep(5,6))
age = rep(3:8, 6)
mouse.data = data.frame(response, dose, age, weight)
head(mouse.data)

hist(mouse.data$response)
pairs(mouse.data)
round(cor(mouse.data), 3)

lm.3vars= lm(response ~ dose +
age + weight, data=
mouse.data)
summary(lm.3vars)

lm.dose = lm(response ~ dose, data= mouse.data)
summary(lm.dose)

lm.dose.age = lm(response ~ dose + age, data = mouse.data)
summary(lm.dose.age)

lm.dose.weight = lm(response ~ dose + weight, data=
mouse.data)
summary(lm.dose.weight)

# ---- 12.2 What have we learned, and what should we do in future experiments? ----

# ---- 12.3 How to deal with correlated variables ----

# ---- 12.4 Final example: the OncotypeDx breast cancer recurrence test ----

#==============================================================================
# 13 Look (out) for interactions: when subgroups respond differently to treatment
#==============================================================================

# ---- 13.1 Example: how interactions can make p-values worse ----

CR.data = data.frame(CR = c("Yes", "Yes", "Yes", "Yes", "Yes",
"Yes", "Yes", "Yes", "No", "No", "No", "No", "No", "No", "No",
"No"),
Learning.score = c(65, 68, 72, 60, 69, 71, 62, 70, 64, 67, 62,
61, 68, 63, 64, 65))
CR.data

library(ggplot2)
mouse.CR.graph = ggplot(CR.data, aes(y=Learning.score, x=CR))
+
geom_point() + ggtitle("Mouse CR data") +
scale_y_continuous(breaks=0:80*10, limits=c(50, 80))
plot(mouse.CR.graph)
with(CR.data, t.test(Learning.score ~ CR, var.equal = TRUE))

# t-test using lm
mouse.lm.CR=with(CR.data, lm(Learning.score ~ CR))
summary(mouse.lm.CR)

CR.Age.data = data.frame(CR = c("Yes", "Yes", "Yes", "Yes",
"Yes", "Yes", "Yes", "Yes","No", "No", "No", "No", "No", "No",
"No", "No"),
Learning.score = c(65, 68, 72, 60, 69, 71, 62, 70, 64, 67, 62,
61, 68, 63, 64, 65),
Age = c("Young", "Old", "Old", "Young", "Old", "Old", "Young",
"Old", "Young", "Young", "Old", "Young", "Old", "Old", "Old",
"Old"))
CR.Age.data

mouse.lm.CR.Age=with(CR.Age.data, lm(Learning.score ~ CR +
Age))
summary(mouse.lm.CR.Age)

mouse.interaction.graph = ggplot(CR.Age.data,
aes(y=Learning.score, x=CR,
group=Age, col=Age,
shape=Age)) +
geom_point() +
ggtitle("Interaction of CR with Age \nin their effect on Learning score") +
scale_y_continuous(breaks=0:80*10, limits=c(50, 80))
plot(mouse.interaction.graph)
with(CR.Age.data, interaction.plot(CR, Age,
Learning.score,
ylim = c(60, 75)))
mouse.lm.CR.Age.interaction = with(CR.Age.data,
lm(Learning.score ~ CR * Age))
summary(mouse.lm.CR.Age.interaction)

young.mouse.data = subset(CR.Age.data, Age=="Young")
young.mouse.CR.graph = ggplot(young.mouse.data,
aes(y=Learning.score, x=CR)) +
geom_point() +
ggtitle("Young mouse CR data") +
scale_y_continuous(breaks=0:80*10, limits=c(50, 80))
plot(young.mouse.CR.graph)
young.mouse.lm
= lm(Learning.score ~ CR,
data=young.mouse.data)
summary(young.mouse.lm)

old.mouse.data
= subset(CR.Age.data, Age=="Old")
old.mouse.CR.graph = ggplot(old.mouse.data,
aes(y=Learning.score, x=CR)) +
geom_point() +
ggtitle("Old mouse CR data") +
scale_y_continuous(breaks=0:80*10, limits=c(50, 80))
plot(old.mouse.CR.graph)
old.mouse.lm
= lm(Learning.score ~ CR,
data=old.mouse.data)
summary(old.mouse.lm)

# ---- 13.2 What have we seen and learned? ----

# ---- 13.3 Multiple hypothesis testing and interactions ----

choose(10,2) # Gives 45.

choose(10,3) # Gives 120.

# ---- 13.4 When should you look for interactions? ----

#==============================================================================
# 14 Use two primary analyses each with alpha 0.025 instead of one with 0.05
#==============================================================================

# ---- 14.1 Power for a test of a single gene using alpha of 0.05 ----

power.t.test(delta=1, power = 0.8, sd = 2, sig.level
= 0.05)

# ---- 14.2 Power for tests of two genes using alpha of 0.025 each ----

power.t.test(n = 64, delta=1, sd = 2, sig.level
= 0.025)

# ---- 14.3 Bonferroni and Holm's adjustments for multiple hypothesis testing ----

#==============================================================================
# 15 Avoid one-variable-at-a-time experiments. Use factorial experiment design.
#==============================================================================

# ---- 15.1 An intuitive view of experiment design: baking cookies ----

# ---- 15.2 The problems with "One variable at a time" ----

# ---- 15.3 Are interactions important in medicine and in biological research? ----

# ---- 15.4 Factorial experiment design: an alternative to OVAT ----

# ---- 15.5 Example: experiment design for 3 factors ----

# ---- 15.6 Factorial designs have more power than OVAT designs ----

# ---- 15.7 What about replication in OVAT versus factorial designs? ----

# ---- 15.8 Statistical analysis of factorial designs ----

temperature = factor(rep(c("Low","Low","High","High"), 4))
salt.conc
= factor(rep(c("Low","High","Low","High"), 4))
buffer.type = factor(rep(c("A","A","A","A","B","B","B","B"),
2))
yield
= c(7.88, 6.54, 16.12, 11.14, 9.26, 10.43, 13.92,
8.47,
7.63, 6.11, 15.45, 11.72, 9.80,
7.22, 11.89,
14.57)
yield.data
= data.frame(temperature, salt.conc, buffer.type,
yield)
yield.data

# Graph the results using plot.design.
plot.design(yield ~ buffer.type + temperature + salt.conc,
data = yield.data, ylim = c(7, 14))

yield.data.lm = lm(yield ~ buffer.type + temperature +
salt.conc, data = yield.data)
summary(yield.data.lm)

# ---- 15.9 When may OVAT be preferable? ----

# ---- 15.10 Recommended YouTube videos ----

# ---- 15.11 Recommended free textbook PDF ----

#==============================================================================
# 16 Rescue failed experiments using fractional factorial designs
#==============================================================================

# ---- 16.1 Example: How to bake good cookies: A fractional factorial experiment ----

cookie.data = data.frame(flour
=
c("brown","white","brown","white","brown","white","brown","white"),
temp
=
c("low","low","high","high","low","low","high","high"),
soda
=
c("low","low","low","low","high","high","high","high"),
time
=
c("long","short","short","long","long","short","short","long"),
water
=
c("high","low","high","low","low","high","low","high"),
sift
= c("sift","sift","no","no","no","no","sift","sift"),
butter =
c("low","high","high","low","high","low","low","high"),
yield
= c(47, 74, 84, 62, 53, 78, 87, 60))
cookie.data[1:7] = lapply(cookie.data[1:7], factor)
cookie.data

plot.design(cookie.data)

with(cookie.data, aggregate(yield ~ time, FUN = mean))

lm.cookie = lm(yield ~ ., data = cookie.data)
summary(lm.cookie)

lm.cookie.2vars = lm(yield ~ time + temp, data = cookie.data)
summary(lm.cookie.2vars)

lm.cookie.2vars.interaction = lm(yield ~ time * temp, data =
cookie.data)
summary(lm.cookie.2vars.interaction)

# ---- 16.2 Additional resources ----

#==============================================================================
# 17 The wrong way and the right way to get p < 0.05: multiple hypothesis testing
#==============================================================================

# ---- 17.1 The probability of a false positive result with multiple comparisons ----

# ---- 17.2 Methods to adjust p-values to control for multiple hypothesis testing ----

# ---- 17.3 Hierarchical hypothesis testing ----

# ---- 17.4 False discovery rate: an alternative to controlling false positives ----

#==============================================================================
# 18 Is P < 0.05 the right goal?
#==============================================================================

#==============================================================================
# 19 What if you still get p > 0.05?
#==============================================================================

# ---- 19.1 A possibly false negative result ----

# ---- 19.2 A true negative result ----

#==============================================================================
# 20 Don't know R? Use Claude to run these analyses for you
#==============================================================================

# ---- 20.1 A short protocol for getting good results ----

# ---- 20.2 A complete example, start to finish ----

# ---- 20.3 Example queries for the methods in this book ----

# ---- 20.4 Asking for model checking and diagnostics ----

# ---- 20.5 Asking for alternative analyses and sensitivity checks ----

# ---- 20.6 Asking for a plain-English interpretation ----

# ---- 20.7 An AI assistant is a collaborator, not an oracle ----

#==============================================================================
# References
#==============================================================================

# Load library
library(tidyverse)
library(ggthemes)
library(stargazer)

# Load data
df <- read_delim("input/Empirical_Problems_Data.dat", delim = " ", col_names = FALSE, trim_ws = TRUE) %>%
  set_names("wage", "female", "black", "educ", "exp", "married")

mypal <- palette_pander(5)

# Q1
## (a)
p.1.a.1 <- ggplot(df) +
  geom_density(aes(x = wage), fill = mypal[1]) +
  theme_few(base_size = 12, base_family = "sans")

p.1.a.2 <- ggplot(df) +
  geom_density(aes(x = log(wage)), fill = mypal[2]) +
  theme_few(base_size = 12, base_family = "sans")
  
ggsave(filename = "p.1.a.1.pdf", p.1.a.1, path = "output/figs")
ggsave(filename = "p.1.a.2.pdf", p.1.a.2, path = "output/figs")

## (b)
m.1.b.1 <- lm(wage ~ educ + exp + I(exp^2), df)
resid1 <- m.1.b.1$residuals
m.1.b.2 <- lm(log(wage) ~ educ + exp + I(exp^2), df)
resid2 <- m.1.b.2$residuals

p.1.b.1 <- ggplot(df) +
  geom_density(aes(x = resid1), fill = mypal[1]) +
  theme_few(base_size = 12, base_family = "sans")

p.1.b.2 <- ggplot(df) +
  geom_density(aes(x = resid2), fill = mypal[2]) +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "p.1.b.1.pdf", p.1.b.1, path = "output/figs")
ggsave(filename = "p.1.b.2.pdf", p.1.b.2, path = "output/figs")

## (c)

## (d)
m.1.d <- lm(log(wage) ~ black, df)
stargazer(m.1.d, type = "text", out = "output/tabs/m.1.d.tex", float = FALSE, single.row = TRUE)

## (e)
m.1.e <- lm(log(wage) ~ black + educ + poly(exp, 2), df)
stargazer(m.1.e, type = "text", out = "output/tabs/m.1.e.tex", float = FALSE, single.row = TRUE)

# Q2
## (a)
m.2.a.1 <- lm(log(wage) ~ female + educ + poly(exp, 2), df)
stargazer(m.2.a.1, type = "text", out = "output/tabs/m.2.a.1.tex", float = FALSE, single.row = TRUE)

m.2.a.2 <- lm(log(wage) ~ female + I(female * black) + educ + poly(exp, 2), df)
stargazer(m.2.a.2, type = "text", out = "output/tabs/m.2.a.2.tex", float = FALSE, single.row = TRUE)

## (b)
m.2.b.1 <- lm(log(wage) ~ married * female + educ + poly(exp , 2), df)
stargazer(m.2.b.1, type = "text", out = "output/tabs/m.2.b.1.tex", float = FALSE, single.row = TRUE)

resid1 <- m.2.b.1$residuals

m.2.b.2 <- lm(log(wage) ~ educ + poly(exp , 2), df)
resid2 <- m.2.b.2$residuals

rss1 <- sum(resid1^2)
rss2 <- sum(resid2^2)
df1 <- m.2.b.2$df.residual - m.2.b.1$df.residual
df2 <- m.2.b.1$df.residual

ftest <- ((rss2 - rss1)/df1)/(rss1/df2)
pvalue <- pf(ftest1, df1, df2, lower.tail = FALSE)

print(c(ftest, pvalue))

# Q3

# Q4

# Q5
## (c)

calculate_fittedvalue <- function(sex, biochem, physiol, genetic, pediatr, medicin, cert, clin, prate, exper, assistn, associa) {
  0.0638 * sex + -0.8648 * biochem + -1.0326 * physiol + -0.7110 * pediatr + -0.3660 * medicin + 0.1854 * cert + 0.1591 * clin +
    -0.0238 * prate + 0.0292 * exper + -0.000335 * exper^2 + -0.1908 * assistn + -0.0811 * associa + 12.0696
  }

# female
female <- calculate_fittedvalue(sex = 0, biochem = 0, physiol = 0, genetic = 0, pediatr = 0,medicin = 0,
                      cert = 1, clin = 1, prate = 0, exper = 5, assistn = 1, associa = 0)

# male
male <- calculate_fittedvalue(sex = 1, biochem = 0, physiol = 0, genetic = 0, pediatr = 0, medicin = 0, 
                      cert = 1, clin = 1, prate = 0, exper = 5, assistn = 1, associa = 0)

exp(male) - exp(female)

## (d)
df1 <- 2
df2 <- 261 - 14

ftest = (df2/df1) * (0.9352 - 0.9237)/(1 - 0.9352)
pvalue <- pf(ftest, df1, df2, lower.tail = FALSE)

print(c(ftest, pvalue))


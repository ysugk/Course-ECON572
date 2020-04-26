## Data import
mat <- R.matlab::readMat("input/Part IV/data.mat")
data <- mat$data
x <- data[ ,1] # (nX1) vector (age)
y <- data[ ,6] # (nX1) vector (happiness)
w <- data[ ,7] # (nX1) vector (income)

w_max <- quantile(w, 0.95)
w_min <- quantile(w, 0.05)
dummy <- (w <= w_max) & (w >= w_min)

y <- y[dummy]
x <- x[dummy]
w <- w[dummy]
w <- log(w)

## 1. (a) histogram for x and y
library(ggplot2)
library(ggthemes)
mypal <- palette_pander(8)

p1 <- ggplot() +
  geom_histogram(aes(x, after_stat(count / sum(count))), fill = mypal[1], color = "black", binwidth = 3) +
  labs(y = "frequency", x = "x", title = "Historgram of x") +
  theme_few(base_size = 12, base_family = "sans")

p2 <- ggplot() +
  geom_histogram(aes(y, after_stat(count / sum(count))), fill = mypal[2], color = "black", binwidth = 0.2) +
  labs(y = "frequency", x = "y", title = "Historgram of y") +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.1.a.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.a.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)

## 1. (b) naive density estimator for x and y
N <- 100 # number of estimation points in the domain of x and y
x_0 <- seq(min(x), max(x), length.out = N)
y_0 <- seq(min(y), max(y), length.out = N)

uniform <- function(u) abs(u) <= 1/2;
f_x_uniform <- vector("numeric", N)
f_y_uniform <- vector("numeric", N)
h_x <- sd(x) * (length(x)^(-1/5))
h_y <- sd(y) * (length(y)^(-1/5))

for (k in 1:N) {
 z_x <- uniform((x - x_0[k])/h_x)
 f_x_uniform[k] <- (1/(length(x) * h_x)) * sum(z_x)
 
 z_y <- uniform((y - y_0[k])/h_y)
 f_y_uniform[k] <- (1/(length(y) * h_y)) * sum(z_y)
}

p1 <- ggplot() +
  geom_line(aes(x_0, f_x_uniform)) +
  labs(x = "x", y = "f(x)", title = "Naive density estimation for x") +
  theme_few(base_size = 12, base_family = "sans")

p2 <- ggplot() +
  geom_line(aes(y_0, f_y_uniform)) +
  labs(x = "y", y = "f(y)", title = "Naive density estimation for y") +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.1.b.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.b.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)

## 1. (c) Kernel density estimator for x and y
gauss <- function(u) exp((-1/2) * u^2) / sqrt(2*pi);
epanechnikov <- function(u) (3/4) * (1 - u^2) * (abs(u) <= 1);
parzen <- function(u) (1 - 6*u^2 + 6*abs(u)^3) * (abs(u) <= 1/2) + (2*(1 - abs(u))^3) * (1/2 < abs(u)) * (abs(u) <= 1);
bartlett <- function(u) (1 - abs(u)) * (abs(u) <= 1);

f_x_gauss <- vector("numeric", N)
f_y_gauss <- vector("numeric", N)

f_x_epanechnikov <- vector("numeric", N)
f_y_epanechnikov <- vector("numeric", N)

f_x_parzen <- vector("numeric", N)
f_y_parzen <- vector("numeric", N)

f_x_bartlett <- vector("numeric", N)
f_y_bartlett <- vector("numeric", N)

for (k in 1:N) {
  z_x <- gauss((x - x_0[k])/h_x);
  f_x_gauss[k] <- (1/(length(x) * h_x)) * sum(z_x);
  
  z_x <- epanechnikov((x - x_0[k])/h_x);
  f_x_epanechnikov[k] <- (1/(length(x) * h_x)) * sum(z_x);
  
  z_x <- parzen((x - x_0[k])/h_x);
  f_x_parzen[k] <- (1/(length(x) * h_x)) * sum(z_x);
  
  z_x <- bartlett((x - x_0[k])/h_x);
  f_x_bartlett[k] <- (1/(length(x) * h_x)) * sum(z_x);

  z_y <- gauss((y - y_0[k])/h_y);
  f_y_gauss[k] <- (1/(length(y) * h_y)) * sum(z_y);
  
  z_y <- epanechnikov((y - y_0[k])/h_y);
  f_y_epanechnikov[k] <- (1/(length(y) * h_y)) * sum(z_y);
  
  z_y <- parzen((y - y_0[k])/h_y);
  f_y_parzen[k] <- (1/(length(y) * h_y)) * sum(z_y);
  
  z_y <- bartlett((y - y_0[k])/h_y);
  f_y_bartlett[k] <- (1/(length(y) * h_y)) * sum(z_y);
}

# gaussian
p1 <- ggplot() +
  geom_line(aes(x_0, f_x_gauss)) +
  labs(x = "x", y = "f(x)", title = "Gaussian density estimation for x") +
  theme_few(base_size = 12, base_family = "sans")

p2 <- ggplot() +
  geom_line(aes(y_0, f_y_gauss)) +
  labs(x = "y", y = "f(y)", title = "Gaussian density estimation for y") +
  theme_few(base_size = 12, base_family = "sans")

# epanechnikov
p3 <- ggplot() +
  geom_line(aes(x_0, f_x_epanechnikov)) +
  labs(x = "x", y = "f(x)", title = "Epanechnikov density estimation for x") +
  theme_few(base_size = 12, base_family = "sans")

p4 <- ggplot() +
  geom_line(aes(y_0, f_y_epanechnikov)) +
  labs(x = "y", y = "f(y)", title = "Epanechnikov density estimation for y") +
  theme_few(base_size = 12, base_family = "sans")

# parzen
p5 <- ggplot() +
  geom_line(aes(x_0, f_x_parzen)) +
  labs(x = "x", y = "f(x)", title = "Parzen density estimation for x") +
  theme_few(base_size = 12, base_family = "sans")

p6 <- ggplot() +
  geom_line(aes(y_0, f_y_parzen)) +
  labs(x = "y", y = "f(y)", title = "Parzen density estimation for y") +
  theme_few(base_size = 12, base_family = "sans")

# bartlett
p7 <- ggplot() +
  geom_line(aes(x_0, f_x_bartlett)) +
  labs(x = "x", y = "f(x)", title = "Bartlett density estimation for x") +
  theme_few(base_size = 12, base_family = "sans")

p8 <- ggplot() +
  geom_line(aes(y_0, f_y_bartlett)) +
  labs(x = "y", y = "f(y)", title = "Bartlett density estimation for y") +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.1.c.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.c.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.c.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.c.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.c.5.pdf", plot = p5, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.c.6.pdf", plot = p6, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.c.7.pdf", plot = p7, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.1.c.8.pdf", plot = p8, path = "output/figs/Part IV/", width = 8, height = 8)

## 2. (a) & (b)
m_gauss <- vector(mode = "numeric", N)
m_epanechnikov <- vector(mode = "numeric", N)
m_parzen <- vector(mode = "numeric", N)
m_bartlett <- vector(mode = "numeric", N)

for (k in 1:N) {
  z_x <- gauss((x - x_0[k])/h_x);
  m_gauss[k] <- sum(z_x * y) / sum(z_x);
  
  z_x <- epanechnikov((x - x_0[k])/h_x);
  m_epanechnikov[k] <- sum(z_x * y) / sum(z_x);
  
  z_x <- parzen((x - x_0[k])/h_x);
  m_parzen[k] <- sum(z_x * y) / sum(z_x);
  
  z_x <- bartlett((x - x_0[k])/h_x);
  m_bartlett[k] <- sum(z_x * y) / sum(z_x);
}

X <- cbind(rep(1, length(x)), x)
b <- solve(t(X) %*% X) %*% t(X) %*% y

p1 <- ggplot() +
  geom_line(aes(x = x_0, y = m_gauss), color = mypal[1]) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed") +
  labs(x = "x", y = "m(x)", title = "Gaussian regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p2 <- ggplot() +
  geom_line(aes(x = x_0, y = m_epanechnikov), color = mypal[2]) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed") +
  labs(x = "x", y = "m(x)", title = "Epanechnikov regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p3 <- ggplot() +
  geom_line(aes(x = x_0, y = m_parzen), color = mypal[3]) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed") +
  labs(x = "x", y = "m(x)", title = "Parzen regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p4 <- ggplot() +
  geom_line(aes(x = x_0, y = m_bartlett), color = mypal[4]) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed") +
  labs(x = "x", y = "m(x)", title = "Bartlett regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.2.b.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.b.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.b.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.b.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)

## 2. (c)
m_gauss <- vector(mode = "numeric", N)
m_epanechnikov <- vector(mode = "numeric", N)
m_parzen <- vector(mode = "numeric", N)
m_bartlett <- vector(mode = "numeric", N)

for (k in 1:N) {
  z_x <- (x - x_0[k])/h_x;
  k_x <- gauss(z_x);
  k1_x <- x_0[k] * gauss(z_x);
  k2_x <- x_0[k]^2 * gauss(z_x);
  m_gauss[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/h_x;
  k_x <- epanechnikov(z_x);
  k1_x <- x_0[k] * epanechnikov(z_x);
  k2_x <- x_0[k]^2 * epanechnikov(z_x);
  m_epanechnikov[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/h_x;
  k_x <- parzen(z_x);
  k1_x <- x_0[k] * parzen(z_x);
  k2_x <- x_0[k]^2 * parzen(z_x);
  m_parzen[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/h_x;
  k_x <- bartlett(z_x);
  k1_x <- x_0[k] * bartlett(z_x);
  k2_x <- x_0[k]^2 * bartlett(z_x);
  m_bartlett[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
}

p1 <- ggplot() +
  geom_line(aes(x = x_0, y = m_gauss), color = mypal[1]) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed") +
  labs(x = "x", y = "m(x)", title = "Gaussian regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p2 <- ggplot() +
  geom_line(aes(x = x_0, y = m_epanechnikov), color = mypal[2]) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed") +
  labs(x = "x", y = "m(x)", title = "Epanechnikov regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p3 <- ggplot() +
  geom_line(aes(x = x_0, y = m_parzen), color = mypal[3]) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed") +
  labs(x = "x", y = "m(x)", title = "Parzen regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p4 <- ggplot() +
  geom_line(aes(x = x_0, y = m_bartlett), color = mypal[4]) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed") +
  labs(x = "x", y = "m(x)", title = "Bartlett regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.2.c.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.c.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.c.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.c.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)

## 2. (d)
# CV function
cv <- function(h, kernel){
  m_loo <- vector("numeric", length(x))
  
  for (k in 1:length(x)) {
    z_x <- (x[-k] - x[k])/h
    k_x <- kernel(z_x)
    m_loo[k] <- sum(k_x * y[-k]) / sum(k_x)
  }
  
  sum((m_loo - y)^2)
}

optim_gauss <- optimize(cv, c(0, 10), tol = 0.01, kernel = gauss) # hcv = 2.073
optim_epanechnikov <- optimize(cv, c(0, 10), tol = 0.01, kernel = epanechnikov) # hcv = 4.631
optim_parzen <- optimize(cv, c(0, 10), tol = 0.01, kernel = parzen) # hcv = 7.042
optim_bartlett <- optimize(cv, c(0, 10), tol = 0.01, kernel = bartlett) # hcv = 5.000

library(stargazer)
tab.2.d <- data.frame(row.names = c("Gauss", "Epanechnikov", "Parzen", "Bartlett"),
                      h_n = h_x,
                      h_cv = c(optim_gauss$minimum, optim_epanechnikov$minimum, optim_parzen$minimum, optim_bartlett$minimum))

stargazer(tab.2.d, type = "text", summary = FALSE, float = FALSE, digits = 3, align = TRUE, out = "output/tabs/Part IV/tag.2.d.tex")

mcv_gauss <- vector(mode = "numeric", N)
mcv_epanechnikov <- vector(mode = "numeric", N)
mcv_parzen <- vector(mode = "numeric", N)
mcv_bartlett <- vector(mode = "numeric", N)

for (k in 1:N) {
  z_x <- (x - x_0[k])/optim_gauss$minimum;
  k_x <- gauss(z_x);
  k1_x <- x_0[k] * gauss(z_x);
  k2_x <- x_0[k]^2 * gauss(z_x);
  mcv_gauss[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/optim_epanechnikov$minimum;
  k_x <- epanechnikov(z_x);
  k1_x <- x_0[k] * epanechnikov(z_x);
  k2_x <- x_0[k]^2 * epanechnikov(z_x);
  mcv_epanechnikov[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/optim_parzen$minimum;
  k_x <- parzen(z_x);
  k1_x <- x_0[k] * parzen(z_x);
  k2_x <- x_0[k]^2 * parzen(z_x);
  mcv_parzen[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/optim_bartlett$minimum;
  k_x <- bartlett(z_x);
  k1_x <- x_0[k] * bartlett(z_x);
  k2_x <- x_0[k]^2 * bartlett(z_x);
  mcv_bartlett[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
}

p1 <- ggplot() +
  geom_line(aes(x = x_0, y = m_gauss), color = mypal[1], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_gauss), color = mypal[2]) +
  labs(x = "x", y = "m(x)", title = "Gauss regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p2 <- ggplot() +
  geom_line(aes(x = x_0, y = m_epanechnikov), color = mypal[3], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_epanechnikov), color = mypal[4]) +
  labs(x = "x", y = "m(x)", title = "Epanechnikov regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p3 <- ggplot() +
  geom_line(aes(x = x_0, y = m_parzen), color = mypal[5], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_parzen), color = mypal[6]) +
  labs(x = "x", y = "m(x)", title = "Parzen regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p4 <- ggplot() +
  geom_line(aes(x = x_0, y = m_bartlett), color = mypal[7], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_bartlett), color = mypal[8]) +
  labs(x = "x", y = "m(x)", title = "Bartlett regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.2.d.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.d.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.d.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.2.d.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)

## 3. (b)
cv_w <- function(h, kernel){
  m_loo <- vector("numeric", length(x))
  
  for (k in 1:length(x)) {
    z_x <- (x[-k] - x[k])/h
    k_x <- kernel(z_x)
    m_loo[k] <- sum(k_x * w[-k]) / sum(k_x)
  }
  
  sum((m_loo - w)^2)
}

optim_gauss_w <- optimize(cv_w, c(0, 10), tol = 0.01, kernel = gauss) # hcv = 2.44
optim_epanechnikov_w <- optimize(cv_w, c(0, 10), tol = 0.01, kernel = epanechnikov) # hcv = 5.91
optim_parzen_w <- optimize(cv_w, c(0, 10), tol = 0.01, kernel = parzen) # hcv = 8.41
optim_bartlett_w <- optimize(cv_w, c(0, 10), tol = 0.01, kernel = bartlett) # hcv = 6.27

tab.3.b <- data.frame(row.names = c("Gauss", "Epanechnikov", "Parzen", "Bartlett"),
                      h_n = h_x,
                      h_cv = c(optim_gauss_w$minimum, optim_epanechnikov_w$minimum, optim_parzen_w$minimum, optim_bartlett_w$minimum))

stargazer(tab.3.b, type = "text", summary = FALSE, float = FALSE, digits = 3, align = TRUE, out = "output/tabs/Part IV/tag.3.b.tex")

x_0 <- unique(x)
N <- length(x_0)

m_gauss <- vector(mode = "numeric", N)
m_epanechnikov <- vector(mode = "numeric", N)
m_parzen <- vector(mode = "numeric", N)
m_bartlett <- vector(mode = "numeric", N)

m_gauss_w <- vector(mode = "numeric", N)
m_epanechnikov_w <- vector(mode = "numeric", N)
m_parzen_w <- vector(mode = "numeric", N)
m_bartlett_w <- vector(mode = "numeric", N)

for (k in 1:N) {
  z_x <- gauss((x - x_0[k])/h_x);
  m_gauss[k] <- sum(z_x * y) / sum(z_x);
  m_gauss_w[k] <- sum(z_x * w) / sum(z_x);
  
  z_x <- epanechnikov((x - x_0[k])/h_x);
  m_epanechnikov[k] <- sum(z_x * y) / sum(z_x);
  m_epanechnikov_w[k] <- sum(z_x * w) / sum(z_x);
  
  z_x <- parzen((x - x_0[k])/h_x);
  m_parzen[k] <- sum(z_x * y) / sum(z_x);
  m_parzen_w[k] <- sum(z_x * w) / sum(z_x);
  
  z_x <- bartlett((x - x_0[k])/h_x);
  m_bartlett[k] <- sum(z_x * y) / sum(z_x);
  m_bartlett_w[k] <- sum(z_x * w) / sum(z_x);
}

mcv_gauss <- vector(mode = "numeric", N)
mcv_epanechnikov <- vector(mode = "numeric", N)
mcv_parzen <- vector(mode = "numeric", N)
mcv_bartlett <- vector(mode = "numeric", N)

mcv_gauss_w <- vector(mode = "numeric", N)
mcv_epanechnikov_w <- vector(mode = "numeric", N)
mcv_parzen_w <- vector(mode = "numeric", N)
mcv_bartlett_w <- vector(mode = "numeric", N)

for (k in 1:N) {
  z_x <- gauss((x - x_0[k])/optim_gauss$minimum);
  mcv_gauss[k] <- sum(z_x * y) / sum(z_x);
  z_x <- gauss((x - x_0[k])/optim_gauss_w$minimum);
  mcv_gauss_w[k] <- sum(z_x * w) / sum(z_x);
  
  z_x <- epanechnikov((x - x_0[k])/optim_epanechnikov$minimum);
  mcv_epanechnikov[k] <- sum(z_x * y) / sum(z_x);
  z_x <- epanechnikov((x - x_0[k])/optim_epanechnikov_w$minimum);
  mcv_epanechnikov_w[k] <- sum(z_x * w) / sum(z_x);
  
  z_x <- parzen((x - x_0[k])/optim_parzen$minimum);
  mcv_parzen[k] <- sum(z_x * y) / sum(z_x);
  z_x <- gauss((x - x_0[k])/optim_parzen_w$minimum);
  mcv_parzen_w[k] <- sum(z_x * w) / sum(z_x);
  
  z_x <- bartlett((x - x_0[k])/optim_bartlett$minimum);
  mcv_bartlett[k] <- sum(z_x * y) / sum(z_x);
  z_x <- bartlett((x - x_0[k])/optim_bartlett_w$minimum);
  mcv_bartlett_w[k] <- sum(z_x * w) / sum(z_x);
}

p1 <- ggplot() +
  geom_line(aes(x = x_0, y = m_gauss), color = mypal[1], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_gauss), color = mypal[2]) +
  labs(x = "x", y = "E[y|x]", title = "Gauss regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p2 <- ggplot() +
  geom_line(aes(x = x_0, y = m_epanechnikov), color = mypal[3], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_epanechnikov), color = mypal[4]) +
  labs(x = "x", y = "E[y|x]", title = "Epanechnikov regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p3 <- ggplot() +
  geom_line(aes(x = x_0, y = m_parzen), color = mypal[5], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_parzen), color = mypal[6]) +
  labs(x = "x", y = "E[y|x]", title = "Parzen regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p4 <- ggplot() +
  geom_line(aes(x = x_0, y = m_bartlett), color = mypal[7], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_bartlett), color = mypal[8]) +
  labs(x = "x", y = "E[y|x]", title = "Bartlett regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(4, 7)) +
  theme_few(base_size = 12, base_family = "sans")

p5 <- ggplot() +
  geom_line(aes(x = x_0, y = m_gauss_w), color = mypal[1], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_gauss_w), color = mypal[2]) +
  labs(x = "x", y = "E[w|x]", title = "Gauss regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(6, 10)) +
  theme_few(base_size = 12, base_family = "sans")

p6 <- ggplot() +
  geom_line(aes(x = x_0, y = m_epanechnikov_w), color = mypal[3], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_epanechnikov_w), color = mypal[4]) +
  labs(x = "x", y = "E[w|x]", title = "Epanechnikov regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(6, 10)) +
  theme_few(base_size = 12, base_family = "sans")

p7 <- ggplot() +
  geom_line(aes(x = x_0, y = m_parzen_w), color = mypal[5], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_parzen_w), color = mypal[6]) +
  labs(x = "x", y = "E[w|x]", title = "Parzen regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(6, 10)) +
  theme_few(base_size = 12, base_family = "sans")

p8 <- ggplot() +
  geom_line(aes(x = x_0, y = m_bartlett_w), color = mypal[7], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_bartlett_w), color = mypal[8]) +
  labs(x = "x", y = "E[w|x]", title = "Bartlett regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(6, 10)) +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.3.b.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.b.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.b.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.b.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.b.5.pdf", plot = p5, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.b.6.pdf", plot = p6, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.b.7.pdf", plot = p7, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.b.8.pdf", plot = p8, path = "output/figs/Part IV/", width = 8, height = 8)

## 3. (c)
e1_gauss <- vector("numeric", length(y))
e2_gauss <- vector("numeric", length(y))

e1_epanechnikov <- vector("numeric", length(y))
e2_epanechnikov <- vector("numeric", length(y))

e1_parzen <- vector("numeric", length(y))
e2_parzen <- vector("numeric", length(y))

e1_bartlett <- vector("numeric", length(y))
e2_bartlett <- vector("numeric", length(y))

e1_cv_gauss <- vector("numeric", length(y))
e2_cv_gauss <- vector("numeric", length(y))

e1_cv_epanechnikov <- vector("numeric", length(y))
e2_cv_epanechnikov <- vector("numeric", length(y))

e1_cv_parzen <- vector("numeric", length(y))
e2_cv_parzen <- vector("numeric", length(y))

e1_cv_bartlett <- vector("numeric", length(y))
e2_cv_bartlett <- vector("numeric", length(y))

for (k in 1:length(y)) {
  e1_gauss[k] <- y[k] - m_gauss[x[k] == x_0]
  e2_gauss[k] <- w[k] - m_gauss_w[x[k] == x_0]
  
  e1_epanechnikov[k] <- y[k] - m_epanechnikov[x[k] == x_0]
  e2_epanechnikov[k] <- w[k] - m_epanechnikov_w[x[k] == x_0]
  
  e1_parzen[k] <- y[k] - m_parzen[x[k] == x_0]
  e2_parzen[k] <- w[k] - m_parzen_w[x[k] == x_0]
  
  e1_bartlett[k] <- y[k] - m_bartlett[x[k] == x_0]
  e2_bartlett[k] <- w[k] - m_bartlett_w[x[k] == x_0]
  
  e1_cv_gauss[k] <- y[k] - mcv_gauss[x[k] == x_0]
  e2_cv_gauss[k] <- w[k] - mcv_gauss_w[x[k] == x_0]
  
  e1_cv_epanechnikov[k] <- y[k] - mcv_epanechnikov[x[k] == x_0]
  e2_cv_epanechnikov[k] <- w[k] - mcv_epanechnikov_w[x[k] == x_0]
  
  e1_cv_parzen[k] <- y[k] - mcv_parzen[x[k] == x_0]
  e2_cv_parzen[k] <- w[k] - mcv_parzen_w[x[k] == x_0]
  
  e1_cv_bartlett[k] <- y[k] - mcv_bartlett[x[k] == x_0]
  e2_cv_bartlett[k] <- w[k] - mcv_bartlett_w[x[k] == x_0]
}

calculate_beta <- function(e1, e2){
  sum(e1 * e2) / sum(e2^2)
}

beta_rot <- c(calculate_beta(e1_gauss, e2_gauss),
              calculate_beta(e1_epanechnikov, e2_epanechnikov),
              calculate_beta(e1_parzen, e2_parzen),
              calculate_beta(e1_bartlett, e2_bartlett))

beta_cv <- c(calculate_beta(e1_cv_gauss, e2_cv_gauss),
             calculate_beta(e1_cv_epanechnikov, e2_cv_epanechnikov),
             calculate_beta(e1_cv_parzen, e2_cv_parzen),
             calculate_beta(e1_cv_bartlett, e2_cv_bartlett))

tab.3.c <- data.frame(row.names = c("gauss", "epanechnikov", "parzen", "bartlett"),
                      beta_rot,
                      beta_cv)

stargazer(tab.3.c, type = "text", summary = FALSE, float = FALSE, digits = 3, align = TRUE, out = "output/tabs/Part IV/tag.3.c.tex")

## 3. (d)
m_gauss_pl <- vector(mode = "numeric", N)
m_epanechnikov_pl <- vector(mode = "numeric", N)
m_parzen_pl <- vector(mode = "numeric", N)
m_bartlett_pl <- vector(mode = "numeric", N)

for (k in 1:N) {
  z_x <- gauss((x - x_0[k])/h_x);
  m_gauss_pl[k] <- sum(z_x * (y - beta_rot[1] * w)) / sum(z_x);
  
  z_x <- epanechnikov((x - x_0[k])/h_x);
  m_epanechnikov_pl[k] <- sum(z_x * (y - beta_rot[2] * w)) / sum(z_x);
  
  z_x <- parzen((x - x_0[k])/h_x);
  m_parzen_pl[k] <- sum(z_x * (y - beta_rot[3] * w)) / sum(z_x);
  
  z_x <- bartlett((x - x_0[k])/h_x);
  m_bartlett_pl[k] <- sum(z_x * (y - beta_rot[4] * w)) / sum(z_x);
}

mcv_gauss_pl <- vector(mode = "numeric", N)
mcv_epanechnikov_pl <- vector(mode = "numeric", N)
mcv_parzen_pl <- vector(mode = "numeric", N)
mcv_bartlett_pl <- vector(mode = "numeric", N)

for (k in 1:N) {
  z_x <- gauss((x - x_0[k])/optim_gauss$minimum);
  mcv_gauss_pl[k] <- sum(z_x * (y - beta_cv[1] * w)) / sum(z_x);
  
  z_x <- epanechnikov((x - x_0[k])/optim_epanechnikov$minimum);
  mcv_epanechnikov_pl[k] <- sum(z_x * (y - beta_cv[2] * w)) / sum(z_x);
  
  z_x <- parzen((x - x_0[k])/optim_parzen$minimum);
  mcv_parzen_pl[k] <- sum(z_x * (y - beta_cv[3] * w)) / sum(z_x);
  
  z_x <- bartlett((x - x_0[k])/optim_bartlett$minimum);
  mcv_bartlett_pl[k] <- sum(z_x * (y - beta_cv[4] * w)) / sum(z_x);
}

p1 <- ggplot() +
  geom_line(aes(x = x_0, y = m_gauss), color = mypal[1], linetype = "dashed") +
  geom_line(aes(x = x_0, y = m_gauss_pl), color = mypal[2]) +
  labs(x = "x", y = "m(x)", title = "Gauss regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(0, 8)) +
  theme_few(base_size = 12, base_family = "sans")

p2 <- ggplot() +
  geom_line(aes(x = x_0, y = m_epanechnikov), color = mypal[3], linetype = "dashed") +
  geom_line(aes(x = x_0, y = m_epanechnikov_pl), color = mypal[4]) +
  labs(x = "x", y = "m(x)", title = "Epanechnikov regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(0, 8)) +
  theme_few(base_size = 12, base_family = "sans")

p3 <- ggplot() +
  geom_line(aes(x = x_0, y = m_parzen), color = mypal[5], linetype = "dashed") +
  geom_line(aes(x = x_0, y = m_parzen_pl), color = mypal[6]) +
  labs(x = "x", y = "m(x)", title = "Parzen regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(0, 8)) +
  theme_few(base_size = 12, base_family = "sans")

p4 <- ggplot() +
  geom_line(aes(x = x_0, y = m_bartlett), color = mypal[7], linetype = "dashed") +
  geom_line(aes(x = x_0, y = m_bartlett_pl), color = mypal[8]) +
  labs(x = "x", y = "m(x)", title = "Bartlett regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(0, 8)) +
  theme_few(base_size = 12, base_family = "sans")

p5 <- ggplot() +
  geom_line(aes(x = x_0, y = mcv_gauss), color = mypal[1], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_gauss_pl), color = mypal[2]) +
  labs(x = "x", y = "m(x)", title = "Gauss regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(0, 8)) +
  theme_few(base_size = 12, base_family = "sans")

p6 <- ggplot() +
  geom_line(aes(x = x_0, y = mcv_epanechnikov), color = mypal[3], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_epanechnikov_pl), color = mypal[4]) +
  labs(x = "x", y = "m(x)", title = "Epanechnikov regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(0, 8)) +
  theme_few(base_size = 12, base_family = "sans")

p7 <- ggplot() +
  geom_line(aes(x = x_0, y = mcv_parzen), color = mypal[5], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_parzen_pl), color = mypal[6]) +
  labs(x = "x", y = "m(x)", title = "Parzen regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(0, 8)) +
  theme_few(base_size = 12, base_family = "sans")

p8 <- ggplot() +
  geom_line(aes(x = x_0, y = mcv_bartlett), color = mypal[7], linetype = "dashed") +
  geom_line(aes(x = x_0, y = mcv_bartlett_pl), color = mypal[8]) +
  labs(x = "x", y = "m(x)", title = "Bartlett regression function") +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(0, 8)) +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.3.d.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.d.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.d.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.d.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.d.5.pdf", plot = p5, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.d.6.pdf", plot = p6, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.d.7.pdf", plot = p7, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.3.d.8.pdf", plot = p8, path = "output/figs/Part IV/", width = 8, height = 8)

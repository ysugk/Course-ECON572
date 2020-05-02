## Data import
mat <- R.matlab::readMat("input/Part IV/PS2/data.mat")
data <- mat$data
z <- data[ ,1] # (nX1) vector (age)
y <- data[ ,6] # (nX1) vector (happiness)
x <- data[ ,7] # (nX1) vector (income)
sex <- data[ ,8]
emp <- data[ ,9]
educ <- data[ , 10]

x_max <- quantile(x, 0.95)
x_min <- quantile(x, 0.05)
dummy <- (x <= x_max) & (x >= x_min)

y <- y[dummy]
x <- x[dummy]
z <- z[dummy]
sex <- sex[dummy]
emp <- emp[dummy]
educ <- educ[dummy]
x <- log(x)

# rule-of-thumb h
h_x <- sd(x) * (length(x)^(-1/5))

# kernel functions
gauss <- function(u) exp((-1/2) * u^2) / sqrt(2*pi);
epanechnikov <- function(u) (3/4) * (1 - u^2) * (abs(u) <= 1);
parzen <- function(u) (1 - 6*u^2 + 6*abs(u)^3) * (abs(u) <= 1/2) + (2*(1 - abs(u))^3) * (1/2 < abs(u)) * (abs(u) <= 1);
bartlett <- function(u) (1 - abs(u)) * (abs(u) <= 1);

## 1. (a), (b)

# cross-validation
# local constant
cv_lc <- function(h, x, y,  kernel){
  m_loo <- vector("numeric", length(x))
  
  for (k in 1:length(x)) {
    z_x <- (x[-k] - x[k])/h
    k_x <- kernel(z_x)
    m_loo[k] <- sum(k_x * y[-k]) / sum(k_x)
  }
  
  sum((m_loo - y)^2)
}

# local linear
cv_ll <- function(h, x, y,  kernel){
  m_loo <- vector("numeric", length(x))
  
  for (k in 1:length(x)) {
    z_x <- (x[-k] - x[k])/h
    k_x <- kernel(z_x)
    k1_x <- x[-k] * kernel(z_x)
    k2_x <- x[-k]^2 * kernel(z_x)
    m_loo[k] <- ((sum(k_x * y[-k]) - sum(k1_x) * sum(k1_x * y[-k]))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  }
  
  sum((m_loo - y)^2)
}

# h set by minimizing cv
# local constant
# optim_lc_gauss <- optimize(f = cv_lc, x = x, y = y, kernel = gauss, lower = 0, upper = 2) # 0.121
h_lc_gauss <- 0.121
# optim_lc_epanechnikov <- optimize(f = cv_lc, x = x, y = y, kernel = epanechnikov, lower = 0, upper = 2) # 0.284
h_lc_epanechnikov <- 0.284
# optim_lc_parzen <- optimize(f = cv_lc, x = x, y = y, kernel = parzen, lower = 0, upper = 2) # 0.399
h_lc_parzen <- 0.399
# optim_lc_bartlett <- optimize(f = cv_lc, x = x, y = y, kerCnel = bartlett, lower = 0, upper = 2) # 0.182
h_lc_bartlett <- 0.182

# local linear
# optim_ll_gauss <- optimize(f = cv_ll, x = x, y = y, kernel = gauss, lower = 0, upper = 2) # 0.120
h_ll_gauss <- 0.120
# optim_ll_epanechnikov <- optimize(f = cv_ll, x = x, y = y, kernel = epanechnikov, lower = 0, upper = 2) # 0.284
h_ll_epanechnikov <- 0.284
# optim_ll_parzen <- optimize(f = cv_ll, x = x, y = y, kernel = parzen, lower = 0, upper = 2) # 0.399
h_ll_parzen <- 0.399
# optim_ll_bartlett <- optimize(f = cv_ll, x = x, y = y, kernel = bartlett, lower = 0, upper = 2) # 0.182
h_ll_bartlett <- 0.182

# estimate m
x_0 <- seq(from = 6, to = 10, by = h_x/2)
N <- length(x_0)

mlc_rot_gauss <- vector(mode = "numeric", N)
mlc_rot_epanechnikov <- vector(mode = "numeric", N)
mlc_rot_parzen <- vector(mode = "numeric", N)
mlc_rot_bartlett <- vector(mode = "numeric", N)

mll_rot_gauss <- vector(mode = "numeric", N)
mll_rot_epanechnikov <- vector(mode = "numeric", N)
mll_rot_parzen <- vector(mode = "numeric", N)
mll_rot_bartlett <- vector(mode = "numeric", N)

mlc_cv_gauss <- vector(mode = "numeric", N)
mlc_cv_epanechnikov <- vector(mode = "numeric", N)
mlc_cv_parzen <- vector(mode = "numeric", N)
mlc_cv_bartlett <- vector(mode = "numeric", N)

mll_cv_gauss <- vector(mode = "numeric", N)
mll_cv_epanechnikov <- vector(mode = "numeric", N)
mll_cv_parzen <- vector(mode = "numeric", N)
mll_cv_bartlett <- vector(mode = "numeric", N)

for (k in 1:N) {
  # gauss
  z_x <- (x - x_0[k])/h_x;
  k_x <- gauss(z_x);
  k1_x <- x_0[k] * gauss(z_x);
  k2_x <- x_0[k]^2 * gauss(z_x);
  mlc_rot_gauss[k] <- sum(k_x * y) / sum(k_x);
  mll_rot_gauss[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/h_lc_gauss;
  k_x <- gauss(z_x);
  k1_x <- x_0[k] * gauss(z_x);
  k2_x <- x_0[k]^2 * gauss(z_x);
  mlc_cv_gauss[k] <- sum(k_x * y) / sum(k_x);
  
  z_x <- (x - x_0[k])/h_ll_gauss;
  k_x <- gauss(z_x);
  k1_x <- x_0[k] * gauss(z_x);
  k2_x <- x_0[k]^2 * gauss(z_x);
  mll_cv_gauss[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  # epanechnikov
  z_x <- (x - x_0[k])/h_x;
  k_x <- epanechnikov(z_x);
  k1_x <- x_0[k] * epanechnikov(z_x);
  k2_x <- x_0[k]^2 * epanechnikov(z_x);
  mlc_rot_epanechnikov[k] <- sum(k_x * y) / sum(k_x);
  mll_rot_epanechnikov[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/h_lc_epanechnikov;
  k_x <- epanechnikov(z_x);
  k1_x <- x_0[k] * epanechnikov(z_x);
  k2_x <- x_0[k]^2 * epanechnikov(z_x);
  mlc_cv_epanechnikov[k] <- sum(k_x * y) / sum(k_x);
  
  z_x <- (x - x_0[k])/h_ll_epanechnikov;
  k_x <- epanechnikov(z_x);
  k1_x <- x_0[k] * epanechnikov(z_x);
  k2_x <- x_0[k]^2 * epanechnikov(z_x);
  mll_cv_epanechnikov[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  # parzen
  z_x <- (x - x_0[k])/h_x;
  k_x <- parzen(z_x);
  k1_x <- x_0[k] * parzen(z_x);
  k2_x <- x_0[k]^2 * parzen(z_x);
  mlc_rot_parzen[k] <- sum(k_x * y) / sum(k_x);
  mll_rot_parzen[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/h_lc_parzen;
  k_x <- parzen(z_x);
  k1_x <- x_0[k] * parzen(z_x);
  k2_x <- x_0[k]^2 * parzen(z_x);
  mlc_cv_parzen[k] <- sum(k_x * y) / sum(k_x);

  
  z_x <- (x - x_0[k])/h_ll_parzen;
  k_x <- parzen(z_x);
  k1_x <- x_0[k] * parzen(z_x);
  k2_x <- x_0[k]^2 * parzen(z_x);
  mll_cv_parzen[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  # bartlett
  z_x <- (x - x_0[k])/h_x;
  k_x <- bartlett(z_x);
  k1_x <- x_0[k] * bartlett(z_x);
  k2_x <- x_0[k]^2 * bartlett(z_x);
  mlc_rot_bartlett[k] <- sum(k_x * y) / sum(k_x);
  mll_rot_bartlett[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
  
  z_x <- (x - x_0[k])/h_lc_bartlett;
  k_x <- bartlett(z_x);
  k1_x <- x_0[k] * bartlett(z_x);
  k2_x <- x_0[k]^2 * bartlett(z_x);
  mlc_cv_bartlett[k] <- sum(k_x * y) / sum(k_x);
  
  z_x <- (x - x_0[k])/h_ll_bartlett;
  k_x <- bartlett(z_x);
  k1_x <- x_0[k] * bartlett(z_x);
  k2_x <- x_0[k]^2 * bartlett(z_x);
  mll_cv_bartlett[k] <- ((sum(k_x * y) - sum(k1_x) * sum(k1_x * y))/sum(k2_x))/((sum(k_x) - sum(k1_x)^2)/sum(k2_x));
}

# ols
X <- cbind(rep(1, length(x)), x)
b <- solve(t(X) %*% X) %*% t(X) %*% y

# draw plot
library(tidyverse)
library(ggthemes)

# local constant
p1 <- tibble::tibble(x_0,
               "Gaussian" = mlc_rot_gauss, 
               "Epanechnikov" = mlc_rot_epanechnikov, 
               "Parzen" = mlc_rot_parzen, 
               "Bartlett" = mlc_rot_bartlett) %>%
  tidyr::gather("Kernel", "value", -x_0) %>%
  ggplot() +
  geom_line(aes(x = x_0, y = value, color = Kernel)) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed", color = "blue") +
  labs(x = "x", y = "m(x)", title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p2 <- tibble::tibble(x_0,
                     "Gaussian" = mlc_cv_gauss, 
                     "Epanechnikov" = mlc_cv_epanechnikov, 
                     "Parzen" = mlc_cv_parzen, 
                     "Bartlett" = mlc_cv_bartlett) %>%
  tidyr::gather("Kernel", "value", -x_0) %>%
  ggplot() +
  geom_line(aes(x = x_0, y = value, color = Kernel)) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed", color = "blue") +
  labs(x = "x", y = "m(x)", title = "Local constant (Cross validation)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

# local linear
p3 <- tibble::tibble(x_0,
                     "Gaussian" = mll_rot_gauss, 
                     "Epanechnikov" = mll_rot_epanechnikov, 
                     "Parzen" = mll_rot_parzen, 
                     "Bartlett" = mll_rot_bartlett) %>%
  tidyr::gather("Kernel", "value", -x_0) %>%
  ggplot() +
  geom_line(aes(x = x_0, y = value, color = Kernel)) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed", color = "blue") +
  labs(x = "x", y = "m(x)", title = "Local linear (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p4 <- tibble::tibble(x_0,
                     "Gaussian" = mll_cv_gauss, 
                     "Epanechnikov" = mll_cv_epanechnikov, 
                     "Parzen" = mll_cv_parzen, 
                     "Bartlett" = mll_cv_bartlett) %>%
  tidyr::gather("Kernel", "value", -x_0) %>%
  ggplot() +
  geom_line(aes(x = x_0, y = value, color = Kernel)) +
  geom_abline(slope = b[2], intercept = b[1], linetype = "dashed", color = "blue") +
  labs(x = "x", y = "m(x)", title = "Local linear (Cross-validation)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.4.a.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.4.a.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.4.a.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.4.a.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)

## 2. (a)
# rule-of-thumb h
h_z <- sd(z) * (length(z)^(-1/5))

z_0 <- seq(from = 20, to = 90, by = h_z/2)
N <- length(z_0)
X <- cbind(rep(1, length(x)), x)

fcm_lc_rot_gauss <- matrix(NA, N, ncol(X))
fcm_lc_rot_epanechnikov <- matrix(NA, N, ncol(X))
fcm_lc_rot_parzen <- matrix(NA, N, ncol(X))
fcm_lc_rot_bartlett <- matrix(NA, N, ncol(X))

fcm_ll_rot_gauss <- matrix(NA, N, ncol(X))
fcm_ll_rot_epanechnikov <- matrix(NA, N, ncol(X))
fcm_ll_rot_parzen <- matrix(NA, N, ncol(X))
fcm_ll_rot_bartlett <- matrix(NA, N, ncol(X))

for (k in 1:N) {
  # gauss
  z_z <- (z - z_0[k])/h_z;
  k_z <- gauss(z_z);
  augX <- cbind(X, X * z_z)
  
  fcm_lc_rot_gauss[k, ] <- solve(t(X) %*% diag(k_z) %*% X, t(X) %*% diag(k_z) %*% y)
  fcm_ll_rot_gauss[k, ] <- solve(t(augX) %*% diag(k_z) %*% augX, t(augX) %*% diag(k_z) %*% y)[1:ncol(X)]
  
  # epanechnikov
  z_z <- (z - z_0[k])/h_z;
  k_z <- epanechnikov(z_z);
  augX <- cbind(X, X * z_z)
  
  fcm_lc_rot_epanechnikov[k, ] <- solve(t(X) %*% diag(k_z) %*% X, t(X) %*% diag(k_z) %*% y) 
  fcm_ll_rot_epanechnikov[k, ] <- solve(t(augX) %*% diag(k_z) %*% augX, t(augX) %*% diag(k_z) %*% y)[1:ncol(X)]
  
  # parzen
  z_z <- (z - z_0[k])/h_z;
  k_z <- parzen(z_z);
  augX <- cbind(X, X * z_z)
  
  fcm_lc_rot_parzen[k, ] <- solve(t(X) %*% diag(k_z) %*% X, t(X) %*% diag(k_z) %*% y)
  fcm_ll_rot_parzen[k, ] <- solve(t(augX) %*% diag(k_z) %*% augX, t(augX) %*% diag(k_z) %*% y)[1:ncol(X)]
  
  # bartlett
  z_z <- (z - z_0[k])/h_z;
  k_z <- bartlett(z_z);
  augX <- cbind(X, X * z_z)
  
  fcm_lc_rot_bartlett[k, ] <- solve(t(X) %*% diag(k_z) %*% X, t(X) %*% diag(k_z) %*% y)
  fcm_ll_rot_bartlett[k, ] <- solve(t(augX) %*% diag(k_z) %*% augX, t(augX) %*% diag(k_z) %*% y)[1:ncol(X)]
}

# local constant
p1 <- tibble::tibble(z_0,
                     "Gaussian" = fcm_lc_rot_gauss[ ,1], 
                     "Epanechnikov" = fcm_lc_rot_epanechnikov[ ,1], 
                     "Parzen" = fcm_lc_rot_parzen[ ,1], 
                     "Bartlett" = fcm_lc_rot_bartlett[ ,1]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(alpha), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p2 <- tibble::tibble(z_0,
                     "Gaussian" = fcm_lc_rot_gauss[ ,2], 
                     "Epanechnikov" = fcm_lc_rot_epanechnikov[ ,2], 
                     "Parzen" = fcm_lc_rot_parzen[ ,2], 
                     "Bartlett" = fcm_lc_rot_bartlett[ ,2]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

# local linear
p3 <- tibble::tibble(z_0,
                     "Gaussian" = fcm_ll_rot_gauss[ ,1], 
                     "Epanechnikov" = fcm_ll_rot_epanechnikov[ ,1], 
                     "Parzen" = fcm_ll_rot_parzen[ ,1], 
                     "Bartlett" = fcm_ll_rot_bartlett[ ,1]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(alpha), title = "Local linear (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p4 <- tibble::tibble(z_0,
                     "Gaussian" = fcm_ll_rot_gauss[ ,2], 
                     "Epanechnikov" = fcm_ll_rot_epanechnikov[ ,2], 
                     "Parzen" = fcm_ll_rot_parzen[ ,2], 
                     "Bartlett" = fcm_ll_rot_bartlett[ ,2]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta), title = "Local linear (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.5.a.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.5.a.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.5.a.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.5.a.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)

## 2. (b)
X <- cbind(rep(1, length(x)), x, x^2)

fcmsq_lc_rot_gauss <- matrix(NA, N, ncol(X))
fcmsq_lc_rot_epanechnikov <- matrix(NA, N, ncol(X))
fcmsq_lc_rot_parzen <- matrix(NA, N, ncol(X))
fcmsq_lc_rot_bartlett <- matrix(NA, N, ncol(X))

fcmsq_ll_rot_gauss <- matrix(NA, N, ncol(X))
fcmsq_ll_rot_epanechnikov <- matrix(NA, N, ncol(X))
fcmsq_ll_rot_parzen <- matrix(NA, N, ncol(X))
fcmsq_ll_rot_bartlett <- matrix(NA, N, ncol(X))

for (k in 1:N) {
  # gauss
  z_z <- (z - z_0[k])/h_z;
  k_z <- gauss(z_z);
  augX <- cbind(X, X * z_z)
  
  fcmsq_lc_rot_gauss[k, ] <- solve(t(X) %*% diag(k_z) %*% X, t(X) %*% diag(k_z) %*% y)
  fcmsq_ll_rot_gauss[k, ] <- solve(t(augX) %*% diag(k_z) %*% augX, t(augX) %*% diag(k_z) %*% y)[1:ncol(X)]
  
  # epanechnikov
  z_z <- (z - z_0[k])/h_z;
  k_z <- epanechnikov(z_z);
  augX <- cbind(X, X * z_z)
  
  fcmsq_lc_rot_epanechnikov[k, ] <- solve(t(X) %*% diag(k_z) %*% X, t(X) %*% diag(k_z) %*% y)
  fcmsq_ll_rot_epanechnikov[k, ] <- solve(t(augX) %*% diag(k_z) %*% augX, t(augX) %*% diag(k_z) %*% y)[1:ncol(X)]
  
  # parzen
  z_z <- (z - z_0[k])/h_z;
  k_z <- parzen(z_z);
  augX <- cbind(X, X * z_z)
  
  fcmsq_lc_rot_parzen[k, ] <- solve(t(X) %*% diag(k_z) %*% X, t(X) %*% diag(k_z) %*% y)
  fcmsq_ll_rot_parzen[k, ] <- solve(t(augX) %*% diag(k_z) %*% augX, t(augX) %*% diag(k_z) %*% y)[1:ncol(X)]
  
  # bartlett
  z_z <- (z - z_0[k])/h_z;
  k_z <- bartlett(z_z);
  augX <- cbind(X, X * z_z)
  
  fcmsq_lc_rot_bartlett[k, ] <- solve(t(X) %*% diag(k_z) %*% X, t(X) %*% diag(k_z) %*% y)
  fcmsq_ll_rot_bartlett[k, ] <- solve(t(augX) %*% diag(k_z) %*% augX, t(augX) %*% diag(k_z) %*% y)[1:ncol(X)]
}

# local constant
p1 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss[ ,1], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov[ ,1], 
                     "Parzen" = fcmsq_lc_rot_parzen[ ,1], 
                     "Bartlett" = fcmsq_lc_rot_bartlett[ ,1]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(alpha), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p2 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss[ ,2], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov[ ,2], 
                     "Parzen" = fcmsq_lc_rot_parzen[ ,2], 
                     "Bartlett" = fcmsq_lc_rot_bartlett[ ,2]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta[1]), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p3 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss[ ,3], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov[ ,3], 
                     "Parzen" = fcmsq_lc_rot_parzen[ ,3], 
                     "Bartlett" = fcmsq_lc_rot_bartlett[ ,3]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta[2]), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

# local linear
p4 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_ll_rot_gauss[ ,1], 
                     "Epanechnikov" = fcmsq_ll_rot_epanechnikov[ ,1], 
                     "Parzen" = fcmsq_ll_rot_parzen[ ,1], 
                     "Bartlett" = fcmsq_ll_rot_bartlett[ ,1]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(alpha), title = "Local linear (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p5 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_ll_rot_gauss[ ,2], 
                     "Epanechnikov" = fcmsq_ll_rot_epanechnikov[ ,2], 
                     "Parzen" = fcmsq_ll_rot_parzen[ ,2], 
                     "Bartlett" = fcmsq_ll_rot_bartlett[ ,2]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta[1]), title = "Local linear (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p6 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_ll_rot_gauss[ ,3], 
                     "Epanechnikov" = fcmsq_ll_rot_epanechnikov[ ,3], 
                     "Parzen" = fcmsq_ll_rot_parzen[ ,3], 
                     "Bartlett" = fcmsq_ll_rot_bartlett[ ,3]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta[2]), title = "Local linear (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.5.b.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.5.b.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.5.b.3.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.5.b.4.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.5.b.5.pdf", plot = p5, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.5.b.6.pdf", plot = p6, path = "output/figs/Part IV/", width = 8, height = 8)

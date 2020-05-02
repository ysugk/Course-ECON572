## Data import
mat <- R.matlab::readMat("input/Part IV/PS2/data.mat")
data <- mat$data
z <- data[ ,1] # (nX1) vector (age)
y <- data[ ,6] # (nX1) vector (happiness)
x <- data[ ,7] # (nX1) vector (income)
emp <- data[ ,9]

x_max <- quantile(x, 0.95)
x_min <- quantile(x, 0.05)
dummy <- (x <= x_max) & (x >= x_min)

y <- y[dummy]
x <- x[dummy]
z <- z[dummy]
emp <- emp[dummy]
x <- log(x)

y0 <- y[emp == 0]
x0 <- x[emp == 0]
z0 <- z[emp == 0]

y1 <- y[emp == 1]
x1 <- x[emp == 1]
z1 <- z[emp == 1]

# rule-of-thumb h
h_x <- sd(x) * (length(x)^(-1/5))

# kernel functions
gauss <- function(u) exp((-1/2) * u^2) / sqrt(2*pi);
epanechnikov <- function(u) (3/4) * (1 - u^2) * (abs(u) <= 1);
parzen <- function(u) (1 - 6*u^2 + 6*abs(u)^3) * (abs(u) <= 1/2) + (2*(1 - abs(u))^3) * (1/2 < abs(u)) * (abs(u) <= 1);
bartlett <- function(u) (1 - abs(u)) * (abs(u) <= 1);

## 1. (a), (b)

# estimate m
x_0 <- seq(from = 6, to = 10, by = h_x/2)
N <- length(x_0)

mlc_rot_gauss_0 <- vector(mode = "numeric", N)
mlc_rot_epanechnikov_0 <- vector(mode = "numeric", N)
mlc_rot_parzen_0 <- vector(mode = "numeric", N)
mlc_rot_bartlett_0 <- vector(mode = "numeric", N)

mlc_rot_gauss_1 <- vector(mode = "numeric", N)
mlc_rot_epanechnikov_1 <- vector(mode = "numeric", N)
mlc_rot_parzen_1 <- vector(mode = "numeric", N)
mlc_rot_bartlett_1 <- vector(mode = "numeric", N)


for (k in 1:N) {
  # gauss
  z_x0 <- (x0 - x_0[k])/h_x;
  k_x0 <- gauss(z_x0);
  mlc_rot_gauss_0[k] <- sum(k_x0 * y0) / sum(k_x0);
  
  z_x1 <- (x1 - x_0[k])/h_x;
  k_x1 <- gauss(z_x1);
  mlc_rot_gauss_1[k] <- sum(k_x1 * y1) / sum(k_x1);
  
  # epanechnikov
  z_x0 <- (x0 - x_0[k])/h_x;
  k_x0 <- epanechnikov(z_x0);
  mlc_rot_epanechnikov_0[k] <- sum(k_x0 * y0) / sum(k_x0);
  
  z_x1 <- (x1 - x_0[k])/h_x;
  k_x1 <- epanechnikov(z_x1);
  mlc_rot_epanechnikov_1[k] <- sum(k_x1 * y1) / sum(k_x1);
  
  # parzen
  z_x0 <- (x0 - x_0[k])/h_x;
  k_x0 <- parzen(z_x0);
  mlc_rot_parzen_0[k] <- sum(k_x0 * y0) / sum(k_x0);
  
  z_x1 <- (x1 - x_0[k])/h_x;
  k_x1 <- parzen(z_x1);
  mlc_rot_parzen_1[k] <- sum(k_x1 * y1) / sum(k_x1);
  
  # bartlett
  z_x0 <- (x0 - x_0[k])/h_x;
  k_x0 <- bartlett(z_x0);
  mlc_rot_bartlett_0[k] <- sum(k_x0 * y0) / sum(k_x0);
  
  z_x1 <- (x1 - x_0[k])/h_x;
  k_x1 <- bartlett(z_x1);
  mlc_rot_bartlett_1[k] <- sum(k_x1 * y1) / sum(k_x1);
}

# ols
X0 <- cbind(rep(1, length(x0)), x0)
b0 <- solve(t(X0) %*% X0, t(X0) %*% y0)

X1 <- cbind(rep(1, length(x1)), x1)
b1 <- solve(t(X1) %*% X1, t(X1) %*% y1)


# draw plot
library(tidyverse)
library(ggthemes)

# local constant
p1 <- tibble::tibble(x_0,
                     "Gaussian" = mlc_rot_gauss_0, 
                     "Epanechnikov" = mlc_rot_epanechnikov_0, 
                     "Parzen" = mlc_rot_parzen_0, 
                     "Bartlett" = mlc_rot_bartlett_0) %>%
  tidyr::gather("Kernel", "value", -x_0) %>%
  ggplot() +
  geom_abline(slope = b0[2], intercept = b0[1], linetype = "dashed", color = "blue") +
  geom_line(aes(x = x_0, y = value, color = Kernel)) +
  labs(x = "x", y = "m(x)", title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p2 <- tibble::tibble(x_0,
                     "Gaussian" = mlc_rot_gauss_1, 
                     "Epanechnikov" = mlc_rot_epanechnikov_1, 
                     "Parzen" = mlc_rot_parzen_1, 
                     "Bartlett" = mlc_rot_bartlett_1) %>%
  tidyr::gather("Kernel", "value", -x_0) %>%
  ggplot() +
  geom_abline(slope = b1[2], intercept = b1[1], linetype = "dashed", color = "blue") +
  geom_line(aes(x = x_0, y = value, color = Kernel)) +
  labs(x = "x", y = "m(x)", title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.6.b.1.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.2.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)

## 2. (a)
# rule-of-thumb h
h_z <- sd(z) * (length(z)^(-1/5))

z_0 <- seq(from = 20, to = 90, by = h_z/2)
N <- length(z_0)
X0 <- cbind(rep(1, length(x0)), x0)
X1 <- cbind(rep(1, length(x1)), x1)

fcm_lc_rot_gauss_0 <- matrix(NA, N, ncol(X0))
fcm_lc_rot_epanechnikov_0 <- matrix(NA, N, ncol(X0))
fcm_lc_rot_parzen_0 <- matrix(NA, N, ncol(X0))
fcm_lc_rot_bartlett_0 <- matrix(NA, N, ncol(X0))

fcm_lc_rot_gauss_1 <- matrix(NA, N, ncol(X1))
fcm_lc_rot_epanechnikov_1 <- matrix(NA, N, ncol(X1))
fcm_lc_rot_parzen_1 <- matrix(NA, N, ncol(X1))
fcm_lc_rot_bartlett_1 <- matrix(NA, N, ncol(X1))

for (k in 1:N) {
  # gauss
  z_z0 <- (z0 - z_0[k])/h_z;
  k_z0 <- gauss(z_z0);
  augX0 <- cbind(X0, X0 * z_z0)
  
  fcm_lc_rot_gauss_0[k, ] <- solve(t(X0) %*% diag(k_z0) %*% X0, t(X0) %*% diag(k_z0) %*% y0)
  
  z_z1 <- (z1 - z_0[k])/h_z;
  k_z1 <- gauss(z_z1);
  augX1 <- cbind(X1, X1 * z_z1)
  
  fcm_lc_rot_gauss_1[k, ] <- solve(t(X1) %*% diag(k_z1) %*% X1, t(X1) %*% diag(k_z1) %*% y1)
  
  # epanechnikov
  z_z0 <- (z0 - z_0[k])/h_z;
  k_z0 <- epanechnikov(z_z0);
  augX0 <- cbind(X0, X0 * z_z0)
  
  fcm_lc_rot_epanechnikov_0[k, ] <- solve(t(X0) %*% diag(k_z0) %*% X0, t(X0) %*% diag(k_z0) %*% y0)
  
  z_z1 <- (z1 - z_0[k])/h_z;
  k_z1 <- epanechnikov(z_z1);
  augX1 <- cbind(X1, X1 * z_z1)
  
  fcm_lc_rot_epanechnikov_1[k, ] <- solve(t(X1) %*% diag(k_z1) %*% X1, t(X1) %*% diag(k_z1) %*% y1)
  
  # parzen
  z_z0 <- (z0 - z_0[k])/h_z;
  k_z0 <- parzen(z_z0);
  augX0 <- cbind(X0, X0 * z_z0)
  
  fcm_lc_rot_parzen_0[k, ] <- solve(t(X0) %*% diag(k_z0) %*% X0, t(X0) %*% diag(k_z0) %*% y0)
  
  z_z1 <- (z1 - z_0[k])/h_z;
  k_z1 <- parzen(z_z1);
  augX1 <- cbind(X1, X1 * z_z1)
  
  fcm_lc_rot_parzen_1[k, ] <- solve(t(X1) %*% diag(k_z1) %*% X1, t(X1) %*% diag(k_z1) %*% y1)
  
  # bartlett
  z_z0 <- (z0 - z_0[k])/h_z;
  k_z0 <- bartlett(z_z0);
  augX0 <- cbind(X0, X0 * z_z0)
  
  fcm_lc_rot_bartlett_0[k, ] <- solve(t(X0) %*% diag(k_z0) %*% X0, t(X0) %*% diag(k_z0) %*% y0)
  
  z_z1 <- (z1 - z_0[k])/h_z;
  k_z1 <- bartlett(z_z1);
  augX1 <- cbind(X1, X1 * z_z1)
  
  fcm_lc_rot_bartlett_1[k, ] <- solve(t(X1) %*% diag(k_z1) %*% X1, t(X1) %*% diag(k_z1) %*% y1)
}

# local constant (indicator == 0)
p1 <- tibble::tibble(z_0,
                     "Gaussian" = fcm_lc_rot_gauss_0[ ,1], 
                     "Epanechnikov" = fcm_lc_rot_epanechnikov_0[ ,1], 
                     "Parzen" = fcm_lc_rot_parzen_0[ ,1], 
                     "Bartlett" = fcm_lc_rot_bartlett_0[ ,1]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(alpha), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p2 <- tibble::tibble(z_0,
                     "Gaussian" = fcm_lc_rot_gauss_0[ ,2], 
                     "Epanechnikov" = fcm_lc_rot_epanechnikov_0[ ,2], 
                     "Parzen" = fcm_lc_rot_parzen_0[ ,2], 
                     "Bartlett" = fcm_lc_rot_bartlett_0[ ,2]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

# local constant (indicator == 1)
p3 <- tibble::tibble(z_0,
                     "Gaussian" = fcm_lc_rot_gauss_1[ ,1], 
                     "Epanechnikov" = fcm_lc_rot_epanechnikov_1[ ,1], 
                     "Parzen" = fcm_lc_rot_parzen_1[ ,1], 
                     "Bartlett" = fcm_lc_rot_bartlett_1[ ,1]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(alpha), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p4 <- tibble::tibble(z_0,
                     "Gaussian" = fcm_lc_rot_gauss_1[ ,2], 
                     "Epanechnikov" = fcm_lc_rot_epanechnikov_1[ ,2], 
                     "Parzen" = fcm_lc_rot_parzen_1[ ,2], 
                     "Bartlett" = fcm_lc_rot_bartlett_1[ ,2]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.6.b.3.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.4.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.5.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.6.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)

## 2. (b)
X0 <- cbind(rep(1, length(x0)), x0, x0^2)
X1 <- cbind(rep(1, length(x1)), x1, x1^2)

fcmsq_lc_rot_gauss_0 <- matrix(NA, N, ncol(X0))
fcmsq_lc_rot_epanechnikov_0 <- matrix(NA, N, ncol(X0))
fcmsq_lc_rot_parzen_0 <- matrix(NA, N, ncol(X0))
fcmsq_lc_rot_bartlett_0 <- matrix(NA, N, ncol(X0))

fcmsq_lc_rot_gauss_1 <- matrix(NA, N, ncol(X1))
fcmsq_lc_rot_epanechnikov_1 <- matrix(NA, N, ncol(X1))
fcmsq_lc_rot_parzen_1 <- matrix(NA, N, ncol(X1))
fcmsq_lc_rot_bartlett_1 <- matrix(NA, N, ncol(X1))

for (k in 1:N) {
  # gauss
  z_z0 <- (z0 - z_0[k])/h_z;
  k_z0 <- gauss(z_z0);
  augX0 <- cbind(X0, X0 * z_z0)
  
  fcmsq_lc_rot_gauss_0[k, ] <- solve(t(X0) %*% diag(k_z0) %*% X0, t(X0) %*% diag(k_z0) %*% y0)
  
  z_z1 <- (z1 - z_0[k])/h_z;
  k_z1 <- gauss(z_z1);
  augX1 <- cbind(X1, X1 * z_z1)
  
  fcmsq_lc_rot_gauss_1[k, ] <- solve(t(X1) %*% diag(k_z1) %*% X1, t(X1) %*% diag(k_z1) %*% y1)
  
  # epanechnikov
  z_z0 <- (z0 - z_0[k])/h_z;
  k_z0 <- epanechnikov(z_z0);
  augX0 <- cbind(X0, X0 * z_z0)
  
  fcmsq_lc_rot_epanechnikov_0[k, ] <- solve(t(X0) %*% diag(k_z0) %*% X0, t(X0) %*% diag(k_z0) %*% y0)
  
  z_z1 <- (z1 - z_0[k])/h_z;
  k_z1 <- epanechnikov(z_z1);
  augX1 <- cbind(X1, X1 * z_z1)
  
  fcmsq_lc_rot_epanechnikov_1[k, ] <- solve(t(X1) %*% diag(k_z1) %*% X1, t(X1) %*% diag(k_z1) %*% y1)
  
  # parzen
  z_z0 <- (z0 - z_0[k])/h_z;
  k_z0 <- parzen(z_z0);
  augX0 <- cbind(X0, X0 * z_z0)
  
  fcmsq_lc_rot_parzen_0[k, ] <- solve(t(X0) %*% diag(k_z0) %*% X0, t(X0) %*% diag(k_z0) %*% y0)
  
  z_z1 <- (z1 - z_0[k])/h_z;
  k_z1 <- parzen(z_z1);
  augX1 <- cbind(X1, X1 * z_z1)
  
  fcmsq_lc_rot_parzen_1[k, ] <- solve(t(X1) %*% diag(k_z1) %*% X1, t(X1) %*% diag(k_z1) %*% y1)
  
  # bartlett
  z_z0 <- (z0 - z_0[k])/h_z;
  k_z0 <- bartlett(z_z0);
  augX0 <- cbind(X0, X0 * z_z0)
  
  fcmsq_lc_rot_bartlett_0[k, ] <- solve(t(X0) %*% diag(k_z0) %*% X0, t(X0) %*% diag(k_z0) %*% y0)
  
  z_z1 <- (z1 - z_0[k])/h_z;
  k_z1 <- bartlett(z_z1);
  augX1 <- cbind(X1, X1 * z_z1)
  
  fcmsq_lc_rot_bartlett_1[k, ] <- solve(t(X1) %*% diag(k_z1) %*% X1, t(X1) %*% diag(k_z1) %*% y1)
}

# local constant (indicator == 0)
p1 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss_0[ ,1], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov_0[ ,1], 
                     "Parzen" = fcmsq_lc_rot_parzen_0[ ,1], 
                     "Bartlett" = fcmsq_lc_rot_bartlett_0[ ,1]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(alpha), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p2 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss_0[ ,2], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov_0[ ,2], 
                     "Parzen" = fcmsq_lc_rot_parzen_0[ ,2], 
                     "Bartlett" = fcmsq_lc_rot_bartlett_0[ ,2]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta[1]), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p3 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss_0[ ,3], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov_0[ ,3], 
                     "Parzen" = fcmsq_lc_rot_parzen_0[ ,3], 
                     "Bartlett" = fcmsq_lc_rot_bartlett_0[ ,3]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta[2]), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

# local constant (indicator == 1)
p4 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss_1[ ,1], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov_1[ ,1], 
                     "Parzen" = fcmsq_lc_rot_parzen_1[ ,1], 
                     "Bartlett" = fcmsq_lc_rot_bartlett_1[ ,1]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(alpha), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p5 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss_1[ ,2], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov_1[ ,2], 
                     "Parzen" = fcmsq_lc_rot_parzen_1[ ,2], 
                     "Bartlett" = fcmsq_lc_rot_bartlett_1[ ,2]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta[1]), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

p6 <- tibble::tibble(z_0,
                     "Gaussian" = fcmsq_lc_rot_gauss_1[ ,3], 
                     "Epanechnikov" = fcmsq_lc_rot_epanechnikov_1[ ,3], 
                     "Parzen" = fcmsq_lc_rot_parzen_1[ ,3], 
                     "Bartlett" = fcmsq_lc_rot_bartlett_1[ ,3]) %>%
  filter(z_0 <= 80) %>%
  tidyr::gather("Kernel", "value", -z_0) %>%
  ggplot() +
  geom_line(aes(x = z_0, y = value, color = Kernel)) +
  labs(x = "z", y = expression(beta[2]), title = "Local constant (Rule of thumb)") +
  scale_color_colorblind() +
  theme_few(base_size = 12, base_family = "sans")

ggsave(filename = "fig.6.b.7.pdf", plot = p1, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.8.pdf", plot = p2, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.9.pdf", plot = p3, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.10.pdf", plot = p4, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.11.pdf", plot = p5, path = "output/figs/Part IV/", width = 8, height = 8)
ggsave(filename = "fig.6.b.12.pdf", plot = p6, path = "output/figs/Part IV/", width = 8, height = 8)

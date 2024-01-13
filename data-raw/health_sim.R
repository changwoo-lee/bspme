
## code to prepare `healthdata_simul` dataset goes here
n_y = 2000 # number of health data
data("NO2_Jan2012")
set.seed(1234567)
coords_health = cbind(rnorm(n_y, -95.5, sd = 0.3),rnorm(n_y, 29.8, sd = 0.2))
# simulate health data based on day 10 (Jan 10, 2012)
data_jan10 = NO2_Jan2012[NO2_Jan2012$date == as.POSIXct("2012-01-10"),]

coords_monitor = cbind(data_jan10$lon, data_jan10$lat)
#quilt.plot(coords_monitor, data_jan10$lnNO2)
# visualize
#plot(coords_monitor, pch = 16, col = "purple")
#points(coords_health)
#fields::US(add = T)


sigma <- 1 # standard deviation
a <- 5 # range, in miles
## Matrices needed for prediction

distmat_xx <- fields::rdist.earth(coords_monitor)
distmat_xy <- fields::rdist.earth(coords_monitor, coords_health)
distmat_yy <- fields::rdist.earth(coords_health)

Sigmaxx = fields::Matern(distmat_xx, smoothness = 0.5, range = a, phi = sigma^2)
Sigmaxy = fields::Matern(distmat_xy, smoothness = 0.5, range = a, phi = sigma^2)
Sigmayy = fields::Matern(distmat_yy, smoothness = 0.5, range = a, phi = sigma^2)

## GP mean
pred_mean <- t(Sigmaxy) %*% solve(Sigmaxx, data_jan10$lnNO2)
## GP cov
pred_cov <- Sigmayy - t(Sigmaxy) %*% solve(Sigmaxx,Sigmaxy)

# draw predictive sample
set.seed(1)
true_exposure = mvnfast::rmvn(1, pred_mean, pred_cov)
true_exposure = as.vector(true_exposure)
# par(mfrow = c(1,3))
quilt.plot(coords_health, pred_mean, main = "mean"); US(add= T)
quilt.plot(coords_health, sqrt(diag(pred_cov)), main = "sd"); US(add= T)
quilt.plot(coords_health, true_exposure, main = "draw from MVN"); US(add= T)

# simulated covariate
z = rnorm(n_y)

# generate health data
set.seed(1234)
Y = 1 -2*true_exposure + 3*z + rnorm(n_y, sd = 1)
#summary(Y)
#var(1-2*true_exposure + 3*z) / var(Y)
#summary(lm(Y ~ true_exposure + z))
linpred = -1 + 2*true_exposure - 3*z
Ybinary = rbinom(n_y, 1, prob = 1/(1+exp(-linpred)))
#summary(Ybinary)
#summary(glm(Ybinary ~ true_exposure + z), family = "binomial")

health_sim = data.frame(Y = Y,
                        Ybinary = Ybinary,
                        #X_mean = pred_mean,
                        #X_cov = pred_cov,
                        lon = coords_health[,1],
                        lat = coords_health[,2],
                        Z = z,
                        X_true = true_exposure)


usethis::use_data(health_sim, overwrite = TRUE)



# data
Y <- state.x77[,5]
n <- length(Y)
Z <- cbind(rep(1,n),as.matrix(state.x77[,6:7]))
r <- dim(Z)[2]-1

# least square estimates
beta_hat <- solve(t(Z)%*%Z)%*%t(Z)%*%Y
beta_hat

# R^2 statistic
R_square <- 1 - sum((Y - Z%*%beta_hat)^2)/sum((Y-mean(Y))^2)
R_square

# sigma_hat_square
sigma_hat_square <- sum((Y - Z%*%beta_hat)^2)/(n-r-1)
sigma_hat_square

# estimated covariance of hat{beta}
sigma_hat_square * solve(t(Z)%*%Z)

# t-test for single coefficient
# H_0: beta_j = 0, H_a: beta_j != 0

j <- 1
t_stat <- (beta_hat[j+1] - 0)/sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1])
t_stat

alpha <- 0.05
cval_t <- qt(1-alpha/2, n-r-1)
cval_t

# One-at-a-time confidence interval for beta_j

j <- 1
cat('[',
    beta_hat[j+1] - qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# confidence region based simultaneous confidence intervals 

j <- 2
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# Bonferroni correction based simultaneous confidence intervals

j <- 2
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# F-test
# H_0: beta_1 = beta_2 = 0

C <- matrix(c(0,0,1,0,0,1),2,3)

df_1 <- qr(C)$rank # df_1: rank of matrix R

f_stat <- (t(C%*%beta_hat)%*%solve(C%*%solve(t(Z)%*%Z)%*%t(C))%*%(C%*%beta_hat)/df_1)/sigma_hat_square
f_stat

cval_f <- qf(1-alpha, 2, n-r-1)
cval_f

# (equivalent) F-test by comparing residuals

# fit the reduced model
beta_hat_reduced <- solve(t(Z[,1])%*%Z[,1])%*%t(Z[,1])%*%Y
beta_hat_reduced

f_stat_reduced <- ((sum((Y - Z[,1]%*%beta_hat_reduced)^2) - sum((Y - Z%*%beta_hat)^2))/2)/sigma_hat_square
f_stat_reduced

# confidence interval for z_0^T beta

z_0 <- c(1, 50, 50)

cat('[',
    z_0%*%beta_hat - sqrt(sigma_hat_square)*sqrt(t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ',',
    z_0%*%beta_hat + sqrt(sigma_hat_square)*sqrt(t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ']')

# prediction interval for Y_0 = z_0^T beta + epsilon_0

cat('[',
    z_0%*%beta_hat - sqrt(sigma_hat_square)*sqrt(1+t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ',',
    z_0%*%beta_hat + sqrt(sigma_hat_square)*sqrt(1+t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ']')

# Confidence Region for (beta_1, beta_2)^T

center <- beta_hat[2:3]
es<-eigen(R%*%solve(t(Z)%*%Z)%*%t(R))
e1<-es$vec %*% diag(sqrt(es$val))
r1<-sqrt(df_1*cval_f*sigma_hat_square)
theta<-seq(0,2*pi,len=250)
v1<-cbind(r1*cos(theta), r1*sin(theta))
pts<-t(center - (e1%*%t(v1)))
plot(pts,type="l",main="Confidence Region for (beta_1, beta_2)^T",xlab="beta_1",ylab="beta_2",asp=1,
     xlim = c(-0.3,0.1),ylim=c(-0.1,0.04))
segments(0,center[2],center[1],center[2],lty=2) # highlight the center
segments(center[1],0,center[1],center[2],lty=2)
arrows(-0.3,0,0.1,0)
arrows(0,-0.1,0,0.05)

th2<-c(0,pi/2,pi,3*pi/2,2*pi)   #adding the axis
v2<-cbind(r1*cos(th2), r1*sin(th2))
pts2<-t(center-(e1%*%t(v2)))
segments(pts2[3,1],pts2[3,2],pts2[1,1],pts2[1,2],lty=3)  
segments(pts2[2,1],pts2[2,2],pts2[4,1],pts2[4,2],lty=3)

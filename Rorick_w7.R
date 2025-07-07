#10.1) Refer to Example 10.1 and Figure 10.1. Suppose that we want to test 
#H0 : F= G, where F is the distribution of weight for the casein feed 
#group and G is the distribution of weight for the sunflower feed group 
#of the chickwts data. A test can be based on the two-sample Kolmogorov-Smirnov 
#statistic as shown in Example 10.1. Display a histogram of the permutation 
#replicates of the Kolmogorov-Smirnov two-sample test statistic for this 
#test. Is the test significant at α= 0.10?
attach(chickwts)
x <- sort(weight[feed == "casein"])
y <- sort(weight[feed == "sunflower"])
detach(chickwts)
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:length(z)
reps <- numeric(R) #storage for replicates
t0 <- t.test(x, y)$statistic
for (i in 1:R) {
    k <- sample(K, size = length(x), replace = FALSE) #generate indices k for the first sample
    x1 <- z[k]
    y1 <- z[-k] #complement of x1
    reps[i] <- t.test(x1, y1)$statistic
}
p <- mean(c(t0, reps) >= t0)
hist(reps, main = "T Statistic Histogram", freq = FALSE, xlab = "T", breaks = "scott")
points(t0, 0, cex = 1, pch = 16)


#10.3) Implement the two-sample Cramér-von Mises test for equal distributions 
#as a permutation test using (10.14). Apply the test to the data in Examples 10.1 and 10.2.
library(twosamples)
R=1000
K=24
D = numeric(R)
z <- c(x, y)
D0 <- cvm_stat(x, y)
for (i in 1:R) {
    #generate indices k for the first sample
    k <- sample(1:K, 12, replace = FALSE)
    x1 <- z[k]
    y1 <- z[-k]      #complement of x1
    D[i] <- cvm_stat(x1, y1)
}
p <- mean(c(D0, D) >= D0)
p


#10.4) An rth Nearest Neighbors test statistic for equal distributions: 
#Write a function (for the statistic argument of the boot function) to 
#compute the test statistic Tn,r (10.6). The function syntax should be 
#Tnr(z, ix, sizes, nn) with the data matrix z as its first argument, and an 
#index vector ix as the second argument. The vector of sample sizes sizes 
#and the number of nearest neighbors nn should be the third and fourth 
#arguments. (See the ann function in package yaImpute and Example 10.6.)
library(boot)
library(yaImpute)
NN.idx <- function(x, tree.type="kd", k=NROW(x)) ## function to return the matrix of indices NN_j of nearest neighbors
{
    x <- as.matrix(x)
    k <- min(c(k+1, NROW(x)))
    NN <- yaImpute::ann(ref=x, target=x, tree.type="kd", k=k, verbose=FALSE)
    idx <- NN$knnIndexDist[,1:k]
    nn.idx <- idx[,-1]   #first NN is in column 2
    row.names(nn.idx) <- idx[,1]   #give row names, without this line, it is fine
    nn.idx
}

## function to compute the NN statistic T(n,r)
Tnr <- function(z, ix=1:NROW(z), sizes, nn) 
{
    z <- as.matrix(z)
    n1 <- sizes[1]
    n2 <- sizes[2]
    n <- n1 + n2
    z <- as.matrix(z[ix, ])
    nn.idx <- NN.idx(z, k=nn)
    block1 <- nn.idx[1:n1, ]
    block2 <- nn.idx[(n1+1):n, ]
    i1 <- sum(block1 < n1 + .5)
    i2 <- sum(block2 > n1 + .5)
    return((i1 + i2) / (nn * n))
}


#10.5) The iris data is a four-dimensional distribution with measurements 
#on three species of iris flowers. Using your function Tnr of Exercise 10.4 
#and the boot function, apply your nearest neighbors statistic (r = 2) to test 
#H0 : F = G, where F is the distribution of the iris setosa species, and G is 
#the distribution of the iris virginica species. Repeat the test with r = 3 and 
#r = 4
data(iris)
x <- iris[iris$Species == 'setosa',]
y <- iris[iris$Species == 'virginica',]
z <- as.matrix(rbind(x[, 1:4], y[, 1:4]))
r = 4
N <- c(nrow(x), nrow(y))
boot.obj <- boot(data = z, statistic = Tnr, sim = "permutation", R = 999, sizes = N, nn = r)
tb <- c(boot.obj$t, boot.obj$t0)
p <- mean(tb >= boot.obj$t0)
p
hist(tb, freq=FALSE, main="r = 4", xlab="replicates of T(n,r) statistic")
points(boot.obj$t0, 0, cex=1, pch=16)




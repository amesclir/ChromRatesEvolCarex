#Setting up the tree
library(ape)
#opening tree
phy <- read.tree("Carex2n.tree")
phy
length(phy$tip.label)
is.ultrametric(phy)
plot(phy, show.tip.label=F)
#it is not ultrametric, see this post
#http://blog.phytools.org/2016/08/fixing-ultrametric-tree-whose-edges-are.html
#library(phytools)
#library(phangorn)
## compute the NNLS ultrametric tree
#phy<-nnls.tree(cophenetic(phy),phy,rooted=TRUE)
## check
#is.ultrametric(phy)
#Now we are making the tree dichotomous
phy <- multi2di(phy)
is.ultrametric(phy)
is.rooted(phy)
#Loading the data
mydata <- read.csv("PPAdata.csv")
mydata
#Pruning the tree
#tips.to.remove <- setdiff(phy$tip.label, mydata$species)
#phy <- drop.tip(phy, tips.to.remove)
#setdiff(phy$tip.label, mydata[,1]) 
#The tree is ready for the analyses!!!!!!
#Setting the character
states <- mydata$R2n 
names(states) = mydata$species
states
#Setting the character variation

sd <- rep(1.49, 755)
names(sd) = mydata$species
sd

#OK. Now the Quasse analyses
library(diversitree)

p <- starting.point.quasse(phy, states)
p
xr <- range(states) + c(-1, 1) * 20 * p["diffusion"]
xr
linear.x <- make.linear.x(xr[1], xr[2])
make <- function(lambda, mu) make.quasse(phy, states, sd, lambda, mu, sampling.f=0.30)
nodrift <- function(f) constrain(f, drift ~ 0)
f.c.c <- make(constant.x, constant.x)
f.l.c <- make(linear.x, constant.x)
f.s.c <- make(sigmoid.x, constant.x) 
f.h.c <- make(noroptimal.x, constant.x) 
 
control <- list(parscale = 0.1, reltol = 0.001)
chrom.mle.c.c <- find.mle(nodrift(f.c.c), p, lower = 0, control = control, verbose = 0)
p.c <- chrom.mle.c.c$par
p.c
p.l.c <- c(p.c[1], l.m = 0, p.c[2:3])
p.l.c
p.s.c <- c(p.c[1], p.c[1], mean(states), 1, p.c[2:3])
p.s.c
p.h.c <- c(p.c[1], p.c[1], mean(states), 1, p.c[2:3])
p.h.c
chrom.mle.d.c.c <- find.mle(f.c.c, coef(chrom.mle.c.c, TRUE), control = control, verbose = 0)
chrom.mle.l.c <- find.mle(nodrift(f.l.c), p.l.c, control = control, verbose = 0)
chrom.mle.d.l.c <- find.mle(f.l.c, coef(chrom.mle.l.c, TRUE), control = control, verbose = 0)
chrom.mle.s.c <- find.mle(nodrift(f.s.c), p.s.c, control = control, verbose = 0)
chrom.mle.d.s.c <- find.mle(f.s.c, coef(chrom.mle.s.c, TRUE), control = control, verbose = 0)
chrom.mle.h.c <- find.mle(nodrift(f.h.c), p.h.c, control = control, verbose = 0)
chrom.mle.d.h.c <- find.mle(f.h.c, coef(chrom.mle.h.c, TRUE), control = control, verbose = 0)
chrom.anova <- anova(chrom.mle.c.c, OU.constant.constant = chrom.mle.d.c.c, BM.linear.constant = chrom.mle.l.c, OU.linear.constant = chrom.mle.d.l.c, BM.sigmoid.constant = chrom.mle.s.c, OU.sigmoid.constant = chrom.mle.d.s.c, BM.hump.constant = chrom.mle.h.c, OU.hump.constant = chrom.mle.d.h.c)
chrom.anova

#Let's explore more models
f.c.l <- make(constant.x, linear.x)
f.c.s <- make(constant.x, sigmoid.x) 
f.c.h <- make(constant.x, noroptimal.x) 
p.c.l <- c(p.c[1:2], m.m = 0, p.c[3]) 
p.c.l
p.c.s <- c(p.c[1:2], p.c[2], mean(states), 1, p.c[3])
p.c.s
p.c.h <- c(p.c[1:2], p.c[2], mean(states), 1, p.c[3])
p.c.h
chrom.mle.c.l <- find.mle(nodrift(f.c.l), p.c.l, control = control, verbose = 0)
chrom.mle.d.c.l <- find.mle(f.c.l, coef(chrom.mle.c.l, TRUE), control = control, verbose = 0)
chrom.mle.c.s <- find.mle(nodrift(f.c.s), p.c.s, control = control, verbose = 0)
chrom.mle.d.c.s <- find.mle(f.c.s, coef(chrom.mle.c.s, TRUE), control = control, verbose = 0)
chrom.mle.c.h <- find.mle(nodrift(f.c.h), p.c.h, control = control, verbose = 0)
chrom.mle.d.c.h <- find.mle(f.c.h, coef(chrom.mle.c.h, TRUE), control = control, verbose = 0)
chrom.anova2 <- anova(chrom.mle.c.c, OU.constant.constant = chrom.mle.d.c.c, BM.linear.constant = chrom.mle.l.c, BM.constant.linear = chrom.mle.c.l, OU.linear.constant = chrom.mle.d.l.c, OU.constant.linear = chrom.mle.d.c.l, BM.sigmoid.constant = chrom.mle.s.c, OU.sigmoid.constant = chrom.mle.d.s.c, BM.constant.sigmoid = chrom.mle.c.s, OU.constant.sigmoid = chrom.mle.d.c.s, BM.hump.constant = chrom.mle.h.c, OU.hump.constant = chrom.mle.d.h.c, BM.constant.hump = chrom.mle.c.h, OU.constant.hump = chrom.mle.d.c.h)
chrom.anova2

#Let's explore more complex models
f.h.l <- make(noroptimal.x, linear.x)
f.h.s <- make(noroptimal.x, sigmoid.x) 
f.l.l <- make(linear.x, linear.x)
f.l.s <- make(linear.x, sigmoid.x) 
f.s.l <- make(sigmoid.x, linear.x)
f.s.s <- make(sigmoid.x, sigmoid.x) 
p.h.l <- c(p.c[1], p.c[1], mean(states), 1,p.c[2], m.m = 0, p.c[3])
p.h.l
p.h.s <- c(p.c[1], p.c[1], mean(states), 1,p.c[2], p.c[2], mean(states), 1, p.c[3])
p.h.s
p.l.l <- c(p.c[1], l.m = 0, p.c[2], m.m = 0, p.c[3])
p.l.l
p.l.s <- c(p.c[1], l.m = 0, p.c[2], p.c[2], mean(states), 1, p.c[3])
p.l.s 
p.s.l <- c(p.c[1], p.c[1], mean(states), 1, p.c[2], m.m = 0, p.c[3])
p.s.l  
p.s.s <- c(p.c[1], p.c[1], mean(states), 1, p.c[2], p.c[2], mean(states), 1, p.c[3])
p.s.s
chrom.mle.h.l <- find.mle(nodrift(f.h.l), p.h.l, control = control, verbose = 0)
chrom.mle.d.h.l <- find.mle(f.h.l, coef(chrom.mle.h.l, TRUE), control = control, verbose = 0)
chrom.mle.h.s <- find.mle(nodrift(f.h.s), p.h.s, control = control, verbose = 0)
chrom.mle.d.h.s <- find.mle(f.h.s, coef(chrom.mle.h.s, TRUE), control = control, verbose = 0)
chrom.mle.l.l <- find.mle(nodrift(f.l.l), p.l.l, control = control, verbose = 0)
chrom.mle.d.l.l <- find.mle(f.l.l, coef(chrom.mle.l.l, TRUE), control = control, verbose = 0)
chrom.mle.l.s <- find.mle(nodrift(f.l.s), p.l.s, control = control, verbose = 0)
chrom.mle.d.l.s <- find.mle(f.l.s, coef(chrom.mle.l.s, TRUE), control = control, verbose = 0)
chrom.mle.s.l <- find.mle(nodrift(f.s.l), p.s.l, control = control, verbose = 0)
chrom.mle.d.s.l <- find.mle(f.s.l, coef(chrom.mle.s.l, TRUE), control = control, verbose = 0)
chrom.mle.s.s <- find.mle(nodrift(f.s.s), p.s.s, control = control, verbose = 0)
chrom.mle.d.s.s <- find.mle(f.s.s, coef(chrom.mle.s.s, TRUE), control = control, verbose = 0)
chrom.anova3 <- anova(chrom.mle.c.c, OU.constant.constant = chrom.mle.d.c.c, BM.linear.constant = chrom.mle.l.c, BM.constant.linear = chrom.mle.c.l, OU.linear.constant = chrom.mle.d.l.c, OU.constant.linear = chrom.mle.d.c.l, BM.sigmoid.constant = chrom.mle.s.c, OU.sigmoid.constant = chrom.mle.d.s.c, BM.constant.sigmoid = chrom.mle.c.s, OU.constant.sigmoid = chrom.mle.d.c.s, BM.hump.constant = chrom.mle.h.c, OU.hump.constant = chrom.mle.d.h.c, BM.constant.hump = chrom.mle.c.h, OU.constant.hump = chrom.mle.d.c.h, BM.hump.linear = chrom.mle.h.l, BM.hump.sigmoid = chrom.mle.h.s, BM.linear.linear = chrom.mle.l.l, BM.linear.sigmoid = chrom.mle.l.s, BM.sigmoid.linear = chrom.mle.s.l, BM.sigmoid.sigmoid = chrom.mle.s.s, OU.hump.linear = chrom.mle.d.h.l, OU.hump.sigmoid = chrom.mle.d.h.s, OU.linear.linear = chrom.mle.d.l.l, OU.linear.sigmoid = chrom.mle.d.l.s, OU.sigmoid.linear = chrom.mle.d.s.l, OU.sigmoid.sigmoid = chrom.mle.d.s.s)
chrom.anova3

save(chrom.mle.c.c, chrom.mle.d.c.c, chrom.mle.l.c, chrom.mle.c.l, chrom.mle.d.l.c, chrom.mle.d.c.l, chrom.mle.s.c, chrom.mle.d.s.c, chrom.mle.c.s, chrom.mle.d.c.s, chrom.mle.h.c, chrom.mle.d.h.c, chrom.mle.c.h, chrom.mle.d.c.h, chrom.mle.h.l, chrom.mle.h.s, chrom.mle.l.l, chrom.mle.l.s, chrom.mle.s.l, chrom.mle.s.s, chrom.mle.d.h.l, chrom.mle.d.h.s, chrom.mle.d.l.l, chrom.mle.d.l.s, chrom.mle.d.s.l, chrom.mle.d.s.s, file = "QuaSSE.RData")

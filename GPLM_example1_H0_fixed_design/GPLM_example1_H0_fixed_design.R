
## clear history
rm(list = ls(all = TRUE))
graphics.off()

## install and load packages
libraries = c("mgcv")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
    install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Epanechnikov kernel function
epkernel = function(x) {
    temp = 0.75 * (1 - x * x)
    temp[x <= -1] = 0
    temp[x >= 1]  = 0
    temp
}

# smoother matrix in Huang and Chen (2008) Epanechnikov kernel normlizing &
# symmatric & sum to one
proj.ep.symmetric = function(x, xgrid, h) {
    gridlength = length(xgrid)
    datalength = length(x)
    gridint    = xgrid[2] - xgrid[1]
    localproj  = matrix(0, ncol = datalength, nrow = datalength)
    bigK       = matrix(0, ncol = gridlength, nrow = datalength)
    for (i in 1:datalength) {
      bigK[i, ] = (epkernel(((x[i] - xgrid)/h))/h)
    }
    adjvect = as.numeric(1/bigK %*% rep(gridint, gridlength))
    bigKK   = (diag(adjvect, nrow = datalength, ncol = datalength) %*% bigK)
    
    for (i in 1:gridlength) {
        kweight   = diag(bigKK[, i], ncol = datalength)
        bigX      = cbind(rep(1, length(x)), (x - xgrid[i]))
        localH    = bigX %*% (solve(t(bigX) %*% kweight %*% bigX)) %*% t(bigX) %*% kweight
        localproj = localproj + kweight %*% localH
    }
    localproj  = localproj * gridint
    defree     = sum(diag(localproj))
    list(Hstar = localproj, defree = defree)
}
## analysis of deviance function 
## using Epanechnikov kernel, xgrid must be equally spaced testing no effect of x
anodev.logit.chisq = function(y, x, xgrid, h) {
    gridlength = length(xgrid)
    gridint    = xgrid[2] - xgrid[1]
    datalength = length(x)
    bigK       = matrix(0, ncol = gridlength, nrow = datalength)
    for (i in 1:datalength) {
      bigK[i, ] = (epkernel(((x[i] - xgrid)/h))/h)
    }
    adjvect = as.numeric(1/bigK %*% rep(gridint, gridlength))
    bigKK   = (diag(adjvect, nrow = datalength, ncol = datalength) %*% bigK)
    
    global.glm = glm(y ~ 1, family = binomial(link = "logit"))
    
    d.dev   = rep(0, gridlength)
    d.coeff = matrix(0, nrow = gridlength, ncol = 2)
    
    for (i in 1:gridlength) {
        x.local  = x - xgrid[i]
        kweights = bigKK[, i]
        
        d.glm        = glm(y ~ x.local, family = binomial(link = "logit"), weights = kweights)
        d.dev[i]     = d.glm$dev
        d.coeff[i, ] = d.glm$coeff
        
    }
    # i loop
    intD      = sum(d.dev) * gridint
    teststat  = -(intD - global.glm$dev)
    list(intD = intD, teststat = teststat)
    
}

# a=0 under H0
a = 0

# sample size 50, 100, 200
samplen = 50
# number of simulations
Nsim = 5000
# fixed equally-spaced design
xt = seq(0, 1, length = samplen)

# try 7 values of bandwidth, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, and 0.3
hlength = 7

# initialization
chi.teststat = matrix(0, nrow = Nsim, ncol = hlength)
pvalue       = matrix(0, nrow = Nsim, ncol = hlength)
degfee       = matrix(0, nrow = Nsim, ncol = hlength)

# integrated deviance
intDall = matrix(0, nrow = Nsim, ncol = hlength)

# AICc values for 7 values of bandwidth
AICc = matrix(0, nrow = Nsim, ncol = hlength)

# bandwidth select by AICc and its p-value
AICch      = rep(0, Nsim)
AICcpvalue = rep(0, Nsim)

# Horowitz and Spokoiny (2001) set the initial value as the smallest bandwidth
adaptiveh      = rep(0.1, Nsim)
adaptivepvalue = rep(0.99999, Nsim)

# for gam in mgcv
gam.p0  = rep(0, Nsim)
gam.p1  = rep(0, Nsim)
gam.pn  = rep(0, Nsim)
gam.edf = rep(0, Nsim)
gam.ts  = rep(0, Nsim)


# a=0 n=50 fixed design xtequally-spaced
set.seed(2015012710)

# a=0 n=100 fixed design x equally-spaced set.seed(2012061310)

# a=0 n=200 fixed design x equally-spaced set.seed(2015012810)

tgrid = seq(0, 1, 0.005)  #length 201

for (j in 1:Nsim) {
    beta0 = -1
    # under H0
    eta       = beta0 + a * cos(2 * pi * xt)
    p0        = exp(eta)/(1 + exp(eta))
    y0        = as.integer(rbinom(samplen, 1, p0))
    t.order   = order(xt)
    d4        = data.frame(xt[t.order], y0[t.order])
    names(d4) = c("xt", "y0")
    hchoice   = c(0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3)
    AICcmin   = 10
    for (hi in 1:hlength) {
        h = hchoice[hi]
        
        # call the function
        chi.teststat[j, hi] = anodev.logit.chisq(d4$y0, d4$xt, tgrid, h)$teststat
        intDall[j, hi]      = anodev.logit.chisq(d4$y0, d4$xt, tgrid, h)$intD
        
        # get the smoother matrix
        d4.H = proj.ep.symmetric(d4$xt, tgrid, h)
        
        # get the trace of smoother matrix
        degfee[j, hi] = d4.H$defree
        
        # calculate AICc and find AICcmin
        AICc[j, hi] = log(intDall[j, hi]/samplen) + 2 * (degfee[j, hi] + 1)/(samplen - 
            degfee[j, hi] - 2)
        if (AICc[j, hi] < AICcmin) {
            AICcmin  = AICc[j, hi]
            AICch[j] = h
        }
        # proposed chi-square test statistic p-value
        pvalue[j, hi] = 1 - pchisq(chi.teststat[j, hi], (degfee[j, hi] - 1))
        
        # find the smallest pvalue based on Horowitz and Spokoiny (2001)
        if (pvalue[j, hi] < adaptivepvalue[j]) {
            adaptivepvalue[j] = pvalue[j, hi]
            adaptiveh[j]      = h
        }
        
    }
    # end hi bandiwdth
    findaicc      = (hchoice == AICch[j])
    AICcpvalue[j] = pvalue[j, findaicc]
    
    # mgcv package
    b.gam      = gam(d4$y0 ~ s(d4$xt), method = "REML", family = binomial(link = "logit"))
    gam.p0[j]  = summary(b.gam, p.type = 0)$s.pv
    gam.p1[j]  = summary(b.gam, p.type = 1)$s.pv
    gam.pn[j]  = summary(b.gam, p.type = -1)$s.pv
    gam.edf[j] = summary(b.gam)$edf
    gam.ts[j]  = summary(b.gam)$chi.sq
    
    cat(j, "/")
    
}
# j loop


testresults1 = round(matrix(c(hchoice, colSums(pvalue <= 0.05)/Nsim), nrow = 2, ncol = 7, 
    byrow = TRUE), 6)
testresults1
testresults2 = round(matrix(c(mean(AICch), sd(AICch), sum(AICcpvalue < 0.05)/Nsim, 
    mean(adaptiveh), sd(adaptiveh), sum(adaptivepvalue < 0.05)/Nsim), nrow = 2, byrow = TRUE), 
    6)
testresults2

degresults = round(matrix(c(mean(degfee[, 1]), sd(degfee[, 1]), mean(degfee[, 2]), 
    sd(degfee[, 2]), mean(degfee[, 3]), sd(degfee[, 3]), mean(degfee[, 4]), sd(degfee[, 
        4]), mean(degfee[, 5]), sd(degfee[, 5]), mean(degfee[, 6]), sd(degfee[, 6]), 
    mean(degfee[, 7]), sd(degfee[, 7])), nrow = 2, byrow = FALSE), 6)
degresults

print(table(AICch))
print(table(adaptiveh))

print("gam")
print(sum(gam.p0 <= 0.05))/Nsim
print(sum(gam.p1 <= 0.05))/Nsim
print(sum(gam.pn <= 0.05))
print(summary(gam.edf))
print(sqrt(var(gam.edf)))

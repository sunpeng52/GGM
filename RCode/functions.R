banded.m <- function(n,b) {
        ## generate banded matrix 
        diff.ij <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                               (1:n - 1))
        (1-diff.ij/b)*(1 > (diff.ij/b))
}


block.m <- function(p,cluster,rho){
        
        index = 1:p
        
        id.max = seq(p/cluster,p,p/cluster)
        id.max[-length(id.max)]
        
        base = matrix(rho,p/cluster,p/cluster)+diag(1-rho,p/cluster)
        ressss = diag(cluster) %x% base
        
        for(i in id.max[-length(id.max)]){
                for(j in i+1:(p/cluster)){
                        ressss[i,j] = rho
                        ressss[j,i] = rho
                }
        }
        
        return(ressss)
        
}

pooled.S <- function(S, N){
        ## compute pooled sample covariance
        m <- length(N)
        p <- dim(S[[1]])[1]
        temp <- matrix(0,p,p)
        for(i in 1:m)
                temp <- temp + (N[i]-1)*S[[i]]
        return(temp/(sum(N)-m))
}

SoftThreshold <- function(x, lam) {
        # note: this works also if lam is a matrix of the same size as x.
        sign(x) * (abs(x) - lam) * (abs(x) > lam)
}

LogDet <- function(M){
        # note: determinant() is better than log(det()) when det is very close to zero
        log.det <- determinant(M, logarithm=TRUE)
        if (log.det$sign == -1)
                return(NA)
        as.numeric(log.det$mod)
        #log(det(M))
}


ComputeLikelihood <- function(R, G, N, S, mu){
        m <- length(N)
        p <- nrow(R)
        ind.diag <- 1 + 0:(p-1)*(p+1)
        S.bar <- lapply(mu, function(x) matrix(x,p)%*%t(matrix(x,p)))
        res <- 0
        for(i in 1:m)
                res <- res - (1/2) * (LogDet(R+N[i]*G) + (N[i]-1)*LogDet(R) +
                                              sum(solve(R, (N[i]-1)*S[[i]])[ind.diag]) +
                                              sum(solve(R+N[i]*G,N[i]*S.bar[[i]])[ind.diag]))
        return(res)
}

ComputeBIC <- function(R, G, N, S, mu, p){
        
        -2 * ComputeLikelihood(R, G, N, S, mu) +
                log(length(N))*(sum(G!=0)+p)/2+log(sum(N))*(sum(R!=0)+p)/2
        
}

ComputePenalty <- function(R, G, lambda.r, lambda.g){
        sum(lambda.r*abs(R)) + sum(lambda.g*abs(G))
}


ComputeObjective <- function(R, G, N, S, mu, lambda.r, lambda.g) {
        # the original non-convex problem's objective function
        -2 * ComputeLikelihood(R, G, N, S, mu) + ComputePenalty(R, G, lambda.r, lambda.g)
}

ComputeGradientOfgg <- function(R, G, G0, N, mu){
        m <- length(N)
        p <- nrow(R)
        S.bar <- lapply(mu, function(x,p) matrix(x,p)%*%t(matrix(x,p)),p=p)
        res <- 0
        for(i in 1:m)
                res <- res +N[i]*(solve(R+N[i]*G0) 
                                  - solve(R+N[i]*G)%*%(N[i]*S.bar[[i]])%*%solve(R+N[i]*G))
        return(res)
}

ComputeGradientOfgr <- function(R, R0, G, N, S, mu){
        m <- length(N)
        p <- nrow(R)
        S.bar <- lapply(mu, function(x,p) matrix(x,p)%*%t(matrix(x,p)),p=p)
        res <- 0
        for(i in 1:m)
                res <- res + (solve(R0+N[i]*G) 
                              - solve(R+N[i]*G)%*%(N[i]*S.bar[[i]])%*%solve(R+N[i]*G)) + (N[i]-1)*(
                                      solve(R0)-solve(R)%*%S[[i]]%*%solve(R)
                                  )
        return(res)
}


GGDescent <- function(R, R0, G, G0, N, S, mu, lambda.r, lambda.g, del.r, del.g, 
                      nsteps, step.size, tol.r = 1e-3, tol.g = 1e-3){
        #  tol:  convergence threshold.  Stops when mean(abs(Sigma-Sigma.last)) < tol
        
        ttt <- step.size
        
        converged.g <- FALSE
        converged.r <- FALSE
        
        
        ## update G
        
        obj.g <- ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g)
        G.last <- G
        for (i in seq(nsteps)) {
                
                grad.gg <- ComputeGradientOfgg(R, G, G0, N, mu)
                grad.gg <- (grad.gg + t(grad.gg)) / 2 # make sure this stays symmetric
                
                
                G <- ProxADMM(G - ttt * grad.gg, del.g, 1, P=lambda.g*ttt, rho=.1)$X
                # check for convergence:
                if (mean(abs(G - G.last)) < tol.g)
                        converged.g <- TRUE
                
                obj.g <- c(obj.g, ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g))
                
                if (obj.g[i+1] > obj.g[i]){
                        G <- G.last
                        break
                }else{
                        G.last <- G
                }
                
                if(converged.g)
                        break
        }
        
        
        ## update R
        
        obj.r <- ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g)
        R.last <- R
        
        for (i in seq(nsteps)) {
                
                grad.gr <- ComputeGradientOfgr(R, R0, G, N, S, mu)
                grad.gr <- (grad.gr + t(grad.gr)) / 2 # make sure this stays symmetric
                
                R <- ProxADMM(R - ttt * grad.gr, del.r, 1, P=lambda.r*ttt, rho=.1)$X
                # check for convergence:
                if (mean(abs(R - R.last)) < tol.r)
                        converged.r <- TRUE
                
                obj.r <- c(obj.r, ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g))
                
                if (obj.r[i+1] > obj.r[i]){
                        R <- R.last
                        break
                }else{
                        R.last <- R  
                }
                
                if (converged.r)
                        break
        }
        
        list(R=R,G=G)
}

GGDescent2 <- function(R, R0, G, G0, N, S, mu, lambda.r, lambda.g, del.r, del.g, 
                      nsteps, step.size, tol.r = 1e-3, tol.g = 1e-3){
        #  tol:  convergence threshold.  Stops when mean(abs(Sigma-Sigma.last)) < tol
        
        ttt <- step.size
        
        converged.g <- FALSE
        converged.r <- FALSE
        
        
        ## update G
        
        obj.g <- ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g)
        G.last <- G
        for (i in seq(nsteps)) {
                
                grad.gg <- ComputeGradientOfgg(R, G, G0, N, mu)
                grad.gg <- (grad.gg + t(grad.gg)) / 2 # make sure this stays symmetric
                
                
                G <- ProxADMM2(G - ttt * grad.gg, del.g, 1, P=lambda.g*ttt, rho=.1)$Z
                
                
                # check for convergence:
                if (mean(abs(G - G.last)) < tol.g)
                        converged.g <- TRUE
                
                obj.g <- c(obj.g, ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g))
                
                if (obj.g[i+1] > obj.g[i]){
                        G <- G.last
                        break
                }else{
                        G.last <- G
                }
                
                if(converged.g)
                        break
        }
        
        
        ## update R
        
        obj.r <- ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g)
        R.last <- R
        
        for (i in seq(nsteps)) {
                
                grad.gr <- ComputeGradientOfgr(R, R0, G, N, S, mu)
                grad.gr <- (grad.gr + t(grad.gr)) / 2 # make sure this stays symmetric
                
                R <- ProxADMM2(R - ttt * grad.gr, del.r, 1, P=lambda.r*ttt, rho=.1)$Z
                # check for convergence:
                if (mean(abs(R - R.last)) < tol.r)
                        converged.r <- TRUE
                
                obj.r <- c(obj.r, ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g))
                
                if (obj.r[i+1] > obj.r[i]){
                        R <- R.last
                        break
                }else{
                        R.last <- R  
                }
                
                if (converged.r)
                        break
        }
        
        list(R=R,G=G)
}


pooled.S <- function(S, N){
        m <- length(N)
        p <- dim(S[[1]])[1]
        temp <- matrix(0,p,p)
        for(i in 1:m)
                temp <- temp + (N[i]-1)*S[[i]]
        return(temp/(sum(N)-m))
}

ProxADMM2 <- function(A, del, lam, P, rho=.1, tol=1e-6, maxiters=100, verb=FALSE) {
        # Minimize_X { (1/2)||X - A||_F^2 + lam||P*X||_1} s.t. X >= del * I
        #
        # ADMM approach
        
        # first, check if simple soft-thesholding works... if so, skip the ADMM!
        soft <- SoftThreshold(A, lam * P)
        minev <- min(eigen(soft, symmetric=T, only.values=T)$val)
        if (minev >= del) {
                return(list(X=soft, Z=soft, obj=ComputeProxObjective(soft, A, lam, P)))
        }
        
        p <- nrow(A)
        obj <- NULL
        
        # initialize X, Y
        X <- soft
        Y <- matrix(0, p, p)
        
        # main loop
        for (i in seq(maxiters)) {    
                # update Z:
                B <- (X + (1/rho)*Y)
                if (min(eigen(B, symmetric=T, only.values=T)$val) < del) {
                        # note: even though eigen is called twice, only.values=T is
                        #       much faster, making this worthwhile.
                        eig <- eigen(B, symmetric=T)
                        Z <- eig$vec %*% diag(pmax(eig$val, del)) %*% t(eig$vec)
                }
                else {
                        Z <- B
                }
                
                
                # update X:
                X <- (rho/(1+rho))*SoftThreshold((A-Y)/rho + Z, lam * P / rho)
                
                # check for convergence:
                obj <- c(obj, ComputeProxObjective(X, A, lam, P))
                if (verb)
                        cat(" ", obj[i], fill=T)
                if (i > 1)
                        if (obj[i] > obj[i - 1] - tol) {
                                if (verb)
                                        cat(" ADMM converged after ", i, " steps.", fill=T)
                                break
                        }
                
                # update Y:
                Y <- Y + rho * (X - Z)
        }
        
        list(X=X, Z=Z, obj=obj)
}

ComputeProxObjective <- function(X, A, lam, P) {
        sum((X-A)^2) / 2 + lam * sum(abs(P*X))
}

group.spcov <- function(R, G, N, S, mu, lambda.r, lambda.g, step.size, 
                        del.r, del.g, 
                        n.outer.steps = 1e4, n.inner.steps = 1e4, 
                        tol.outer = 1e-4, thr.inner.r = 1e-2, thr.inner.g = 1e-2){
        
        objective <- ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g)
        for (i in seq(n.outer.steps)){
                R0 <- R
                G0 <- G
                gg.temp <- GGDescent(R=R, R0=R0, G=G, G0=G0, 
                                     N=N, S=S, mu=mu, 
                                     lambda.r=lambda.r, lambda.g=lambda.g, 
                                     del.r=del.r, del.g=del.g, 
                                     nsteps=n.inner.steps, 
                                     step.size=step.size, 
                                     tol.r=thr.inner.r, tol.g = thr.inner.g)
                R <- gg.temp$R
                G <- gg.temp$G
                objective <- c(objective, ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g))
                
                if((objective[i + 1] > objective[i] - tol.outer) ){
                        break
                }
        }
        list(R=gg.temp$R, G=gg.temp$G)
}

group.spcov2 <- function(R, G, N, S, mu, lambda.r, lambda.g, step.size, 
                        del.r, del.g, 
                        n.outer.steps = 1e4, n.inner.steps = 1e4, 
                        tol.outer = 1e-4, thr.inner.r = 1e-2, thr.inner.g = 1e-2){
        
        objective <- ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g)
        for (i in seq(n.outer.steps)){
                R0 <- R
                G0 <- G
                gg.temp <- GGDescent2(R=R, R0=R0, G=G, G0=G0, 
                                     N=N, S=S, mu=mu, 
                                     lambda.r=lambda.r, lambda.g=lambda.g, 
                                     del.r=del.r, del.g=del.g, 
                                     nsteps=n.inner.steps, 
                                     step.size=step.size, 
                                     tol.r=thr.inner.r, tol.g = thr.inner.g)
                R <- gg.temp$R
                G <- gg.temp$G
                objective <- c(objective, ComputeObjective(R, G, N, S, mu, lambda.r, lambda.g))
                
                if((objective[i + 1] > objective[i] - tol.outer) ){
                        break
                }
        }
        list(R=gg.temp$R, G=gg.temp$G)
}

gg <- function(G,G0,R,N,mu,lambda.g){
        m <- length(N)
        p <- nrow(R)
        S.bar <- lapply(mu, function(x,p) matrix(x,p)%*%t(matrix(x,p)),p=p)
        res <- 0
        for(i in 1:m)
                res <- res + (solve(R0+N[i]*G) 
                              - solve(R+N[i]*G)%*%(N[i]*S.bar[[i]])%*%solve(R+N[i]*G)) + (N[i]-1)*(
                                      solve(R0)-solve(R)%*%S[[i]]%*%solve(R)
                              )
        return(res)
}

admm.cov <- function (A, del, lam, P, rho = 0.1, tol.abs = 1e-8, tol.rel = 1e-8, maxiters = 1000, 
                      vary = TRUE, verb = FALSE) {
        
        ## initialization
        soft <- SoftThreshold(A, lam * P)
        min.eigen <- min(eigen(soft, symmetric = T, only.values = T)$values)
        if (min.eigen >= del) {
                return(list(Sigma = soft, Theta = soft, 
                            obj = ComputeProxObjective(soft, A, lam, P)))
        }
        
        p <- nrow(A)
        
        obj <- NULL
        Sigma.path <- NULL
        Theta.path <- NULL
        Lambda.path <- NULL
        
        Theta <- soft
        Lambda <- matrix(0, p, p)
        ## Z <- soft
        ## Y <- matrix(0, p, p)
        for (i in seq(maxiters)) {
                
                ## Sigma update:
                #B <- (A + rho * Z - Y)/(1 + rho)
                Temp <- (A + rho * Theta - Lambda)/(1+rho)
                #if (min(eigen(B, symmetric = T, only.values = T)$val) < 
                #    del) {
                if (min(eigen(Temp, symmetric = T, only.values = T)$values) < del){
                        eigen.sys <- eigen(Temp, symmetric = T)
                        Sigma <- eigen.sys$vectors %*% diag(pmax(eigen.sys$values, del)) %*% t(eigen.sys$vectors)
                }
                else {
                        Sigma <- Temp
                }
                
                obj <- c(obj, ComputeProxObjective(Sigma, A, lam, P))
                Sigma.path <- c(Sigma.path, Sigma)
                
                ## Theta update:
                Theta <- SoftThreshold(Sigma + Lambda/rho, lam * P/rho)
                Theta.path <- c(Theta.path, Theta)
                
                ## Lambda update:
                Lambda <- Lambda + rho * (Sigma - Theta)
                Lambda.path <- c(Lambda.path, Lambda)
                
                if (verb) 
                        cat(" ", obj[i], fill = T)
                if(i > 1){
                        tol.pri <- p * tol.abs + tol.rel * max(sqrt(sum(Sigma^2)), sqrt(sum(Theta^2)))
                        tol.dual <- p * tol.abs + tol.rel * sqrt(sum(Lambda^2))
                        
                        r.norm <- sqrt(sum((Sigma - Theta)^2))
                        s.norm <- sqrt(sum((rho * (Theta.path[i] - Theta.path[i-1]))^2))
                        if(r.norm <= tol.pri & s.norm <= tol.dual){
                                if (verb) 
                                        cat(" ADMM converged after ", i, " steps.", 
                                            fill = T)
                                break
                        }
                }
                
                if(vary & i > 1){
                        if(r.norm > 10 * s.norm){
                                rho <- 2 * rho
                        }else if(s.norm > 10 * r.norm){
                                rho <- rho/2
                        }else{
                                rho <- rho
                        }
                }
                
        }
        list(Sigma = Sigma, Theta = Theta, obj = obj)
}


# auxiliary functions
# author: sebastian daza


# simulation functions
# scaling model 
simScaling = function(E, a0 = 0, a1 = 0.5, b0 = 0.8, b1 = 0.2, h2 = 0.4) {
    N = length(E)
    G = rnorm(N,0,1)
    h = sqrt(h2)
    e = sqrt(1-h^2)
    sigma = exp(b0*e + b1*e*E) 
    y = rnorm(N, a0 + a1*E + b0*h*G + b1*h*E*G, sigma)
    df = data.frame(E = E, y = y, g = G, ys = scale(y))
}

# interaction model
simInteraction = function(E, a0 = 0.0, a1 = 0.8, a2 = 0.5, a3 = 0.3, bs0 = 0.4, bs1 = 0.24) {
    N = length(E)
    G = rnorm(N,0,1)
    sigma = exp(bs0 + bs1 * E) 
    y = rnorm(N, a0 + a1*E + a2*G + a3*G*E, sigma)
    df = data.frame(E = E, y = y, g = G, ys = scale(y))
}

# interaction + scaling model
simScalingInteraction = function(E, a0 = 0.0,  a1 = 0.5, a2 = 0.10, b0 = 0.8, b1 = 0.2, h2 = 0.5) {
    N = length(E)
    G = rnorm(N,0,1)
    h = sqrt(h2)
    e = sqrt(1-h^2)
    sigma = exp(b0*e + b1*e*E) 
    y = rnorm(N, a0 + a1*E + a2*G*E + b0*h*G + b1*h*E*G, sigma)
    df = data.frame(E = E, y = y, g = G, ys = scale(y))
}

# domingue's original function
simDom = function(E,b0 = .8, b1 = .2, b2 = 0, b3 = .05, h = sqrt(.6), a =.5, sigma = 1, scaling = TRUE) {
    N = length(E)
    G = rnorm(N,0,1)
    eps = rnorm(N,0,sigma)
    if (scaling){
        e = sqrt(1-h^2)
        ystar = h*G+e*eps
        y = a*E+(b0 + b1*E)*ystar
    } else {
        y = b1*G+b2*E+b3*G*E+eps
    }
    df = data.frame(E=E,y=y,g=G, ys = scale(y))
}

# sim alternative
simAlternative = function(E, a0 = 0.0, a1 = 0.8, bs0 = 0.50, bs1 = 0.35, 
    h2 = 0.3, e2 = 0.2, he2 = 0.15) {
    
    n = length(E)
    G = rnorm(n)

    hd2 = 1 - (h2 + e2 + he2)

    h = sqrt(h2)
    e = sqrt(e2)
    he = sqrt(he2)
    hd = sqrt(hd2)
    
    ys = h * G + e * E + he * G * E

    sigma = exp(bs0 * hd  + bs1 * hd * E)
    y = rnorm(n, a0 + a1 * ys, sigma)
    #y = a + (b0 + b1 * E) * ys + eps
    df = data.frame(E = E, y = y, g = G, ys = scale(y))
}


# decomposition function scaliing vs  interaction model 
decompR = function(model ,p0 = "g", p1 = "g:e", l0 = "sigma_intercept",
      l1 = "sigma_e") {

    par = list("p0" = p0, "p1" =  p1, "l0" = l0, "l1" =  l1)
    par = lapply(par, function(x) tolower(paste0("b_", x)))
    s = posterior_samples(model)
    names(s) = tolower(names(s))
    sl0 = s[[par$l0]]
    sl1 = s[[par$l1]]
    sp0 = s[[par$p0]]
    sp1 = s[[par$p1]]
    v =  ( (sl0*sp1 - sl1*sp0)/sl0 ) / sp1
    v
}


# scaling test
scalingTest = function(model, 
    h0 ="g * sigma_E = g:E * sigma_Intercept", alpha = 0.05) {
    h = hypothesis(model, h0, alpha = alpha)
    return(h)
}


plotScalingTest = function(test, file = NULL) {
    v = test$samples
    m = mean(v$H1)
    Q2.5 = quantile(v$H1, 0.025)
    Q97.5 = quantile(v$H1, 0.975)
    subtitle = paste0("Mean = ", round(m, 3), " [", round(Q2.5, 3), ", ", round(Q97.5, 3), "]")    
    p =   ggplot(v, aes(x = H1))  +
        geom_density(alpha = 0.4, color = "#2c7fb8", fill = "#2c7fb8") + theme_minimal() + 
        labs(title = TeX("Posterior distribution $\\pi_0\\lambda_1 - \\pi_1\\lambda_0$"), 
        subtitle = subtitle,
        x = " ", y = "Density")

    if (!is.null(file)) {
        savepdf(file)
        print(p)
        dev.off()
    }
}


plotDecomp = function(model, file = NULL) {
    x = decompR(model)
    m = mean(x)
    Q2.5 = quantile(x, 0.025)
    Q97.5 = quantile(x, 0.975)
    subtitle = paste0("Mean = ", round(m, 3), " [", round(Q2.5, 3), ", ", round(Q97.5, 3), "]")    
    p = ggplot(data.frame(x = x), aes(x = x)) + 
        geom_density(alpha = 0.4, color = "#2c7fb8", fill = "#2c7fb8") + theme_minimal() + 
        labs(title = "Posterior distribution of proportion departuring from the scaling model", 
            subtitle = subtitle, x = " ", y = "Density")

    if (!is.null(file)) {
        savepdf(file)
        print(p)
        dev.off()
    }

}


plot_multi_histogram <- function(df, feature, label_column) {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text = label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", 
        linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
    plt + guides(fill=guide_legend(title=label_column))
}


savepdf = function(file, width = 16, height = 10, mgp = c(2.2,0.45,0),
    tcl = -0.4, mar = c(3.3,3.6,1.1,1.1)) {
    fname = paste0(file, ".pdf")
    pdf(fname, width=width / 2.54, height = height / 2.54,
        pointsize = 10)
    par(mgp = mgp, tcl = tcl, mar = mar)
}


specify_decimal = function(x, k) trimws(format(round(x, k), nsmall=k))


countmis = function(dat, vars = NULL, pct = TRUE, exclude.complete = TRUE) {

    if (is.null(vars)) {
        vars = names(dat)
    }

    mis = sort(sapply(dat[, vars, with = FALSE],
        function(x) sum(is.na(x))), decreasing = TRUE)

    if (exclude.complete == TRUE) {
         mis = mis[mis > 0]
    }

    if (pct == FALSE)
        { return(mis) }
    else if ( pct == TRUE ) {
        return( round(mis / nrow(dat), 3))
    }
    return(mis)
}


bayes_r2 = function(y, ypred) {
    e = -1 * sweep(ypred, 2, y)
    var_ypred = matrixStats::rowVars(ypred)
    var_e = matrixStats::rowVars(e)
    r2 = as.matrix(var_ypred / (var_ypred + var_e))
    return(c("m" = median(r2), "l" = as.numeric(quantile(r2, probs = 0.025)), 
        "h" = as.numeric(quantile(r2, probs = 0.975))))
}

bayes_r2_group = function(y, ypred, group) {
    ugroup = sort(unique(group))
    e = -1 * sweep(ypred, 2, y)

    v = NULL
    for (i in ugroup) {
        f = which(group %in% i)
        var_e = matrixStats::rowVars(e[, f])
        var_ypred = matrixStats::rowVars(ypred[, f])
        r2 = as.matrix(var_ypred / (var_ypred + var_e))
        v = rbind(v, c("group" = i, "m" = median(r2), "l" = as.numeric(quantile(r2, probs = 0.025)), 
            "h" = as.numeric(quantile(r2, probs = 0.975))))
    }
    return(data.table(v))
}

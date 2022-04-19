# auxiliary functions
# author: sebastian daza

library(haven)
library(latex2exp)
library(data.table)
library(brms)
library(tidybayes)
library(texreg)
library(ggplot2)
library(scalingGxE)
library(patchwork)
library(foreach)
library(doParallel)


# simulation functions

# scaling model 
simScaling = function(E, a0 = 0, a1 = 0.5, b0 = 0.8, b1 = 0.2, h2 = 0.4) {
    N = length(E)
    G = rnorm(N,0,1)
    h = sqrt(h2)
    e = sqrt(1-h^2)
    sigma = exp(b0*e + b1*e*E) 
    y = rnorm(N, a0 + a1*E + b0*h*G + b1*h*E*G, sigma)
    df = data.table(E = E, y = y, g = G)
    df[, qE := cut(E, quantile(E, probs = 0:10/10),
        labels = FALSE, include.lowest = TRUE)]
    df[, zy := scale(y)][, zg := scale(g)]
    df[, gzy := scale(y), qE]
    df[, gzg := scale(g), qE]
}


 # interaction, heteroscedasticity
simHC = function(E) {
    N = length(E)
    G = rnorm(N, 0, 1)
    sigma = exp(0.2 + 0.5 * E) 
    y = rnorm(N, E * 0.4 +  0.5 * G +  0.2 * E * G, sigma)
    df = data.table(E = E, y = y, g = G)
    df[, qE := cut(E, quantile(E, probs = 0:10/10),
        labels = FALSE, include.lowest = TRUE)]
    df[, zy := scale(y)][, zg := scale(g)]
    df[, gzy := scale(y), qE]
    df[, gzg := scale(g), qE]
}


# interaction model, heteroscedasticity
simInteraction = function(E, a0 = 0.0, a1 = 0.8, a2 = 0.5, a3 = 0.3, bs0 = 0.4, bs1 = 0.24) {
    N = length(E)
    G = rnorm(N,0,1)
    sigma = exp(bs0 + bs1 * E) 
    y = rnorm(N, a0 + a1 * E + a2 * G + a3 * G * E, sigma)
    df = data.table(E = E, y = y, g = G)
    df[, qE := cut(E, quantile(E, probs = 0:10/10),
        labels = FALSE, include.lowest = TRUE)]
    df[, zy := scale(y)][, zg := scale(g)]
    df[, gzy := scale(y), qE]
    df[, gzg := scale(g), qE]
}

# test functions

# scaling test
scalingTest = function(model, 
    h0 ="g * sigma_E = g:E * sigma_Intercept", alpha = 0.05) {
    h = hypothesis(model, h0, alpha = alpha)
    return(h)
}


scalingTestC = function(p0, p1, l0, l1) {
    return(p0 * l1 - p1 * l0)
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
    } else {
        print(p)
    }
}

# decomposition function scaling vs interaction model 
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


# plot decomposition function
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
    } else {
        print(p)
    }

}


slopeCorPlot = function(data, group, file = NULL, ncol = 3, w = 25, h = 20) {
    plots = list()
    v = sort(unique(data[[group]]))
    print(v)
    for (i in v) {
    slope = specify_decimal(coef(lm(y ~ g, data = data[get(group) == i]))[2], 2)
    corr = specify_decimal(cor(data[get(group) == i, .(y, g)])[1, 2], 2)
    plots[[i]] = ggplot(data[get(group) == i], aes(g, y)) + geom_point(size = 0.5, 
        alpha = 0.1, color = "#2b8cbe") +
        geom_smooth(method = "lm", color = "#e34a33", alpha = 0.2, size = 0.3) + 
        labs(title = paste0("E", i), 
            subtitle = TeX(paste0("$\\rho$=", corr, ", $\\beta$=", slope))) +
        theme_minimal()
    }
    if (!is.null(file)) {
        savepdf(file, w, h)
        print(wrap_plots(plots, ncol = 3))
        dev.off()
    }
    else { return(plots) }
}


# multiclass function to save ouptut
multiResultClass = function(distributional = NULL, slope = NULL, correlation = NULL) {
    me = list(distributional = distributional, slope = slope, correlation = correlation)
    class(me) = append(class(me), "multiResultClass")
    return(me)
}


testConvergence = function(models) {
    rhats = NULL
    for (i in seq_along(models)) { 
        rhats  = c(rhats, as.vector(unlist(models[[i]]$rhats)))
    }
    print(sum(rhats > 1.05 | rhats < 0.95))
}


# custom function to run models using clusters
runModels = function(data, clusters = 10) {
    
    cl = makeCluster(clusters)
    registerDoParallel(cl)
    n = length(data)
    output = foreach(i = 1:n, 
        .packages = c("brms", "data.table")) %dopar% {

        source("src/utils.R")
        results = multiResultClass()  
        temp = data[[i]]
        f = bf(y ~ g + E + g * E, sigma ~ 1 + E)
        a = brm(f, 
            data = temp, 
            family = brmsfamily("gaussian", link_sigma = "log"), 
            chains = 1, backend = "cmdstanr", cores = 1,
            control = list(adapt_delta = 0.99))
        f = bf(zy ~ zg + (zg|qE))
        b = brm(f, 
            data = temp, 
            chains = 1, backend = "cmdstanr", cores = 1, 
            control = list(adapt_delta = 0.99))
        f = bf(gzy ~ gzg + (gzg|qE))
        c = brm(f, 
            data = temp, 
            chains = 1, backend = "cmdstanr", cores = 1, 
            control = list(adapt_delta = 0.99))

        results$distributional = a
        results$slope = b
        results$correlation = c
        return(results)
    }
    stopCluster(cl)

    a = combine_models(mlist = extractModelList(output, "distributional"), check_data = FALSE)
    b = combine_models(mlist = extractModelList(output, "slope"), check_data = FALSE)
    c = combine_models(mlist = extractModelList(output, "correlation"), check_data = FALSE)
    rm(output)
    return(list("distributional" = a, "slope" = b, "correlation" = c))
}


# extract objects from list
extractModelList = function(list, name, threshold = 1.1) {
    output = list()
    n = length(list)
    for (i in 1:n) {
        m = list[[i]][[name]]
        if (sum(na.omit(rhat(m)) > threshold )) {
            print(paste0("Warning: rhat > ", threshold))
        }   
        output[[i]] = m
    }
    return(output)
    
}


# plot varying coefficients
brmsVaryingCoefPlot = function(model, coef, dimension, ci = 0.95, file = NULL, 
    title = NULL, subtitle = NULL, x = "x", y = "y", caption = NULL, 
    angle = 0.9, hjust = 0.5, vjust = 0.5,
    yintercept = 0, return_data = FALSE) {


    scoef = sapply(c(coef, dimension), str2lang)
    ccoef = gsub("b_", "", coef)
    sdim = gsub("\\[.+", "", dimension)
    cgroup = gsub("(.+\\[)([a-zA-Z0-9]+)(,.+)", "\\2", dimension)
    ff = str2lang(paste0("median = ", coef, " + ", sdim))
    t = tidybayes::spread_draws(model, !!!scoef)
    s = tidybayes::median_qi(t, !!ff, .width = ci)
    setnames(s, paste0("median = ", coef, " + ", sdim), "median")

    p = ggplot(s, aes_string(cgroup, "median", group = 1)) + 
            geom_line(color='#2b8cbe', size = 0.4)  +
            geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = '#a6bddb', alpha=0.2)  +
                    labs(title = title, subtitle = subtitle, x = x, y = y, caption = caption) +
            theme_minimal() +
            theme(axis.text.x=element_text(angle = angle, hjust = hjust, vjust = vjust)) + 
            geom_hline(yintercept = yintercept, size=0.5, color='red', alpha=0.8, linetype = 'dotted')
    
    if (!is.null(file)) {
        savepdf(file)
            print(p)
         dev.off()
    } else {
       print(p)
    }
    if (return_data) { return(s) }
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


# texreg function
extractBRMS = function(model, r2 = FALSE) {
    ff = summary(model)$fixed
    rf = try(summary(model)$random)
    coefnames = rownames(ff)
    coefs = ff[, 1]
    se = ff[, 2]
    ci.low = ff[, 3]
    ci.up = ff[, 4]

    if (!is.null(rf)) {
        for (i in names(rf)) {
            coefnames = c(coefnames, rownames(rf[[i]]))
            coefs = c(coefs, rf[[i]][, 1])
            se = c(se, rf[[i]][, 2])
            ci.low = c(ci.low, rf[[i]][, 3])
            ci.up = c(ci.up, rf[[i]][, 4])
        }
    }

    gof = numeric()
    gof.names = character()
    gof.decimal = logical()

    n = stats::nobs(model)
    gof = c(gof, n)
    gof.names = c(gof.names, "Num. obs.")
    gof.decimal = c(gof.decimal, FALSE)


    if (r2) {
        rs = brms::bayes_R2(model)[1]
        gof = c(gof, rs)
        gof.names = c(gof.names, "R$^2$")
        gof.decimal = c(gof.decimal, TRUE)
    }

    tr = texreg::createTexreg(
        coef.names = coefnames,
        coef = coefs,
        se = se,
        ci.low = ci.low,
        ci.up = ci.up,
        gof.names = gof.names,
        gof = gof,
        gof.decimal = gof.decimal)
    return(tr)
}


# other functions

getMax = function(x) {
  x = na.omit(x)
  if (length(x) == 0) {
    return(NA_real_)
  } else {
    return(max(x))
  }
}


getMin = function(x) {
  x = na.omit(x)
  if (length(x) == 0) {
    return(NA_real_)
  } else {
    return(min(x))
  }
}


getFirstValue = function(x) {
  return(head(na.omit(x), 1))
}


getLastValue = function(x) {
  return(tail(na.omit(x), 1))
}


# interaction + scaling model
simScalingInteraction = function(E, a0 = 0.0,  a1 = 0.5, a2 = 0.10, b0 = 0.8, b1 = 0.2, h2 = 0.5) {
    N = length(E)
    G = rnorm(N,0,1)
    h = sqrt(h2)
    e = sqrt(1-h^2)
    sigma = exp(b0*e + b1*e*E) 
    y = rnorm(N, a0 + a1*E + a2*G*E + b0*h*G + b1*h*E*G, sigma)
    df[, qE := cut(E, quantile(E, probs = 0:10/10),
        labels = FALSE, include.lowest = TRUE)]
    df[, zy := scale(y)][, zg := scale(g)]
    df[, gzy := scale(y), qE]
    df[, gzg := scale(g), qE]
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
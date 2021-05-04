# auxiliary functions
# author: sebastian daza


plot_multi_histogram <- function(df, feature, label_column) {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text = label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", 
        linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
    plt + guides(fill=guide_legend(title=label_column))
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
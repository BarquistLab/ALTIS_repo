#' Fraction of samples that fall into the interval [-int0, int0], which corresponds to the null hypothesis H0.
#' The p value corresponds to the probability of the null hypothesis.
p_value <- function(x, int0=0, h0=0) {
  x_med <- median(x)
  ifelse(x_med>h0, length(which(x<int0+h0))/length(x), length(which(x>-int0+h0))/length(x))
}

p_value_alt <- function(x, logFC_h0) {
  logFC_med <- median(x)
  df <- rbind(data.frame(type="data", logFC=x),
              data.frame(type="h0", logFC=logFC_h0))
  d1dens <- with(df, density(logFC[type == "data"], 
                             from = min(logFC), 
                             to = max(logFC),
                             n = 1000))
  d2dens <- with(df, density(logFC[type == "h0"], 
                             from = min(logFC),
                             to = max(logFC),
                             n = 1000))
  if (logFC_med<0) {
    ind_max <- which(d1dens$y/d2dens$y<1) %>% min()
    return(1 - sum(d1dens$y[1:ind_max] - d2dens$y[1:ind_max]) / sum(d1dens$y))
  } else {
    ind_min <- which(d1dens$y/d2dens$y<1) %>% max()
    return(1 - sum(d1dens$y[ind_min:length(d1dens$y)] - d2dens$y[ind_min:length(d1dens$y)]) / sum(d1dens$y))
  }
}

logit <- function(x) {log(x/(1-x))}
logit_m1 <- function(x) {exp(x)/(1+exp(x))}

#' data frame containing the required diagnostic information for the
#' bayesplot package
extract_np <- function(mfit_diag, n_iter=1000) {
  
  n_chains <- nrow(mfit_diag)/n_iter
  mfit_diag$Iteration <- rep(1:n_iter, n_chains)
  mfit_diag$Chain <- rep(1:n_chains, each=n_iter)
  np <- pivot_longer(mfit_diag, -c("Iteration", "Chain"), names_to="Parameter", values_to="Value")
  colnames(np) <- c("Iteration", "Chain", "Parameter", "Value")
  np <- np[,c(1,3,4,2)]
  
  np <- subset(np, Parameter != "lp__")
  np$Parameter <- factor(np$Parameter, levels=c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__"))
  
  return(np)
}

plot_pairs <- function(plot_name, mfit, np, pars) {
  
  color_scheme_set("darkgray")
  p <- mcmc_pairs(mfit, np = np, pars=pars,
                  np_style = pairs_style_np(div_color="green", div_alpha=0.8))
  print(p)
  ggsave(plot_name, p, width = 10, height=10)
}

plot_trace <- function(plot_name, mfit, np, pars) {
  
  color_scheme_set("viridis")
  p <- mcmc_trace(mfit, np = np, pars=pars,
                  facet_args = list(ncol = 1, strip.position = "left"))
  print(p)
  ggsave(plot_name, p, width = 10, height=10)
}

# calculates how much the data has increased the probability of the alternative hypothesis
Bayes_factor <- function(logFC, mu=0, sigma=1, h0=0) {
  h1 <- median(logFC)
  post_dens <- density(logFC, from=min(logFC, h0), to=max(logFC, h0), n=1000) %>% approxfun()
  post_h1 <- post_dens(h1)
  post_h0 <- post_dens(h0)
  prior_h1 <- dnorm(h1, mu, sigma)
  prior_h0 <- dnorm(h0, mu, sigma)
  return(post_h0*prior_h1/post_h1/prior_h0)
}

Bayes_factor_mix <- function(logFC, mu=0, sigma1=1, sigma2=5, p, h0=0) {
  h1 <- median(logFC)
  post_dens <- density(logFC, from=min(logFC, h0), to=max(logFC, h0), n=1000) %>% approxfun()
  post_h1 <- post_dens(h1)
  post_h0 <- post_dens(h0)
  prior_h1 <- p*dnorm(h1, mu, sigma1) + (1-p)*dnorm(h1, mu, sigma2)
  prior_h0 <- p*dnorm(h0, mu, sigma1) + (1-p)*dnorm(h0, mu, sigma2)
  return(post_h0*prior_h1/post_h1/prior_h0)
}

calc_FDR_z <- function(z, p0, mean=0, sd=1) {
  N_g <- length(z)
  sapply(z, function(zi) {
    P <- which(z<=zi) %>% length()
    FP <- 2*p0*pnorm(zi, mean=mean, sd=sd)*N_g
    FP/P
  })
}

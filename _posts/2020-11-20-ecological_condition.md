---
title: "Effects of forest certification on the ecological condition of streams - a Bayesian approach"
date: 2020-11-20
tags: [Bayesian statistics]
header:
  image: "/images/eco_cond/banner_eco_cond2.jpg"
excerpt: "Bayesian statistics, Rivers, Ecological condition"
mathjax: "true"
---
**Motivation**

As a learning exercise, I will be revisiting some of my papers and repeating the statistical analysis with Bayesian methods. I'm still learning Bayesian statistics so do not assume that what I do is correct. If you do find any errors, please let me know.

I'm going to start with the following paper. You can get the data for this analysis [here](https://github.com/fsdias/blog_data/tree/main/ecol_cond).

Dias, Filipe S., Miguel N. Bugalho, Patricia M. Rodríguez-González, António Albuquerque, and J. Orestes Cerdeira. “Effects of Forest Certification on the Ecological Condition of Mediterranean Streams.” Edited by Angela Strecker. Journal of Applied Ecology 52, no. 1 (February 2015): 190–98. [https://doi.org/10.1111/1365-2664.12358](https://doi.org/10.1111/1365-2664.12358)

In this paper, we assessed the effects of Forest Stewardship Council (FSC) certification on the ecological condition of streams in Mediterranean evergreen oak woodlands. We used the Stream Visual Assessment Protocol (SVAP) to compare the ecological condition of stream reaches located in areas with 3 and 5 years of certification, non-certified regions, and in least-disturbed streams. SVAP scores range between 0 and 10. Streams with a higher SVAP score are in better ecological condition than those with lower SVAP scores. We expected stream reaches in certified estates to be in better ecological condition than those that are not.

**Exploratory data analysis**


<img src="{{ site.url }}{{ site.baseurl }}/images/eco_cond/study_area.png" alt="linearly separable data">

This figure shows the location of sample sites. We visited streams located in six FSC certified estates. Three estates had been certified for three years, and the three other estates for five years. We also visited rivers in two least-disturbed regions for reference.

We start by plotting SVAP scores against management regime, that is 1) no certification, 2) three years of certification, 3) five years of certification and 4) low disturbance.

<img src="{{ site.url }}{{ site.baseurl }}/images/eco_cond/svap_boxplot.png" alt="linearly separable data">


SVAP scores in stream reaches with five years of certification seem higher than those on non-certified areas, but there's some overlap. Stream reaches with three years of certification seem to have similar SVAP scores to those of non-certified reaches.


**Modeling**

Before we start the analysis, we need to consider two potential sources of dependence in the data. Look at Fig.1 and consider inset 3). The SVAP scores inside the estate are probably not independent from the scores outside the estate. Why? Because the samples were collected in the *same* river. We need to consider this dependence. But there's more, SVAP scores in reaches separated by 100 meters are more likely to be similar than those separated by 300 or 400 meters, regardless of FSC certification's effect. We also need to account for this.

To sum things up, we need to fit a model that 1) allows us to discern the effect of FSC certification on SVAP scores, 2) that considers the dependence caused by surveying reaches in the *same river*, and 3) that considers the dependence caused by geographical proximity (i.e., spatial autocorrelation).

We model SVAP scores with the Bayesian equivalent of a Mixed-effects ANOVA with spatial autocorrelation.

We assume SVAP scores follow a Normal distribution parameterized by the mean "mu" and the standard deviation "sigma". We express "mu" as a function of three additive varying effects "regime", "estate" and "k" as follows: 

```r
SVAP ~ Normal(mu,sigma)
mu = regime[regime_id]+estate[estate_id] + k[obs]
```

"regime" represents the effect of management regimes, that is non-certified, 3 years of certification, 5 years of certification or least-disturbed. "estate" represents the effects of the estates and least disturbed regions. "k" is a varying intercept that is estimated in light of geographic distance. The variable "k" accounts for the spatial autocorrelation and was estimated from a Gaussian process in which we assume the covariance between any two samples declines exponentially with the squared distance between them. Concerning the autocorrelation, I'm following what I read on [Statistical Rethinking](https://xcelab.net/rm/statistical-rethinking/) from Richard McElreath. I got the corresponding Stan code [here](https://vincentarelbundock.github.io/rethinking2/14.html).

In order to fit the model in Stan we we first need to calculate distance matrix for the location of SVAP measurements. The points shapefile is available here.


```r
library(raster)
library(sp)
points<-shapefile("points.shp")
dist<-spDists(points,points)
colnames(dist)<-data$id_code
rownames(dist)<-data$id_code
dist<-dist/1000
```

Then we create a list with all the data we need to fit the model.

```r
data_list<-list(
  N=nrow(data),
  N_estate=length(unique(data$estate_name)),
  N_regime=length(unique(data$regime_code)),
  svap=data$svap,
  regime_id=data$regime_id,
  estate_id=data$estate_id,
  obs=seq(1:nrow(data)),
  dist=dist
) 
```
Now, it's time run Stan. Here's the code in full.


```r
// cov_GPL2 macro from McElreath
functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}
data {
  int<lower=0> N;
  int<lower=0> N_regime;
  int<lower=0> N_estate;
  int<lower=1> obs[N];
  vector[N] svap;
  int<lower=1,upper=N_estate> estate_id[N];
  int<lower=1,upper=N_regime> regime_id[N];
  matrix[N, N] dist;
  }
parameters {
vector[N_regime] regime;
vector[N_estate] estate;
real mu_regime;
real mu_estate;
real<lower=0> sigma_dist;
real<lower=0> sigma_regime;
real<lower=0> sigma_estate;
vector[N] z;

real<lower=0> etasq;
real<lower=0> rhosq;
}

transformed parameters{
vector[N] mu;
//GP stuff
vector[N] k;
matrix[N,N] L_SIGMA;
matrix[N,N] SIGMA;
SIGMA = cov_GPL2(dist, etasq, rhosq, 0.01);
L_SIGMA = cholesky_decompose(SIGMA);
k = L_SIGMA * z;
mu= regime[regime_id]+estate[estate_id]+k[obs];
}

model {
  //GP stuff

  rhosq ~ exponential( 0.5 );
  etasq ~ exponential( 2 );
  z ~ normal( 0 , 1 );

  //priors
   regime~normal(mu_regime,sigma_regime);
   estate~normal(mu_estate,sigma_estate);
   mu_regime~normal(3,2);
   mu_estate~normal(3,2);
   
   sigma_dist~exponential(1);
   sigma_regime~exponential(1);
   sigma_estate~exponential(1);
   
   
   //Likelihood
   svap ~ normal(mu , sigma_dist);
   
}
generated quantities{
real svap_pred[N] = normal_rng(mu,sigma_dist);
real mu_r4_3 = regime[4] -regime[3];
real mu_r3_2 = regime[3] -regime[2];
real mu_r3_1 = regime[3] -regime[1];
real mu_r2_1 = regime[3] -regime[1];
real res[N];

for(i in 1:N){
  res[i] = svap_pred[i]-svap[i];
  
}

}
```

To run it, we use these instructions.


**Run Stan**

```r
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
model<-stan("model_gp.stan",data=data_list,control = list(adapt_delta = 0.99))
```

The run ends with no warnings, which means we can proceed to the validation phase. We extract predicted values and residuals, which we calculated during fitting (see generated quantities block).


**Validation**

```r
svap_pred<-as.matrix(model,pars=c("svap_pred"))
svap_pred_mean<-apply(svap_pred,2,mean)
svap<-data$svap
res<-as.matrix(model,pars=c("res"))
res_mean<-apply(res,2,mean)
estate_id<-data$estate_id
regime_id<-data$regime_id
```

Next, we load the packages that we are going to use during validation.
```r
library(bayesplot)
library(ggplot2)
library(gridExtra)
```

We start by plotting predicted SVAP scores against observed SVAP scores using a density plot. The fit looks reasonable.
```r
ppc_dens_overlay(svap,svap_pred[1:1000,])
```

<img src="{{ site.url }}{{ site.baseurl }}/images/eco_cond/density_plot.png" alt="linearly separable data">



But we should dig deeper by plotting residuals and predicted values against observed values. Remember that both residuals and predicted values have a distribution. Yeah, welcome to Bayesian stats. 

```r
grid.arrange(
ppc_intervals(res_mean,res,svap)+labs(y="Residuals",x="Observed SVAP ")+geom_abline(slope=0)+legend_none(),
ppc_intervals(svap_pred_mean,svap_pred,svap)+labs(y="Predicted SVAP",x="Observed SVAP ")+geom_abline(slope=1)+legend_none(),
ncol=2)
```

<img src="{{ site.url }}{{ site.baseurl }}/images/eco_cond/res_pred_vs_obs.png" alt="linearly separable data">



The result is not perfect, but I don't think it's terrible either. The model does seem to overestimate lower SVAP scores and underestimate higher ones. Still, if we look at the plot on the left, we see the residuals' distribution seems somewhat centered around zero.


Next, we should plot residuals and predicted values against the management regime.

```r
grid.arrange(
ppc_intervals_grouped(svap_pred_mean,svap_pred,svap,group=regime_id)+labs(y="Predicted SVAP",x="Observed SVAP ")+geom_abline(slope=1)+legend_none(),
ppc_intervals_grouped(res_mean,res,svap,group=regime_id)+labs(y="Residuals",x="Observed SVAP ")+geom_abline(slope=0)+legend_none(),
ncol=2)
```

<img src="{{ site.url }}{{ site.baseurl }}/images/eco_cond/res_pred_vs_regime.png" alt="linearly separable data">


We see the fit is worse for regime_id=2, which corresponds to three years of certification. We repeat the plot but with "estate_id."

```r
grid.arrange(
ppc_intervals_grouped(svap_pred_mean,svap_pred,svap,group=estate_id)+labs(y="Predicted SVAP",x="Observed SVAP ")+geom_abline(slope=1)+legend_none(),
ppc_intervals_grouped(res_mean,res,svap,group=estate_id)+labs(y="Residuals",x="Observed SVAP ")+geom_abline(slope=0)+legend_none(),
ncol=2)
```

<img src="{{ site.url }}{{ site.baseurl }}/images/eco_cond/res_pred_vs_estate.png" alt="linearly separable data">



The fit is not great for estates 1, 7, and 8. Note that 7 and 8 correspond to reference regions. We found some stream reaches in worse condition than we would expect in the least-disturbed areas.


**What's the effect of forest certification on the ecological condition of streams?**

To answer this question, we analyze the posterior difference between "mu_regime" values. The plot on the left shows the difference between the "mu_regime"  corresponding to "five years of certification" and the "mu_regime"  for three years of certification. We see the difference is reliably positive. The same happens when compare against the "mu_regime" for "no certification." Interestingly the model does not seem to be able to resolve the difference between "five years of certification"" and "reference sites." However, most of the probability mass is close to zero. This means SVAP scores between these two groups are likely similar.

```r
x<-as.matrix(model,pars=c("mu_r3_2","mu_r3_1","mu_r4_3"))
mcmc_hist(x)+ggtitle(c("5 years vs 3 years | 5 years vs no certification | 5 years vs Reference sites"))
```


<img src="{{ site.url }}{{ site.baseurl }}/images/eco_cond/posterior_diff.png" alt="linearly separable data">


**Conclusion**

With Bayesian methods, we got to the same conclusion I obtained when I used frequentist methods. However, the validation plots of the Bayesian model are not as clean. Maybe the technique we used to account for the spatial autocorrelation (Gaussian process) wasn't as successful in capturing the data's structure as its frequentist counterpart. Who knows. If you have any suggestions, I'm happy to hear them.

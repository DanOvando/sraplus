
# sraplus

[![Build
Status](https://travis-ci.org/DanOvando/sraplus.svg?branch=v2.0)](https://travis-ci.org/DanOvando/sraplus)
<!-- README.md is generated from README.Rmd. Please edit that file -->

**THIS PACKAGE IS UNDER ACTIVE DEVELOPMENT AND SHOULD NOT BE USED FOR
INDIVIDUAL STOCK ASSESSMENT AT THIS TIME**

`sraplus` is a flexible assessment package based around Ovando et
al. 2019. At the most “data limited” end, the model approximates the
behavior of catch-msy, sampling from prior distributions to obtain
parameter values that given a catch history do not crash the population
and satisfy supplied priors on initial and final depletion. At the most
data rich end the model can be fit to an abundance index or
catch-per-unit-effort data, while incorporating priors on recent stock
status based on Fisheries Management Index (FMI) scores or swept-area
ratio data.

## Installing

This is an in-development package hosted on github, so you will need to
do a few things to install it.

1.  open R

2.  If you don’t have the `devtools` package installed yet, run

<!-- end list -->

``` r
install.packages("devtools")
```

You’ll need to be connected to the internet.

3.  Once `devtools` has been installed, you can then install `sraplus`
    by
running

<!-- end list -->

``` r
remotes::install_github("danovando/sraplus", ref = 'v2.0', build_vignettes = TRUE)
```

That’s probably going to ask you to install many other packages, agree
to the prompts.

### Troubleshooting

Make sure you try the install with a fresh R session (go to
“Session\>Restart R” to make sure)

If you run into an error, first off try updating your R packages. From
there….

If your version of R is lower than 3.5, you might want to consider
updating R itself. Updating from 3.51 to 3.52 shouldn’t be any hassle.
BIG WARNING THOUGH, updating from say R 3.1 to 3.5 is a major update,
and you’ll lose all your installed packages in the process. I recommend
following the instructions
[here](https://www.datascienceriot.com/r/upgrade-R-packages/) to deal
with that, but even with that fix it can take a while, so I don’t
recommend doing a major R update if you’re on a deadline. There are also
packages to help you with this process, specifically
[`installR`](https://github.com/talgalili/installr/issues) for Windows
and [`updateR`](https://github.com/AndreaCirilloAC/updateR) for Mac.

From there…

  - On Windows, make sure you have the appropriate version of Rtools
    installed ([here](https://cran.r-project.org/bin/windows/Rtools/)),
    most likely Rtools35 if you have R version 3.3 or higher
      - Make sure that you select the box that says something about
        adding Rtools to the PATH variable
  - On macOS, there might be some issues with the your compiler. If you
    get an error that says something like `clang: error: unsupported
    option '-fopenmp'`, follow the instructions
    [here](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/)

Once you’ve tried those, restart your computer and try running

``` r

install.packages("devtools")

remotes::install_github("danovando/sraplus", ref = 'v2.0')
```

# Example Use

First, load libraries. These are the core libraries needed to run
`sraplus`, all of which should have been installed for you if you didn’t
already have them when you installed `sraplus`.

``` r
library(ggplot2)
library(tidyr)
library(dplyr)
library(sraplus)
library(tmbstan)
```

Once you’ve successfully installed `sraplus` you can take for a test
drive with these examples. We’ll now work through a couple examples,
from simple to complex, of using `sraplus`

## “Catch-only” SIR model

For the first example we’ll use a sampling-importance-resampling (SIR)
algorithm with a “catch-only” version of the model. In this mode, we use
catch heuristics (in the manner of the @anderson2014 adaptation of
Catch-MSY). The catch-heuristics are if catch in the first year is less
than 70% of the maximum reported catch, we set an expectation of
depletion in the initial year of 70% of carrying capacity. Otherwise,
the expectation is 40% of carrying capacity. For final depletion, the
heuristic prior is if catch in the final year is greater than 50% of
maximum catch final depletion is assumed to be 60% of carrying capacity,
otherwise 40% of carrying capacity.

The first step in running `sraplus` is the `sraplus::format_driors`
(data and priors) function. In this case, we will use example data for
Atlantic cod (*Gadus morhua*) that is included in the `sraplus` package.

``` r

example_taxa <- "gadus morhua"

data(cod)

head(cod)
#> # A tibble: 6 x 50
#>   scientificname commonname  year  catch stocklong.x TBmsybest ERmsybest
#>   <chr>          <chr>      <dbl>  <dbl> <chr>           <dbl>     <dbl>
#> 1 Gadus morhua   Atlantic …  1963 118000 Atlantic c…   948996.     0.262
#> 2 Gadus morhua   Atlantic …  1964 145000 Atlantic c…   948996.     0.262
#> 3 Gadus morhua   Atlantic …  1965 199000 Atlantic c…   948996.     0.262
#> 4 Gadus morhua   Atlantic …  1966 241000 Atlantic c…   948996.     0.262
#> 5 Gadus morhua   Atlantic …  1967 288000 Atlantic c…   948996.     0.262
#> 6 Gadus morhua   Atlantic …  1968 294000 Atlantic c…   948996.     0.262
#> # … with 43 more variables: TBmgtbest <dbl>, ERmgtbest <lgl>, TBmsy <dbl>,
#> #   SSBmsy <lgl>, Nmsy <lgl>, MSY <dbl>, Fmsy <dbl>, ERmsy <dbl>,
#> #   TBmgt <lgl>, SSBmgt <dbl>, Fmgt <dbl>, ERmgt <lgl>, TB0 <lgl>,
#> #   SSB0 <lgl>, M <lgl>, TBlim <dbl>, SSBlim <dbl>, Flim <dbl>,
#> #   ERlim <lgl>, b_v_bmsy <dbl>, u_v_umsy <dbl>, effort <lgl>,
#> #   total_biomass <dbl>, ss_biomass <dbl>, tsn <dbl>, areaid <chr>,
#> #   stocklong.y <chr>, region <chr>, inmyersdb <dbl>, myersstockid <lgl>,
#> #   tb_v_tb0 <lgl>, ssb_v_ssb0 <lgl>, delta_year <dbl>,
#> #   missing_gaps <lgl>, n_years <dbl>, has_tb0 <lgl>, has_tb <lgl>,
#> #   first_catch_year <dbl>, pchange_effort <dbl>, cs_effort <dbl>,
#> #   index <dbl>, approx_cpue <dbl>, b_rel <dbl>
```

From there, we’ll pass the catch data, and the years corresponding to
the catch data, to `format_driors`, and tell the model to use catch-msy
style heuristics by setting `use_hueristics = TRUE`. You’ll note that
`sraplus` prints a warning `use_hueristics = TRUE` reminding you that
any estimates of stock status are simply a transformation of your prior
beliefs expressed through the catch heuristics.

``` r


catch_only_driors <- sraplus::format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  use_heuristics = TRUE
)
#> Warning in sraplus::format_driors(taxa = example_taxa, catch = cod$catch, :
#> WARNING: You are using catch heursitics as your stock assessment
```

You can take a look at the information in the `catch_only_driors` object
by using `sraplus::plot_driors`

``` r

sraplus::plot_driors(catch_only_driors)
```

![](README_files/figure-gfm/c-msy-3-1.svg)<!-- -->

From there, we pass the `driors` object to `sraplus::fit_sraplus`, and
plot the results using `sraplus::plot_sraplus`. The `engine` argument
specifies how the model will be fit. When not actually “fitting” to
anything (rather simply sampling from priors that don’t crash the
population), we recommend setting engine to “sir”. The `draws` argument
tells `sraplus` how many draws from the SIR algorithm to generate, and
`n_keep` how many draws to sample from the total `draws` iterations. So
in this case the SIR algorithm will run 1 million iterations, and sample
2000 entries from those million in proportion to their likelihood.

``` r

 catch_only_fit <- sraplus::fit_sraplus(driors = catch_only_driors,
                       engine = "sir",
                       draws = 1e6,
                       n_keep = 2000)
```

Running `fit_sraplus` always produces a list with two objects: `results`
and `fit`. `results` is (mostly) standardized across engines set of
summary outputs from the fitted model. This allows us to easily plot and
compare outputs from models fit using `sir`, TMB, or stan. The `fit`
object contains the actual fitted model, which will of course vary
dramatically depending on what engine was used.

Let’s take a quick look at the `results` object.

``` r

head(catch_only_fit$results)
#> # A tibble: 6 x 6
#>    year variable           mean           sd       lower       upper
#>   <dbl> <chr>             <dbl>        <dbl>       <dbl>       <dbl>
#> 1  1963 b_div_bmsy        0.576       0.123        0.424       0.735
#> 2  1963 b           3395972.    1143302.     2306777.    4850974.   
#> 3  1963 c_div_msy         0.308       0.0947       0.186       0.430
#> 4  1963 depletion         0.389       0.0802       0.290       0.492
#> 5  1963 index_hat_t  163923.     115433.       32414.     311194.   
#> 6  1963 u_div_umsy        0.538       0.130        0.363       0.689
```

`results` is organized as a dataframe tracking different variables over
years. `mean` is the mean value for a given variable, `sd` the estimated
standard deviation around the mean, and `lower` and `upper` being the
10th and 90th percentiles of the estimates.

You can access other variables from the raw `fit` object, though this is
not standardized by engine and so requires knowing how to for example
get outputs out of `stanfit` objects. In the `catch_only_fit`, the `fit`
object is the output of the SIR algorithm.

``` r
head(catch_only_fit$fit)
#>   variable year draw   value
#> 1      b_t 1963    1 3572752
#> 2      b_t 1964    1 3957488
#> 3      b_t 1965    1 4361816
#> 4      b_t 1966    1 3995146
#> 5      b_t 1967    1 3948414
#> 6      b_t 1968    1 4208618
```

From there, we can generate some standard plots of B/Bmsy
(b\_div\_bmsy), Catch/MSY, U/Umsy, and depletion over time using
`plot_sraplus`.

``` r

sraplus::plot_sraplus(catch_only = catch_only_fit, years = catch_only_driors$years)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.svg)<!-- -->

## Fisheries Management Index and Swept Area Ratio

Now suppose that we obtain some FMI and SAR data for this fishery. We
can use these values to provide updated priors on current fishing
mortality rates and stock status (see full report for details on how
this is accomplished). Note that the FMI and SAR values year are
entirely fictional and any resemblance to any real fishery is purely
coincidental\!

You’ll notice that we now add a few more options to format\_driors.
We’ll manually set priors on initial depletion, with a prior of
initial biomass equal to carrying capacity (`initial_b = 1`), with a
standard deviation of 0.2. We’ll explicitly tell the model not to use
catch heuristics (though you don’t always need to specify this, FALSE is
the default). We’ll then pass the `driors` a swept area ratio of 2 (`sar
= 2`), and a named vector of FMI scores (where FMI scores range from 0
to 1). Note that FMI scores should be obtained through the formal FMI
survey process and not made up on the spot. W

``` r

fmi_sar_driors <- format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  initial_b = 1,
  initial_b_sd = 0.2,
  use_heuristics = FALSE,
  sar = 3,
  fmi = c("research" = 0.85,"management" = 0.75, "enforcement" = 0.75, "socioeconomics" = 0.75),
)

plot_driors(fmi_sar_driors)
```

![](README_files/figure-gfm/fmi-sar-1-1.svg)<!-- -->

We’ll then fit and plot our model. Note that you can pass multiple
`sraplus` fits, and can name each run whatever you would like.

``` r
fmi_sar_fit <- fit_sraplus(
  driors = fmi_sar_driors,
  engine = "sir",
  draws = 1e6,
  n_keep = 2000
)

plot_sraplus(fmi_sar = fmi_sar_fit,
             catch_only = catch_only_fit,
             years = fmi_sar_driors$years)
```

![](README_files/figure-gfm/fmi-sar-2-1.svg)<!-- -->

## Abundance Index via Maximum Likelihood

We’ll now try adding in some actual data to fit to. For illustration
purposes (and since we’re no longer using FMI/SAR data which can’t
really be simulated), we’ll use a simulated fishery. Let’s start with a
very simple example,using a simple fishery simulator built into
`sraplus`.

``` r

set.seed(42)
sim <-
  sraplus_simulator(
    sigma_proc = 0,
    sigma_u = 0,
    q_slope = 0,
    r = 0.4,
    years = 25,
    q = 1e-3,
    m = 1.01,
    init_u_umsy = 1
  )

sim$pop %>% 
  select(year, depletion,catch, effort,u) %>% 
  gather(metric, value, -year) %>% 
  ggplot(aes(year, value)) + 
  geom_point() + 
  facet_wrap(~metric, scales = "free_y") + 
  labs(y = "Value", x = "Year") + 
  sraplus::theme_sraplus()
```

![](README_files/figure-gfm/sim-index-1-1.svg)<!-- -->

Now, let’s pretend that we have a perfect index of abundance, which is
just `biomass * 1e-3`. We pass indices to `srplus` inside
`format_driors` as `index` and `index_years`, where `index_years` is a
vector of the same length of `index` specifying which years index data
are available. Well now use Template Model Builder (TMB) to estimate
stock status based on this index of abundance. We’ll add in some priors
on the growth rate and the shape of the Pella-Tomlinson model (1.01
roughly corresponds to a Fox model, where Bmsy/K ~= 0.4). Note that we
now set `engine = "tmb"` to fit the model via maximum likelihood using
TMB.

``` r
index_driors <- format_driors(
  catch = sim$pop$catch,
  years = sim$pop$year,
  index = sim$pop$biomass * 1e-3,
  index_years = sim$pop$year,
  growth_rate = 0.4,
  growth_rate_cv = 0.1,
  shape_prior = 1.01,
  shape_cv = 0.1)

plot_driors(index_driors)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.svg)<!-- -->

``` r

index_fit <- fit_sraplus(driors = index_driors,
                      engine = "tmb",
                      model = "sraplus_tmb")

plot_sraplus(index = index_fit,years = index_driors$years)
```

![](README_files/figure-gfm/unnamed-chunk-10-2.svg)<!-- -->

Looks good, now let’s try something a bit trickier.

## Fit Bayesian CPUE model with `stan`

We’ll now simulate a fishery with random-walk effort dynamics,
increasing catchability, and process error.

``` r

set.seed(42)

sim <-
  sraplus_simulator(
    sigma_proc = 0.05,
    sigma_u = 0.05,
    q_slope = 0.05,
    r = 0.2,
    years = 25,
    q = 1e-3,
    m = 1.01,
    init_u_umsy = 0.75
  )

sim$pop %>% 
  select(year, depletion,catch, effort,u) %>% 
  gather(metric, value, -year) %>% 
  ggplot(aes(year, value)) + 
  geom_point() + 
  facet_wrap(~metric, scales = "free_y") + 
  labs(y = "Value", x = "Year") + 
  sraplus::theme_sraplus()
```

![](README_files/figure-gfm/cpue-fit-1-1.svg)<!-- -->

Now suppose we no longer have a perfect index of abundance, but instead
data on the catch and effort (CPUE\!). But, there are a few problems
with these CPUE data. First, we know from the simulation that q is
increasing over time, so simply treating Catch/Effort as an index of
abundance will be biased (since a unit of effort in the past is not the
same as a unit of effort in the present). Second, we need to account for
diminishing returns from increasing amounts of effort, and in-season
losses to natural mortality. `sraplus` provides some support to these
problems. `sraplus` will estimate a q. If desired it allows the user to
either specify an assumed slope in catchability (where
\(q_t = q_{t-1} \times qslope\)), or to attempt to estimate `qslope`
directly. Second, given values of q and qslope (and estimates of natural
mortality either supplied by the user or drawn from `FishLife`),
`sraplus` uses the Baranov equation to translate effort into an
effective fishing mortality rate.

**One important note**. By default, `sraplus` includes estimation of
process error. When running a simplified CPUE like this, the model can’t
really handle estimating both process error and a q\_slope (since the
persistent trend in the CPUE values caused by the qslope can be soaked
into the process error or the qslope). So, you need to provide a VERY
imformative prior on the q\_slope parameter if you’re going to try and
estimate (i.e. fix the q\_slope parameter), or turn off process error
(inside `fit_sraplus` set `estimate_proc_error = FALSE`) (the
recommended option in this case).

By now the order of operations should be pretty familiar: pass things to
driors, then driors to fit\_sraplus. In this case, instead of passing an
index, we pass effort data, and effort years.

Just to explore functionality of `sraplus`, we’ll fit the model using
Bayesian estimation through stan (`engine = "stan"`). We’ll compare two
versions, one trying to estimate qslope, and one not. Note that we can
pass standard `rstan` options to `fit_sraplus`.

``` r

library(tmbstan)

cpue_driors <- format_driors(taxa = example_taxa,
                           catch = sim$pop$catch,
                           years = sim$pop$year,
                           effort = sim$pop$effort,
                           effort_years = sim$pop$year,
                           growth_rate = 0.4,
                           growth_rate_cv = 0.1,
                           shape_prior = 1.01,
                           q_slope = 0.025,
                           q_slope_cv = 0.25)


cpue_fit <- fit_sraplus(driors = cpue_driors,
                             engine = "stan",
                             model = "sraplus_tmb",
                             adapt_delta = 0.9,
                             max_treedepth = 10,
                             n_keep = 4000,
                             chains = 1, 
                             cores = 1,
                             estimate_qslope = FALSE)
#> Warning: There were 3 divergent transitions after warmup. Increasing adapt_delta above 0.9 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: Examine the pairs() plot to diagnose sampling problems

cpue_qslope_fit <- fit_sraplus(driors = cpue_driors,
                             engine = "stan",
                             model = "sraplus_tmb",
                             adapt_delta = 0.9,
                             max_treedepth = 10,
                             n_keep = 4000,
                             chains = 1, 
                             cores = 1,
                             estimate_qslope = TRUE,
                             estimate_proc_error = FALSE)
#> Warning: There were 205 divergent transitions after warmup. Increasing adapt_delta above 0.9 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: There were 849 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess
```

``` r

plot_sraplus(`CPUE fit no qslope` = cpue_fit, `CPUE fit with qslope` =  cpue_qslope_fit, years = cpue_driors$years)
```

![](README_files/figure-gfm/unnamed-chunk-12-1.svg)<!-- -->

As a final step, we can try adding in some fictional SAR data to our
fake fishery, just to see how it works. We can weight the SAR data using
the `sar_sd` input. Leaving `sar_sd = NA` uses the srandard deviation
from the posterior predictive of the fitted relationahip between SAR and
U/Umsy contained in the model. In other words, setting sar\_sd = NA lets
the data tell you how much weight to assign to the SAR data. You can
however overwrite this if desired and pass a stronger weight to the SAR
data if desired.

``` r
cpue_sar_qslope_driors <- format_driors(taxa = example_taxa,
                           catch = sim$pop$catch,
                           years = sim$pop$year,
                           effort = sim$pop$effort,
                           effort_years = sim$pop$year,
                           growth_rate = 0.4,
                           growth_rate_cv = 0.1,
                           shape_prior = 1.01,
                           q_slope = 0.025,
                           q_slope_cv = 0.25,
                           sar = 2,
                           sar_sd = NA)

cpue_sar_qslope_fit <- fit_sraplus(driors = cpue_sar_qslope_driors,
                             engine = "stan",
                             model = "sraplus_tmb",
                             adapt_delta = 0.9,
                             max_treedepth = 10,
                             n_keep = 4000,
                             chains = 1, 
                             cores = 1,
                             estimate_qslope = TRUE,
                             estimate_proc_error = FALSE)
#> Warning: There were 416 divergent transitions after warmup. Increasing adapt_delta above 0.9 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: There were 57 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess
```

And for good measure one more with SAR data and process error instead of
qslope

``` r

cpue_sar_proc_driors <- format_driors(taxa = example_taxa,
                           catch = sim$pop$catch,
                           years = sim$pop$year,
                           effort = sim$pop$effort,
                           effort_years = sim$pop$year,
                           growth_rate = 0.4,
                           growth_rate_cv = 0.1,
                           shape_prior = 1.01,
                           q_slope = 0,
                           q_slope_cv = 0.25,
                           sar = 4,
                           sar_sd = .05)

cpue_sar_proc_fit <- fit_sraplus(driors = cpue_sar_proc_driors,
                             engine = "stan",
                             model = "sraplus_tmb",
                             adapt_delta = 0.9,
                             max_treedepth = 10,
                             n_keep = 4000,
                             chains = 1, 
                             cores = 1,
                             estimate_qslope = FALSE,
                             estimate_proc_error = TRUE)
#> Warning: There were 353 divergent transitions after warmup. Increasing adapt_delta above 0.9 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: There were 1489 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 2.12, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess
```

``` r

plot_sraplus(`no rocess error and no qslope ` = cpue_fit, 
             `no process error with qslope` =  cpue_qslope_fit, 
             `no process error with qslope and sar` = cpue_sar_qslope_fit,
             `process error and sar` = cpue_sar_proc_fit,
             years = cpue_driors$years)
```

![](README_files/figure-gfm/unnamed-chunk-15-1.svg)<!-- -->

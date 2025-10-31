# Regression-Based Estimation of Causal Effects in the Presence of Selection Bias and Confounding

We consider the problem of estimating the expected causal effect $E[Y|do(X)]$ for a target variable $Y$ when treatment $X$ is set by intervention, focusing on continuous random variables. In settings without selection bias or confounding, $E[Y|do(X)] = E[Y|X]$, which can be estimated using standard regression methods. However, regression fails when systematic missingness induced by selection bias, or confounding distorts the data. Proxy variables unaffected by the selection process can, under certain constraints, be used to correct for selection bias to recover $E[Y|X]$, and hence $E[Y|do(X)]$, reliably. When data is additionally affected by confounding, however, recovering the causal effect from selection-biased data is more challenging and requires access to proxies to both correct for confounding and for the selection mechanism. Assuming access to such proxies from external unbiased observational data, we derive theoretical conditions ensuring identifiability and recoverability of causal effects. We further introduce a linear two-step regression estimator (TSR), which can be extended through adding non-linear basis functions, capable of exploiting proxy variables to adjust for selection bias while accounting for confounding. We show that TSR is consistent with previous estimators when confounding is absent, but achieves a lower variance. Extensive simulation studies validate TSR’s correctness for scenarios that include both selection bias and confounding with proxy variables.

## Project Structure

In the following, we provide the general project structure.

```bash
Template-R-Project/
├──sec4.1.R                  #simulations from section 4.1
├──sec4.2_ex12.R             #simulations from section 4.2  :  examples 1 and 2
├──sec4.2_ex34.R             #              "                  examples 3 and 4
├──sec4.2_ex56.R             #              "                  examples 5 and 6
├──sec4.3.R                  #simulations from section 4.3
├──sec4.4.R                  #simulations from section 4.4
```


## Installation Instructions

R and the R packages ggplot2 and glmnet (ridge regression) need to be installed.

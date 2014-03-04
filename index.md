---
title       : Quantile Regression in the Presence of Monotone Missingness with Sensitivity Analysis
subtitle    : 2014 ENAR Presentation
author      : Minzhao Liu
job         : 'Supervisor: Dr. Mike Daniels. Department of Statistics, University of Florida'
license     : by-nc-sa
framework   : io2012 #{io2012, html5slides, deckjs, shower, dzslides, ...}
deckjs:
   theme: web-2.0
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : github      # zenburn, github, tomorrow
widgets     : [mathjax, bootstrap]            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft, selfcontained}

---

## Outline

- Introduction
- Models
- Real Data Analysis: TOURS
- Summary

---

## Missing Data Mechanism

- Missing data mechanism: $$p(\mathbf r| \mathbf y, \mathbf x, \mathbf \phi(\mathbf \omega))$$
- Missing Complete At Random (MCAR)
$$
    p(\mathbf r |  y_{obs}, y_{mis}, \mathbf x, \mathbf \phi) =   p(\mathbf r | \mathbf x, \mathbf \phi).
$$
- Missing At Random (MAR)
$$
    p(\mathbf r |  y_{obs}, y_{mis}, \mathbf x, \mathbf \phi) =   p(\mathbf r | y_{obs}, \mathbf x, \mathbf \phi).
$$
- Missing Not At Random (MNAR), for $y_{mis} \neq y_{mis}\prime$,
$$p(\mathbf r |  y_{obs}, y_{mis}, \mathbf x, \mathbf \phi) \neq   p(\mathbf r | y_{obs}, y_{mis}\prime, \mathbf x, \mathbf \phi).$$

---

## Notation

- Under monotone dropout, WOLOG, denote $S_i \in \{1, 2, \ldots, J\}$ to be the number of observed $Y_{ij}'s$ for subject $i$,
- $\mathbf Y_i = (Y_{i1}, Y_{i2}, \ldots, Y_{iJ})^{T}$ to be the full data response vector for subject $i$,
- $J$ is the maximum follow up time.
- We assume $Y_{i1}$ is always observed.
- We are interested in the $\tau$-th marginal quantile regression coefficients $\mathbf \gamma_j = (\gamma_{j1}, \gamma_{j2}, \ldots, \gamma_{jp})^T$,
$$
  Pr (Y_{ij} \leq \mathbf x_i^{T} \mathbf \gamma_j ) = \tau, \mbox{ for } j = 1, \ldots, J,
$$
where $\mathbf x_i$ is a $p \times 1$ vector of covariates for subject $i$.
- $$
  p_k(Y) = p(Y | S = k), \quad  p_{\geq k} (Y)  = p(Y | S \geq k)
$$
be the densities of response $\mathbf Y$ given follow-up time $S=k$ and $S \geq k$. And $Pr_k$ be the corresponding probability given $S = k$.

---

## Pattern Mixture Model

- Mixture models factor the joint distribution of response and missingness as
$$
  p (\mathbf y, \mathbf S |\mathbf x, \mathbf \omega) = p (\mathbf y|\mathbf S, \mathbf x, \mathbf \omega) p (\mathbf S | \mathbf x, \mathbf \omega).
$$
- The full-data response distribution is given by
$$
  p (\mathbf y | \mathbf x, \mathbf \omega) = \sum_{S \in \mathcal{S}} p(\mathbf y| \mathbf S, \mathbf x, \mathbf \theta) p (\mathbf S | \mathbf x, \mathbf \phi),
$$
where $\mathcal{S}$ is the sample space for dropout time $S$ and the parameter vector $\mathbf \omega$ is partitioned as $(\mathbf \theta, \mathbf \phi)$.
- Furthermore, the conditional distribution of response within patterns can be decomposed as
$$
  P (Y_{obs}, Y_{mis} | \mathbf S, \mathbf \theta) = P
  (Y_{mis}|Y_{obs}, \mathbf S, \mathbf \theta_E) Pr (Y_{obs} | \mathbf S, \mathbf
  \theta_{y, O}),
$$

- $\mathbf \theta_E$:  extrapolation distribution
- $\mathbf \theta_{y, O}$ : distribution of observed responses

<!-- ### Sensitivity Parameters -->

<!-- A reparameterization $\mathbf \xi(\mathbf \alpha) = (\mathbf \xi_S, \mathbf \xi_M)$: -->

<!-- - $\mathbf \xi_S$ is a nonconstant function of $\mathbf \theta_E$, -->
<!-- -  the observed-data likelihood -->
<!-- $$ -->
<!--       L(\mathbf \xi_S, \mathbf \xi_M | y_{obs}, \mathbf r) -->
<!-- $$ -->
<!-- is constant as a function of $\mathbf \xi_S$, and -->
<!-- -  at a fixed value of $\mathbf \xi_S$, the observed data likelihood is a nonconstant function of $\mathbf \xi_M$, -->

<!-- then $\mathbf \xi_S$ is a sensitivity parameter. -->

---

## Model Settings

- Multivariate normal distributions within pattern
- The marginal quantile regression models as:
$$
  Pr (Y_{ij} \leq \mathbf x_{i}^T \mathbf \gamma_j ) = \tau,
$$
$$
  \begin{array}{l}
      \displaystyle p_k(y_{i1}) = N (\Delta_{i1} +  \mathbf \beta_1^{(k)},
      \sigma_1^{(k)}  ), k = 1, \ldots, J,\\
       \displaystyle p_k(y_{ij}|\mathbf y_{ij^{-}}) =
      \begin{cases}
        \textrm{N} \big (\Delta_{ij} + \mathbf y_{ij^{-}}^T \mathbf
        \beta_{y,j-1},
        \sigma_j \big), & k \geq j ;  \\
        \textrm{N} \big ( \chi(\mathbf x_{i}, \mathbf y_{ij^{-}}),
        \sigma_j \big), & k < j ;  \\
      \end{cases}, \mbox{ for } 2 \leq j \leq J,  \\
      \displaystyle S_{ij} = k|\: \mathbf x_{i} \sim \textrm{Multinomial}(1, \mathbf \phi),
    \end{array}
$$
- $\chi(\mathbf x_{i}, \mathbf y_{ij^{-}})$ is the mean of the unobserved data distribution and allows sensitivity analysis by varying assumptions on $\chi$; for computational reasons we assume that $\chi$ is linear in $y_{ij^{-}}$.
- For example, here we specify
$$
\chi(\mathbf x_{i}, \mathbf y_{ij^{-}}) = \Delta_{ij}  + \mathbf y_{ij^{-}}^T \mathbf \beta_{y,j-1} + h_{0}^{(k)}.
$$

---

## $\Delta$

$\Delta_{ij}$ are subject/time specific intercepts determined by the parameters in the model and
and are determined by the marginal quantile regressions,
$$
  \tau = Pr (Y_{ij} \leq \mathbf x_{i}^T \mathbf \gamma_j ) = \sum_{k=1}^J
  \phi_kPr_k (Y_{ij} \leq \mathbf x_{i}^T \mathbf \gamma_j ) \mbox{  for  } j = 1,
$$
and
$$
\begin{align}
  \tau &= Pr (Y_{ij} \leq \mathbf x_{i}^{T} \mathbf \gamma_j ) =
  \sum_{k=1}^J
  \phi_kPr_k (Y_{ij} \leq \mathbf x_{i}^{T} \mathbf \gamma_j ) \\
  & = \sum_{k=1}^J \phi_k \int\cdots \int Pr_k (Y_{ij} \leq \mathbf
  x_{i}^{T} \mathbf \gamma_j | \mathbf y_{ij^{-}}
  ) p_k (y_{i(j-1)}| \mathbf y_{i(j-1)^{-}})  \nonumber \\
  & \quad \cdots p_k (y_{i2}| y_{i1}) p_k(y_{i1})
  dy_{i(j-1)}\cdots dy_{i1}.  \mbox{  for  } j = 2, \ldots, J .\nonumber
\end{align}
$$

---

## Intuition

- Embed the marginal quantile regressions directly in the model through constraints in the likelihood of pattern mixture models
- The mixture model allows the marginal quantile regression coefficients to differ by quantiles. Otherwise, the quantile lines would be parallel to each other.
- The mixture model also allows sensitivity analysis.
- For identifiability of the observed data distribution, we apply the following constraints,
$$
 \sum_{k=1}^J \beta_{1}^{(k)} = 0.
$$

---

## Missing Data Mechanism and Sensitivity Analysis

- Mixture models are not identified due to insufficient information provided by observed data.
- Specific forms of missingness are needed to induce constraints to identify the distributions for incomplete patterns, in particular, the extrapolation distribution

- In mixture models , MAR holds (Molenberghs et al. 1998; Wang & Daniels, 2011) if and only if, for each $j \geq 2$ and $k < j$:
$$
  p_k(y_j|y_1, \ldots, y_{j-1}) = p_{\geq j}(y_j|y_1, \ldots, y_{j-1}).
$$
- When $2 \leq j \leq J$ and $k < j$, $Y_j$ is not observed, thus $h_0^{(k)}$ can not be identified from the observed data.

---

## Sensitivity Analysis

- $\mathbf \xi_s = (h_0^{(k)})$ is a set of sensitivity parameters
(Daniels & Hogan 2008), where $k =1, ..., J-1$.
- $\mathbf \xi_s = \mathbf \xi_{s0} = \mathbf 0$, MAR holds.
- $\mathbf \xi_s$ is fixed at $\mathbf \xi_s \neq \mathbf \xi_{s0}$, MNAR.
- We can vary $\mathbf \xi_s$ around $\mathbf 0$ to examine the impact of different MNAR mechanisms.

---

## Calculation of $\Delta_{ij}$ ($j = 1$)

$\Delta_{ij}$ depends on subject-specific covariates $\mathbf x_{i}$, thus $\Delta_{ij}$ needs to be calculated for each subject. We now illustrate how to calculate $\Delta_{ij}$ given all the other parameters $\mathbf \xi = (\mathbf \xi_m, \xi_s)$.

**$\Delta_{i1}: $** Expand equation :
$$
\begin{align*}
    \tau = \sum_{k = 1}^J \phi_k \Phi \left( \frac{\mathbf x_{i1}^T
        \mathbf \gamma_1 - \Delta_{i1} -
        \beta_1^{(k)}}{ \sigma_1^{(k)} } \right),
  \end{align*}
$$
  where $\Phi$ is the standard normal CDF. Because the above equation   is continuous and monotone in $\Delta_{i1}$, it can be solved by a   standard numerical root-finding method (e.g. bisection method) with   minimal difficulty.

---

## Calculation of $\Delta_{ij}, 2\leq j \leq J$

Lemma:
$$
\begin{array}{l}
\displaystyle \int \Phi \left( \frac{x-b}{a} \right) d\Phi(x; \mu, \sigma)  =
\begin{cases}
1- \Phi \left( \frac{b-\mu}{\sigma} \big /
\sqrt{\frac{a^2}{\sigma^2}+1} \right) & a > 0, \\
\Phi \left( \frac{b-\mu}{\sigma} \big /
\sqrt{\frac{a^2}{\sigma^2}+1} \right) & a < 0,
\end{cases}
\end{array}
$$

Recursively for the first multiple integral, apply lemma once to obtain:

$$
  \begin{align*}
    Pr_1 (Y_{ij} \leq \mathbf x_{i}^T \mathbf \gamma_j) & =
    \idotsint
    Pr_1 (Y_{ij} \leq \mathbf x_{i}^T\mathbf \gamma_j |\mathbf x_{i}, \mathbf Y_{ij^{-}})\\
    & \quad  dF_1(Y_{i(j-1)}|\mathbf x_{i}, \mathbf Y_{i(j-1)^{-}}) \cdots d F_1 (Y_{i1} |\mathbf x_{i}), \\
    & = \idotsint \Phi \left( \frac{Y_{i(j-2)} - b^{*}}{a^{*}}
    \right) dF_1(Y_{i(j-2)}|\mathbf x_{i}, \mathbf Y_{i(j-2)^{-}})\cdots d F_1 (Y_{i1} | \mathbf x_{i}).
  \end{align*}
$$

---

## MLE

The observed data likelihood for an individual $i$ with follow-up time $S_i = k$ is
$$
\begin{align}
L_i(\mathbf \xi| \mathbf y_i, S_{i} = k) & =
  \phi_kp_k (y_k | y_1, \ldots, y_{k-1})
  p_k (y_{k-1}|y_1, \ldots, y_{k-2}) \cdots p_{k} (y_1), \\
  & = \phi_k p_{\geq k} (y_k | y_1, \ldots, y_{k-1}) p_{\geq k-1}
  (y_{k-1}|y_1, \ldots, y_{k-2}) \cdots p_{k} (y_1), \nonumber
\end{align}
$$
- Use the bootstrap to construct confidence interval and make inferences.

## Goodness of Fit Check

- Check QQ plots of fitted residuals
$$
  \hat{\epsilon}_{ij} =
  \begin{cases}
    (y_{ij} - \hat{\Delta}_{ij} - \hat{\beta}_1^{(k)})/\hat{\sigma}_1^{(k)},& j = 1 \\
    (y_{ij} - \hat{\Delta}_{ij} - \mathbf{y_{ij^{-}}^T
    \hat{\beta}_{y,j-1}})/\hat{\sigma}_j, & j > 1
  \end{cases}.
$$

---
## Real Data Analysis: Tours

- Weights were recorded at baseline ($Y_0$), 6 months ($Y_1$) and 18 months ($Y_2$).
- We are interested in how the distributions of weights at six months and eighteen months change with covariates.
- The regressors of interest include **AGE**, **RACE** (black and white) and **weight at baseline** ($Y_0$).
- Weights at the six months ($Y_1$) were always observed and 13 out of 224 observations (6%) were missing at 18 months ($Y_2$).
- The **AGE** covariate was scaled to 0 to 5 with every increment representing 5 years.
- We fitted regression models for bivariate responses $\mathbf Y_i = (Y_{i1}, Y_{i2})$ for quantiles (10%, 30%, 50%, 70%, 90%).
- We ran 1000 bootstrap samples to obtain 95% confidence intervals.


---

## Results

- For weights of participants at six months, weights of whites are generally 4kg lower than those of blacks for all quantiles, and the coefficients of race are negative and significant.
- Weights of participants are not affected by age since the coefficients are not significant. Differences in quantiles are reflected by the intercept.
- Coefficients of baseline weight show a strong relationship with weights after 6 months.
- For weights at 18 months after baseline, we have similar results.
- Weights at 18 months still have a strong relationship with baseline weights.
- However, whites do not weigh significantly less than blacks at 18 months unlike at 6 months.

---

## Sensitivity Analysis

We also did a sensitivity analysis based on an assumption of MNAR.
- Based on previous studies of pattern of weight regain after lifestyle treatment (Wadden et al. 2001; Perri et al. 2008)
we assume that
$$
  E(Y_2 - Y_1| R=0) = 3.6 \mbox{kg},
$$
which corresponds to 0.3kg regain per month after finishing the initial 6-month program.
- Therefore, we specify  $\chi(\mathbf x_{i}, Y_{i1})$ as
$$
\chi(\mathbf x_{i},  y_{i1}) = 3.6  + y_{i1},
$$

### Results

- There are not large differences for estimates for $Y_2$ under MNAR vs MAR.
- This is partly due to the low proportion of missing data in this study.

---
## Summary

- Developed a marginal quantile regression model for data with monotone missingness.
- Used a pattern mixture model to jointly model the full data response and missingness.
- Estimate marginal quantile regression coefficients instead of conditional on random effects
- Allows for sensitivity analysis which is essential for the analysis of missing data (NAS 2010).
- Allows the missingness to be non-ignorable.
- Recursive integration simplifies computation and can be implemented in high dimensions.

## Future Work

- Sequential multivariate normal distribution for each component in the PMM might be too restrictive
- Simulation results showed that the mis-specification of the error term did have an impact on the extreme quantile regression inferences.
- Working on replacing it with a non-parametric model, for example, a Dirichlet process mixture of normals.

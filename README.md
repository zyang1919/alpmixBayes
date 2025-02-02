# alpmixBayes

**alpmixBayes** is an R package for Bayesian model fitting of alpha mixture survival models. It provides implementations for multiple mixture models, including Weibull-Weibull, Exponential-Weibull, Lognormal-Lognormal, and Exponential-Weibull-Gamma mixtures.

## 📌 Features
- Implements **four Bayesian mixture models**:
  - **Weibull-Weibull Mixture** (`wwmix`)
  - **Exponential-Weibull Mixture** (`ewmix`)
  - **Lognormal-Lognormal Mixture** (`llmix`)
  - **Exponential-Weibull-Gamma Mixture** (`ewgmix`)
- Supports **MCMC sampling** with user-defined priors.
- Provides **credible intervals** for parameter estimates.
- Easily compare different mixture models using `summary()`.

---

## 📥 Installation
Currently, `alpmixBayes` is not on CRAN. You can install it by cloning the repository:

```r
# Install dependencies (if needed)
install.packages(c("coda", "ggplot2"))

# Clone the repository
system("git clone https://github.com/zyang1919/alpmixBayes.git")

# Load functions manually
source("your path/alpmixBayes/alpmixBayes.R")

---

## 🚀 Usage
### **1️⃣ Load the Package**
```r
source("alpmixBayes.R")
```

## 2️⃣ Fit a Model
## Example dataset (replace with your actual data)
```r
data <- rnorm(100)

# Define MCMC settings
mcmc_values <- list(nburn = 1000, niter = 5000, thin = 1)

# Define initial values
init_values <- list(alpha = 1, p = 0.5, th1 = 2, th2 = 2)

# Fit a Weibull-Weibull mixture model
fit <- alpmixBayes(data, mcmc_values, init_values, survmodel = "WW")

# Print results
print(fit)
summary(fit)
```

## 3️⃣ Model Output
## The function returns parameter estimates and 95% credible intervals:
```r
Summary of alpmix model fitting:

Model used: Weibull-Weibull Mixture 

Parameter Estimates (with 95% Credible Intervals):
  Parameter  Estimate  Lower_95  Upper_95
1     alpha 0.9989897 0.9441765 1.0554252
2         p 0.4962427 0.4685414 0.5471943
3       th1 0.4936465 0.4754862 0.5226008
4       th2 1.0061340 0.9837623 1.0317516
```

## 📊 Model Options
You can specify which mixture model to use with the ```survmodel``` argument:
```r
fit_ww  <- alpmixBayes(data, mcmc_values, init_values, survmodel = "WW")   # Weibull-Weibull
fit_ew  <- alpmixBayes(data, mcmc_values, init_values, survmodel = "EW")   # Exponential-Weibull
fit_ll  <- alpmixBayes(data, mcmc_values, init_values, survmodel = "LL")   # Lognormal-Lognormal
fit_ewg <- alpmixBayes(data, mcmc_values, init_values, survmodel = "EWG")  # Exponential-Weibull-Gamma
```
## 📖 Documentation






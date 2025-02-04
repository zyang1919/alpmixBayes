# alpmixBayes

**alpmixBayes** is an R package for Bayesian alpha mixture survival models. It provides Bayesian estimation for alpha mixture models, including Weibull-Weibull mixture, Exponential-Weibull mixture, Lognormal-Lognormal mixture, and Exponential-Weibull-Gamma mixtures.

## 📌 Features
- Implements four Bayesian alpha mixture models under different distributions:
  - Weibull-Weibull Mixture (`wwmix`)
  - Exponential-Weibull Mixture (`ewmix`)
  - Lognormal-Lognormal Mixture (`llmix`)
  - Exponential-Weibull-Gamma Mixture (`ewgmix`)
- Supports MCMC sampling, priors, and initial values options.
- Provides credible intervals for parameter estimates.
- Easily compare different mixture models using `summary()`.
---
### **Arguments**
- **`d`**: A numeric vector of observed survival times. The input data should represent non-negative and non-censered survival times.
- **`mcmc_values`**: A list containing MCMC settings:
  - `nburn`: Number of burning iterations. Default value is 1,000. 
  - `niter`: Collections of MCMC iterations. Default value is 5,000. 
  - `thin`: Thinning interval (e.g., if `thin = 2`, every second iteration is kept). Default value is 1.
- **`init_values`**: A list specifying initial values for parameters. The required parameters depend on the chosen mixture model:
  - For **two-component models** (`"WW"`, `"EW"`, `"LL"`):  
    - `alpha`: Initial value for the mixing power parameter. Default value is 1. 
    - `p`: Initial mixing proportion. Default value is 0.5. 
    - `th1`, `th2`: Initial values for the two scale parameters. Default values are 1 and 1. 
  - For **three-component model** (`"EWG"`):  
    - `alpha`: Initial value for the mixing power parameter. Default value is 1.
    - `p1`, `p2`: Initial mixing proportions.Default value is 0.5. Default values are 0.3 and 0.3. 
    - `th1`, `th2`, `th3`: Initial values for the three scale parameters. Default values are 1, 1, and 1. 
- **`survmodel`**: A string indicating the mixture model to use. The available options are:
  - `"WW"`: Weibull-Weibull Mixture  
  - `"EW"`: Exponential-Weibull Mixture  
  - `"LL"`: Lognormal-Lognormal Mixture  
  - `"EWG"`: Exponential-Weibull-Gamma Mixture
- **`prior`**: A list specifying hyperparameters for Bayesian priors. If `prior = NULL`, default values are used:
  - For **two-component models** (`"WW"`, `"EW"`, `"LL"`):
    - `ap = c(0,0.001)`: Prior for the mixing power, assuming a **Normal(0, 0.001²)** distribution. 
    - `pp = c(1,1)`: Prior for the mixture proportion, assuming a **Beta(1,1)** distribution. 
    - `th1p = c(1,1)`: Prior for `th1`, assuming a **Gamma(1,1)** distribution.  
    - `th2p = c(1,1)`: Prior for `th2`, assuming a **Gamma(1,1)** distribution.  
  - For **three-component model** (`"EWG"`):
    - `ap = c(1,0.001)`: Prior for the mixing power, assuming a **Normal(1, 0.001²)** distribution.   
    - `pp = c(1,1,1)`: Prior for the mixture proportion, assuming a **Dirichlet(1,1)** distribution. 
    - `th1p = c(1,1)`: Prior for `th1`, assuming a **Gamma(1,1)** distribution. 
    - `th2p = c(1,1)`: Prior for `th2`, assuming a **Gamma(1,1)** distribution. 
    - `th3p = c(1,1)`: Prior for `th3`, assuming a **Gamma(1,1)** distribution.

### **Value**
The function returns a **list** containing the following elements:

- **`model`**: A string indicating the mixture model used. This will be one of `"WW"`, `"EW"`, `"LL"`, or `"EWG"`.  
- **`estimates`**: A data frame containing parameter estimates and their **95% credible intervals**. The structure varies depending on the model:
  - For **two-component models** (`"WW"`, `"EW"`, `"LL"`), the data frame includes:
    - `alpha`: Posterior mean estimate for the transformation parameter.
    - `p`: Posterior mean estimate for the mixing proportion.
    - `th1`, `th2`: Posterior mean estimates for the component parameters.
    - `Lower_95`, `Upper_95`: Lower and upper bounds of the 95% credible interval.
  - For **three-component models** (`"EWG"`), the data frame includes:
    - `alpha`: Posterior mean estimate for the transformation parameter.
    - `p1`, `p2`: Posterior mean estimates for the first and second mixture proportions.
    - `th1`, `th2`, `th3`: Posterior mean estimates for the component parameters.
    - `Lower_95`, `Upper_95`: Lower and upper bounds of the 95% credible interval.
- **`mcmc_samples`**: A matrix containing MCMC posterior draws for each parameter (optional, depending on user settings).  
- **`data`**: The original dataset used for model fitting.  
- **`call`**: The matched function call, displaying the arguments used in the function.  
---
## Example
## 📥 Installation
Currently, `alpmixBayes` is not on CRAN. You can install it by cloning the repository:

```r
# Install dependencies (if needed)
install.packages(c("coda", "ggplot2"))

# Clone the repository
system("git clone https://github.com/zyang1919/alpmixBayes.git")

# Load functions manually
source("your path/alpmixBayes/alpmixBayes.R")
```

## 🚀 Usage
### Load the Package
```r
source("alpmixBayes.R")
## Example dataset (replace with your actual data)
data <- rnorm(100)

# Define MCMC settings
mcmc_values <- list(nburn = 1000, niter = 5000, thin = 1)

# Define initial values
init_values <- list(alpha = 1, p = 0.5, th1 = 2, th2 = 2)

# Define prior values
# two-component mixture
prior <- list(ap=c(0,0.001), pp=c(1,1), th1p=c(1,1), th2p=c(1,1))
# three-component mixture
prior <- list(ap=c(1,0.001), pp=c(1,1,1), th1p=c(1,1), th2p=c(1,1), th3p=c(1,1))

# Fit a Weibull-Weibull mixture model
fit <- alpmixBayes(data, mcmc_values, init_values, prior, survmodel = "WW")

# Print results
print(fit)
summary(fit)
```

## Model Output
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
fit_ww  <- alpmixBayes(data, mcmc_values, init_values, prior, survmodel = "WW")   # Weibull-Weibull
fit_ew  <- alpmixBayes(data, mcmc_values, init_values, prior, survmodel = "EW")   # Exponential-Weibull
fit_ll  <- alpmixBayes(data, mcmc_values, init_values, prior, survmodel = "LL")   # Lognormal-Lognormal
fit_ewg <- alpmixBayes(data, mcmc_values, init_values, prior, survmodel = "EWG")  # Exponential-Weibull-Gamma
```

### **Authors**
This package was developed by:

- **Feng Luan** – Lead developer, coding  
- **Duchwan Ryu** – Contributor, coding 
- **Zhexuan Yang** – Contributor, coding, Github  






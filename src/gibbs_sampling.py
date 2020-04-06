import numpy as np
from scipy.stats import invgamma
import seaborn as sns
import matplotlib.pyplot as plt

# hyperparameters
mu_0 = 0
lambda_0 = 1
rho_0 = 0.1
delta_0 = 2.2

# iterations
iters = 10000

# burn-in phase
burn = 500



# data generating
n = 500
mu = .5
sigma = .5
y = np.random.normal(loc= mu, scale= sigma, size= n)


# gibbs sampling

def posterior_mu(y, sigma):
    
    mu_n = (lambda_0**(-1) + n * sigma**(-1))**(-1) * ( lambda_0**(-1) + n * sigma**(-1) * np.mean(y) )
    lambda_n = (lambda_0**(-1) + n * sigma**(-1))**(-1)
    
    return np.random.normal(loc=mu_n, scale=lambda_n)
    
    
def posterior_sigma(y, mu):
    
    ep = (y - mu)
    
    delta = ep.dot(ep)
    
    a = (rho_0 + n) / 2
    b = (delta_0 + delta) / 2
    
    return invgamma.rvs(a=a, scale=b)



mu_t = np.zeros(iters+1)
sigma_t = np.zeros(iters+1)
# initial guess
mu_t[0], sigma_t[0] = 0, 1

for i in range(iters):
    
    mu_t[i+1] = posterior_mu(y, sigma_t[i])
    sigma_t[i+1] = posterior_sigma(y, mu_t[i+1])


# discard burn-in phase


mu_t = mu_t[burn:]
sigma_t = np.sqrt(sigma_t[burn:])


# visualise

sns.jointplot(mu_t, sigma_t)
plt.show()


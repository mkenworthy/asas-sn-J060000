import snake
import numpy as np
import matplotlib.pyplot as plt

# testing an emcee fitting to the lines

# make a test data set

N = 101
t = np.linspace(58000, 58099.5, N)
np.random.seed(77)

t_inflexion = np.linspace(58000,58100, 12)
t_inflexion += np.random.uniform(0,5,size=t_inflexion.size)

f_min = 14

f_max = 15

f_inflexion = np.random.uniform(f_min, f_max, size=np.shape(t_inflexion))

fig1, ax1 = plt.subplots(1, 1, figsize=(12, 4))

ax1.set_xlim(np.min(t_inflexion),np.max(t_inflexion))

ax1.vlines(t_inflexion, f_min, f_max, linestyle='dashed')

ax1.scatter(t_inflexion, f_inflexion,color='orange')


(x_true,y_true) = snake.snake(t, t_inflexion, f_inflexion)

err = 0.1 * np.random.randn(y_true.size)

y_true_noisy = y_true + np.random.randn(y_true.size)*err

ax1.plot(x_true, y_true, 'k', alpha=0.3, lw=5)

ax1.errorbar(x_true, y_true_noisy, yerr=err, fmt='.k', capsize=0)

# try a traditional optimizer

def log_likelihood(theta, x, y, yerr, t_infl):
    (xt, model) = snake.snake(x, t_infl, theta)
    sigma2 = yerr ** 2
    return  -0.5 * np.sum((y - model)**2 / sigma2 + np.log(sigma2))


nll = lambda *args: -log_likelihood(*args)

from scipy.optimize import minimize

initial = f_inflexion + 0.1 * np.random.randn(f_inflexion.size)

soln = minimize(nll, initial, args=(x_true, y_true_noisy, err, t_inflexion))

p = soln.x
print('initial is ',initial)
print('result is  ',p)

print(soln)

(x_fit,y_fit) = snake.snake(t, t_inflexion, p)

ax1.plot(x_fit,y_fit, 'r-')

import emcee

def log_prior(theta):
    if np.all(theta>12) and np.all(theta<17):
        return 0.0
#    return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr, t_infl):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr, t_infl)

pos = soln.x + 1e-2 * np.random.randn(32, f_inflexion.size)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x_true, y_true_noisy, err, t_inflexion))
sampler.run_mcmc(pos, 5000, progress=True);

fig2, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(00, len(samples))

axes[-1].set_xlabel("step number");

# get randomized samples to plot up
flat_samples = sampler.get_chain(discard=1000, thin=15, flat=True)

import corner
fig4, ax4 = plt.subplots(1, 1, figsize=(8, 8))


fig = corner.corner(
    flat_samples, truths=f_inflexion
);

inds = np.random.randint(len(flat_samples), size=100)
for ind in inds:
    sample = flat_samples[ind]
    (x_fit,y_fit) = snake.snake(t, t_inflexion, sample)
    ax1.plot(x_fit,y_fit, "C1", alpha=0.1)



plt.draw()
plt.show()

import numpy as np

def snake(x, x_infl, y_infl):
    '''for points along x, with ends of straight lines marked by
     (x_infl, y_infl), return (xt, y) points
    where xt are points in x that are inside the ends of the straight lines
    x_infl must be strictly ordered in increasing values
    '''

    # reject points lower than min(x_infl) and greater than max(x_infl)
    xt = x[np.where( (x>=np.min(x_infl)) * (x<=np.max(x_infl)))]

    # calculate the constants for y = mx + c
    m = (y_infl[1:] - y_infl[:-1]) / (x_infl[1:] - x_infl[:-1])

    # now we have m...
    # c = y - mx
    c = y_infl[:-1] - m*x_infl[:-1]

    # now for each point in xt we want to know what interval it lands in.
    # lazy way of doing this is

    # BROADCAST to see where each input x coordinate is less than a given inflexion point, resulting in a 2D array
    n_positive = xt - x_infl[:,np.newaxis]

    # you can then count how many times a given x point is less than the list of inflexion points!
    count_positive = np.sum((n_positive>0),axis=0)

    #... so then the index for m and c is then this count minus one.
    ind = count_positive-1

    # calculate 
    y = m[ind]*xt + c[ind]

    return (xt, y)


def log_likelihood(theta, x, y, yerr, t_infl):
    (xt, model) = snake(x, t_infl, theta)
    sigma2 = yerr ** 2
    return  -0.5 * np.sum((y - model)**2 / sigma2 + np.log(sigma2))


def snakemcee(x, y, y_err, xmodel, ymodel):
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

    # default is 32 walkers and 1e-2 error on their starting position
    pos = ymodel + 1e-2 * np.random.randn(32, ymodel.size)
    nwalkers, ndim = pos.shape

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, y_err, xmodel))
    sampler.run_mcmc(pos, 5000, progress=True);

    return sampler


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    # make a test data set

    time = np.linspace(58000, 58099.5, 201)

    t_inflexion = np.linspace(58000,58100, 20)
    t_inflexion = t_inflexion + np.random.uniform(0,5,size=20)

    f_min = 14

    f_max = 15

    m_inflexion = np.random.uniform(f_min, f_max, size=np.shape(t_inflexion))


    fig, ax = plt.subplots(1, 1, figsize=(12, 4)) 

    ax.set_xlim(np.min(t_inflexion),np.max(t_inflexion))

    ax.vlines(t_inflexion,f_min,f_max,linestyle='dashed')
    ax.scatter(t_inflexion, m_inflexion,color='orange')

    (x,y) = snake(time, t_inflexion, m_inflexion)
    ax.plot(x,y) 
    
    plt.draw()
    plt.show()


    # testing an emcee fitting to the lines

    # make a test data set

    N = 101
    n_inflex = 12
    t = np.linspace(58000, 58099.5, N)
    np.random.seed(77)
    
    t_inflexion = np.linspace(58000,58100, n_inflex)
    t_inflexion += np.random.uniform(0,5,size=t_inflexion.size)

    f_min = 14

    f_max = 15

    f_inflexion = np.random.uniform(f_min, f_max, size=np.shape(t_inflexion))

    fig1, ax1 = plt.subplots(1, 1, figsize=(12, 4))

    ax1.set_xlim(np.min(t_inflexion),np.max(t_inflexion))

    ax1.vlines(t_inflexion, f_min, f_max, linestyle='dashed')

    ax1.scatter(t_inflexion, f_inflexion,color='orange')


    (x_true,y_true) = snake(t, t_inflexion, f_inflexion)

    err = 0.1 * np.random.randn(y_true.size)

    y_true_noisy = y_true + np.random.randn(y_true.size)*err

    ax1.plot(x_true, y_true, 'k', alpha=0.3, lw=5)

    ax1.errorbar(x_true, y_true_noisy, yerr=err, fmt='.k', capsize=0)

    # try a traditional optimizer

    nll = lambda *args: -log_likelihood(*args)

    from scipy.optimize import minimize

    initial = f_inflexion + 0.1 * np.random.randn(f_inflexion.size)

    soln = minimize(nll, initial, args=(x_true, y_true_noisy, err, t_inflexion))

    p = soln.x
    print('initial is ',initial)
    print('result is  ',p)

    print(soln)

    (x_fit,y_fit) = snake(t, t_inflexion, p)

    ax1.plot(x_fit,y_fit, 'r-')

    sampler = snakemcee(x_true, y_true_noisy, err, t_inflexion, p)

    fig2, axes = plt.subplots(n_inflex, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    for i in range(n_inflex):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(00, len(samples))

    axes[-1].set_xlabel("step number");

    # get randomized samples to plot up
    flat_samples = sampler.get_chain(discard=1000, thin=15, flat=True)

    import corner
    fig = corner.corner(
        flat_samples, truths=f_inflexion
    );

    inds = np.random.randint(len(flat_samples), size=100)
    for ind in inds:
        sample = flat_samples[ind]
        (x_fit,y_fit) = snake(t, t_inflexion, sample)
        ax1.plot(x_fit,y_fit, "C1", alpha=0.1)


    plt.draw()
    plt.show()





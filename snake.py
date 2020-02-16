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

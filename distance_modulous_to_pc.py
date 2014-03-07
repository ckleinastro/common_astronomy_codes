
def scatter2fract_dist(scatter_mag):
    upper = ((10**((scatter_mag+5)/5))/10)-1
    lower = 1-((10**((-scatter_mag+5)/5))/10)
    return (upper + lower)/2
    

mu = 15.0 # mu is simply m-M, without accounting for exctinction A
mu_err = 0.014
A = 0.0
A_err = 0.0

effective_mu = mu-A
effective_mu_err = (mu_err**2 + A_err**2)**0.5

dist = 10**(effective_mu/5. + 1)
dist_fract_err = scatter2fract_dist(effective_mu_err)
dist_err = dist_fract_err * dist

print "Distance is:", dist, "+/-", dist_err, " pc, or", round(100*dist_err/dist, 3), "%"
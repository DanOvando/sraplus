r = 0.2

m = .04

k = 10000

b0 = 10

years = 1000

plim = 0.1

msy = r*k / 4

test = testpopmodel(r = r, k = k, m = m, b0 = b0, plim = plim, years = years, sigma_proc = 0, catches = rep(10, years))

plot(test, lead(test) - test)

plot(test)

NumericVector popmodel(double r, double k, double m,double b0,double plim,int years, double sigma_proc,NumericVector catches){
  
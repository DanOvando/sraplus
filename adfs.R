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

years = 10

msy = rk/4

test <-
  graddesc(
    r = 0.4,
    m = 2,
    init_deps = 1,
    plim = 0.05,
    years = years,
    sigma_procs = 0,
    catches = rep(0, years),
    final_state = 0.5,
    learn_rate = .1
  )

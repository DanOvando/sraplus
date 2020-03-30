

catch_only_driors <- sraplus::format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  use_heuristics = FALSE,
  terminal_state = .25,
  terminal_state_cv = .5
)


catch_only_fit <- fit_sraplus(driors = catch_only_driors,
                              engine = "sir",
                              draws = 1e5,
                              n_keep = 2000,
                              estimate_proc_error = FALSE, 
                              estimate_shape = TRUE,
                              tune_prior_predictive = TRUE)


plot_sraplus(catch_only_fit)


plot_prior_posterior(catch_only_fit, catch_only_driors)


keepers <- sra_fit$keepers

outs <- stringr::str_detect(names(sra_fit), "_t")

sra_fit$b_t[, keepers] -> a

wtf <- sra_fit$dep_t[nrow(sra_fit$dep_t), ] -> a

state_breaks <- seq(0,2, by = .05)

state_bins <- cut(state_breaks, state_breaks, include.lowest = FALSE, right = FALSE)

edge_p <-  pnorm(log(state_breaks), log(0.5),.2)

p_bin <- lead(edge_p) - (edge_p)

bin_frame <- data.frame(bin = state_bins, p_bin = p_bin) %>% 
  mutate(bin = as.character(bin))


draws <- data.frame(state = wtf) %>% 
  mutate(bin = as.character(cut(state, state_breaks,include.lowest = FALSE, right = FALSE))) %>% 
  left_join(bin_frame, by = "bin") %>% 
  group_by(bin) %>% 
  mutate(weight = unique(p_bin) / length(p_bin)) %>% 
    ungroup() %>% 
  mutate(index = 1:nrow(.))

sample_index <- sample(draws$index, length(keepers), replace = TRUE, prob = draws$weight)

hmm <- draws$state[sample_index]

hist(hmm)

ggplot() + 
  geom_density(data = draws, aes(state)) + 
  geom_density(data = data.frame(a = rlnorm(1000, log(0.5),.2)), aes(a)) + 
  geom_density(data = data.frame(a = hmm), aes(a)) 

hist(sra_fit$k[sample_index])
  



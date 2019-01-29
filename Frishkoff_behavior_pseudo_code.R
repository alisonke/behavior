# Y is the observed number of animals at each site (j, from 1 to J), during each visit (k, from 1 to K), doing each behavior (b, from 1 to B)
# D is the inferred number of individuals *D*oing a behavior during a given site visit

Y[j,k,b] ~ dbin(P[j,k,b], D[j,k, b])

# Allow behavior to modified det. prob. -- may or may not work due to identifiability issues
logit(P[j,k,b]) <- alpha.0[b] + alpha.1*X[j] + ...

alpha.0[b] ~ dnorm(mu.alpha, tau.alpha) # This I imagine would be critical to ensure identifiability

alpha.1 ~ ...

# Or assume that behavior (or the subset of behaviors that we're querying), do not affect detection to avoid potential identifiability issues. My take is that this is an empirical question best answered through simulation.
logit(P[j,k,b]) <- alpha.0 + alpha.1*X[j] + ...


# A multinomial distribution (or something similar) is used to divide up the total number of animals during a sampling occurence. This I think allows the assumption that the underlying number of indiviuals is the same across all site visits is the same (ensuring closure), but allows the actual number of animals doing each behavior to vary (but maintains identifiability --- I THINK --- because the underlying probability of each behavior is assumed to be the same from visit to visit.). I again am not sure this will work, in theory or practice, but whether or not it does could be tested with simulations (and whether the degree to which it is wrong if it doesn't work is material or not).
for (j in 1:J){
	for (k in 1:K){	
		D[j,k,1:B] ~ dmulti(pi[1:b], N[j])
	}
}
# pi is the underlying probabilities of each behavior. I'm pretty sure that multinomial will internally scale pi to sum to one, so I think the behaviors could be parameterized something like:

pi[1] <- 1
for (2:B) {
	pi[b] ~ dnorm(0, 0.001)
}
# The one, issue with the multinomial I think is the N has to be at least 1. I'm sure there's some clever way to get around this though, for instances when no birds exist at a site...

# Then N is determined as normal:
N[j] ~ dpois(lambda[j])
log(lambda[j]) <- betas....


## My thinking is that this formulation allows N not be a derived quantity, and permits the behaviors to stochasitcally vary between visits (while still conforming to N, and benefitting from being constrained by the assumption that the underlying probabilities are constant). I could easily be missing something that would make this whole approach not work, and I can imagine there could be certain areas of parameter space where this model is unidentifiable with realistic amounts of data. But at least right now my take is that it would be worth building some simulations and seeing how well it does. If think all the assumptions check out as formulated... though I could certainly be wrong there. If we build some simulations and it doesn't work, we can scrap the idea with relativly little time wasted. If it does work than that's great! Those are my two cents anyway...



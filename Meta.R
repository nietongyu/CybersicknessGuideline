library(metafor)

# Full Data
data = read.csv("Results.csv")

res_mod_all <- rma(yi = Effect.Size, vi = Variance, model="RE", data = data)
print(res_mod_all$b)
print(confint(res_mod_all))
summary(res_mod_all)

### calculate influence diagnostics
inf <- influence(res_mod_all)
summary(inf)

### plot the influence diagnostics
par(mfrow=c(8,1))
plot(inf)

data = data[-67, ] # upon inspecting inf, we found that row 67 is an outlier

res_mod_all <- rma(yi = Effect.Size, vi = Variance, model="RE", data = data)
print(res_mod_all$b)
print(confint(res_mod_all))
summary(res_mod_all)

res_mod <- rma(yi = Effect.Size, vi = Variance, mods = ~ Design + Female.Participant.Proportion, data = data)
summary(res_mod)

# Delta-SSQ only

data = read.csv("Results_delta_only.csv")

res_mod_delta <- rma(yi = Effect.Size, vi = Variance, model="RE", data = data)
print(res_mod_delta$b)
print(confint(res_mod_delta))
summary(res_mod_delta)

### calculate influence diagnostics
inf <- influence(res_mod_delta)
summary(inf)

### plot the influence diagnostics
par(mfrow=c(8,1))
plot(inf)

# upon inspecting inf, we did not find any outlier

res_mod_delta <- rma(yi = Effect.Size, vi = Variance, model="RE", data = data)
print(res_mod_delta$b)
print(confint(res_mod_delta))
summary(res_mod_delta)

res_mod <- rma(yi = Effect.Size, vi = Variance, mods = ~ Design + Female.Participant.Proportion, data = data)
summary(res_mod)










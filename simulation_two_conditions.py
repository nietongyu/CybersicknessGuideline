import numpy as np

number_of_conds = 2
mean_1 = 50
mean_2 = 30
carryover_mean = .53 # 0.63# the amount of previous session that got carried into the result
carryover_dispersion = 0.001 #* carryover_mean
# mean_1 = 50
# mean_2 = 30
# carryover = -0.21 # the amount of previous session that got carried into the result

sample_size = 30 # number of participants
j = 1 - 3 / (4 * (sample_size - 1) -1) # Hedges’ g

d = 1.43 # 1.28 # cohen's d
sd_1 = mean_1 / d
sd_2 = mean_2 / d


dampen = 0.3
# SDs with carryover
sd_1_carry =  (mean_1 + carryover_mean * mean_2 * (1 / number_of_conds) * dampen) / d
sd_2_carry =  (mean_2 + carryover_mean * mean_1 * (1 / number_of_conds) * dampen) / d


# record statistical significance
total_simulation = 10000
alpha = 0.05
icc = 0.6 #0.7
ind_sig = 0
pair_sig = 0
pair_carry_sig = 0
total_icc_with_carry = 0
total_icc_without_carry = 0
total_effect_size_without_carry = 0
total_effect_size_with_carry = 0



import numpy as np
from scipy.stats import t


# ---------- Helpers ----------

def rowwise_ttest_ind(x, y):
    nx = x.shape[1]
    ny = y.shape[1]

    mx = x.mean(axis=1)
    my = y.mean(axis=1)

    vx = x.var(axis=1, ddof=1)
    vy = y.var(axis=1, ddof=1)

    # Welch's t-test
    t_stat = (mx - my) / np.sqrt(vx/nx + vy/ny)

    df = (vx/nx + vy/ny)**2 / ((vx**2)/((nx**2)*(nx-1)) + (vy**2)/((ny**2)*(ny-1)))

    p = 2 * t.sf(np.abs(t_stat), df)
    return p


def rowwise_ttest_rel(x, y):
    d = x - y
    md = d.mean(axis=1)
    sd = d.std(axis=1, ddof=1)

    t_stat = md / (sd / np.sqrt(d.shape[1]))
    df = d.shape[1] - 1

    p = 2 * t.sf(np.abs(t_stat), df)
    return p


def rowwise_corr(x, y):
    x_center = x - x.mean(axis=1, keepdims=True)
    y_center = y - y.mean(axis=1, keepdims=True)

    num = np.sum(x_center * y_center, axis=1)
    den = np.sqrt(np.sum(x_center**2, axis=1) * np.sum(y_center**2, axis=1))

    return num / den


def effect_size(x):
    return x.mean(axis=1) / x.std(axis=1, ddof=1)


# ---------- Simulation ----------

def run_simulation_vectorized():
    N = total_simulation
    S = sample_size

    # =========================
    # Independent Samples
    # =========================
    s1 = np.random.normal(mean_1, sd_1, size=(N, S))
    s2 = np.random.normal(mean_2, sd_2, size=(N, S))

    p_ind = rowwise_ttest_ind(s1, s2)
    ind_sig = np.sum(p_ind < alpha)

    # =========================
    # Paired (no carryover)
    # =========================
    cov = [[sd_1**2, icc * sd_1 * sd_2],
           [icc * sd_1 * sd_2, sd_2**2]]

    paired = np.random.multivariate_normal(
        [mean_1, mean_2], cov, size=(N, S)
    )

    s1 = paired[:, :, 0]
    s2 = paired[:, :, 1]

    p_pair = rowwise_ttest_rel(s1, s2)
    pair_sig = np.sum(p_pair < alpha)

    eff_no_carry = np.mean((effect_size(s1) + effect_size(s2)) / 2)
    icc_no_carry = np.mean(rowwise_corr(s1, s2))

    # =========================
    # Paired (with carryover)
    # =========================
    cov_carry = [[sd_1_carry**2, icc * sd_1_carry * sd_2_carry],
                 [icc * sd_1_carry * sd_2_carry, sd_2_carry**2]]

    paired = np.random.multivariate_normal(
        [mean_1, mean_2], cov_carry, size=(N, S)
    )

    half = S // 2

    s1_1 = paired[:, :half, 0]
    s1_2 = paired[:, :half, 1]
    s2_1 = paired[:, half:, 0]
    s2_2 = paired[:, half:, 1]

    # Carryover noise (vectorized)
    c1 = np.random.normal(1, carryover_dispersion, size=(N, half)) * carryover_mean
    c2 = np.random.normal(1, carryover_dispersion, size=(N, half)) * carryover_mean

    s1_2 = s1_2 + c1 * s1_1
    s2_1 = s2_1 + c2 * s2_2

    s1_carry = np.concatenate([s1_1, s2_1], axis=1)
    s2_carry = np.concatenate([s1_2, s2_2], axis=1)

    p_carry = rowwise_ttest_rel(s1_carry, s2_carry)
    pair_carry_sig = np.sum(p_carry < alpha)

    eff_carry = np.mean((effect_size(s1_carry) + effect_size(s2_carry)) / 2)
    icc_carry = np.mean(rowwise_corr(s1_carry, s2_carry))

    return {
        "ind_sig": ind_sig,
        "pair_sig": pair_sig,
        "pair_carry_sig": pair_carry_sig,
        "effect_no_carry": eff_no_carry * j,
        "effect_carry": eff_carry * j,
        "icc_no_carry": icc_no_carry,
        "icc_carry": icc_carry
    }

results = run_simulation_vectorized()
print(f'Independent significant: {results["ind_sig"]/total_simulation:.2f}')
print(f'Paired significant: {results["pair_sig"]/total_simulation:.2f}')
print(f'Paired (carryover) significant: {results["pair_carry_sig"]/total_simulation:.2f}')

print("Effect size (no carry):", results["effect_no_carry"])
print("Effect size (carry):", results["effect_carry"])

print("ICC (no carry):", results["icc_no_carry"])
print("ICC (carry):", results["icc_carry"])
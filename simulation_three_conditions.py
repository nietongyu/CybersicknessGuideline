import numpy as np
from scipy.stats import f, t

# ---------- Parameters ----------

number_of_conds = 3
means = np.array([20, 40, 60])  # True means for the 3 conditions
carryover_mean = 0.53
carryover_dispersion = 0.001 
sample_size = 30 # Number of participants (divisible by 6 for balanced square)

d = 1.43 
sds = means / d

dampen = 0.3
# Inflate SDs for the carryover generation based on your original logic
# Using the average of the "other" two conditions to calculate the inflation
sds_carry = []
for i in range(number_of_conds):
    other_means_avg = (np.sum(means) - means[i]) / (number_of_conds - 1)
    effect_to_carry = other_means_avg * (1 / number_of_conds) * carryover_mean
    inflated_sd = (means[i] + (2 * effect_to_carry) * dampen) / d
    sds_carry.append(inflated_sd)
sds_carry = np.array(sds_carry)

# Record statistical significance
total_simulation = 10000
alpha = 0.05
icc = 0.6 
j = 1 - 3 / (4 * (sample_size - 1) - 1)

# ---------- Helpers ----------

def rowwise_rm_anova(*conds):
    """
    Vectorized One-Way Repeated Measures ANOVA across N simulations.
    Accepts arbitrary number of condition arrays, each of shape (Simulations, Subjects).
    """
    k = len(conds)
    n = conds[0].shape[1]
    
    # Stack into shape (Simulations, Subjects, Conditions)
    X = np.stack(conds, axis=2)
    
    grand_mean = X.mean(axis=(1, 2), keepdims=True)
    sub_mean = X.mean(axis=2, keepdims=True)
    cond_mean = X.mean(axis=1, keepdims=True)
    
    # Sum of Squares
    SS_total = np.sum((X - grand_mean)**2, axis=(1, 2))
    SS_sub = k * np.sum((sub_mean - grand_mean)**2, axis=(1, 2))
    SS_cond = n * np.sum((cond_mean - grand_mean)**2, axis=(1, 2))
    SS_err = SS_total - SS_sub - SS_cond
    
    # Degrees of Freedom
    df_cond = k - 1
    df_err = (n - 1) * (k - 1)
    
    # Mean Squares & F-statistic
    MS_cond = SS_cond / df_cond
    MS_err = SS_err / df_err
    F_stat = MS_cond / MS_err
    
    # p-values from survival function
    p = f.sf(F_stat, df_cond, df_err)
    return p

def effect_size(x):
    # Standard score relative to 0 (Pre-experiment baseline assumption)
    return x.mean(axis=1) / x.std(axis=1, ddof=1)

def rowwise_corr_3x(c1, c2, c3):
    """Average pairwise Pearson correlation for 3 variables."""
    def corr(x, y):
        xc = x - x.mean(axis=1, keepdims=True)
        yc = y - y.mean(axis=1, keepdims=True)
        return np.sum(xc * yc, axis=1) / np.sqrt(np.sum(xc**2, axis=1) * np.sum(yc**2, axis=1))
    
    r12 = corr(c1, c2)
    r13 = corr(c1, c3)
    r23 = corr(c2, c3)
    return (r12 + r13 + r23) / 3

# ---------- Simulation ----------

def run_simulation_3_conditions():
    N = total_simulation
    S = sample_size
    
    # =========================
    # Paired (no carryover)
    # =========================
    cov = np.zeros((3, 3))
    for i in range(3):
        for k in range(3):
            cov[i, k] = sds[i]**2 if i == k else icc * sds[i] * sds[k]

    paired = np.random.multivariate_normal(means, cov, size=(N, S))
    
    s1, s2, s3 = paired[:, :, 0], paired[:, :, 1], paired[:, :, 2]
    
    p_pair = rowwise_rm_anova(s1, s2, s3)
    pair_sig = np.sum(p_pair < alpha)
    
    eff_no_carry = np.mean((effect_size(s1) + effect_size(s2) + effect_size(s3)) / 3)
    icc_no_carry = np.mean(rowwise_corr_3x(s1, s2, s3))

    # =========================
    # Paired (with carryover)
    # =========================
    cov_carry = np.zeros((3, 3))
    for i in range(3):
        for k in range(3):
            cov_carry[i, k] = sds_carry[i]**2 if i == k else icc * sds_carry[i] * sds_carry[k]

    paired_carry = np.random.multivariate_normal(means, cov_carry, size=(N, S))
    
    # 6 Sequences for a perfectly balanced 3-condition Latin Square (Williams Design)
    # 0=Cond1, 1=Cond2, 2=Cond3
    sequences = [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
        [0, 2, 1],
        [1, 0, 2],
        [2, 1, 0]
    ]
    
    n_per_seq = S // 6  # 5 participants per sequence
    
    # Receptacles for final sorted scores by Condition
    cond_scores = [np.zeros((N, S)) for _ in range(3)]
    
    for g, seq in enumerate(sequences):
        i1, i2, i3 = seq
        s_idx = slice(g * n_per_seq, (g + 1) * n_per_seq)
        
        # Extract base scores for this group in the order they experience them
        p1 = paired_carry[:, s_idx, i1]
        p2 = paired_carry[:, s_idx, i2]
        p3 = paired_carry[:, s_idx, i3]
        
        # Generate carryover noise
        c1 = np.random.normal(1, carryover_dispersion, size=(N, n_per_seq)) * carryover_mean
        c2 = np.random.normal(1, carryover_dispersion, size=(N, n_per_seq)) * carryover_mean
        
        # Apply sequential carryover
        final_p1 = p1
        final_p2 = p2 + c1 * final_p1
        final_p3 = p3 + c2 * final_p2 + c1 * p1
        
        # Map them back to the correct Condition arrays
        cond_scores[i1][:, s_idx] = final_p1
        cond_scores[i2][:, s_idx] = final_p2
        cond_scores[i3][:, s_idx] = final_p3

    s1_carry, s2_carry, s3_carry = cond_scores
    
    p_carry = rowwise_rm_anova(s1_carry, s2_carry, s3_carry)
    pair_carry_sig = np.sum(p_carry < alpha)
    
    eff_carry = np.mean((effect_size(s1_carry) + effect_size(s2_carry) + effect_size(s3_carry)) / 3)
    icc_carry = np.mean(rowwise_corr_3x(s1_carry, s2_carry, s3_carry))

    return {
        "pair_sig": pair_sig,
        "pair_carry_sig": pair_carry_sig,
        "effect_no_carry": eff_no_carry * j,
        "effect_carry": eff_carry * j,
        "icc_no_carry": icc_no_carry,
        "icc_carry": icc_carry
    }

results = run_simulation_3_conditions()
print(f'Repeated Measures ANOVA significant: {results["pair_sig"]/total_simulation:.2f}')
print(f'Repeated Measures ANOVA (carryover) significant: {results["pair_carry_sig"]/total_simulation:.2f}')
print("-" * 30)
print(f'Avg Effect size (no carry): {results["effect_no_carry"]:.4f}')
print(f'Avg Effect size (carry): {results["effect_carry"]:.4f}')
print(f'ICC (no carry): {results["icc_no_carry"]:.4f}')
print(f'ICC (carry): {results["icc_carry"]:.4f}')
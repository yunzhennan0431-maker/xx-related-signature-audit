import pandas as pd
import numpy as np

main = pd.read_csv('analysis_output/data/phase3_meta_analysis_main.csv')

THEME_GROUPS = {
    'CellDeath': ['ferroptosis','pyroptosis','cuproptosis','disulfidptosis','necroptosis','autophagy','anoikis','PCDgenes','ICD','mitophagyE3'],
    'Metabolism': ['glycolysis','glutamine','cholesterol','mitochondrial','lipid','ERstress','aminoacid','ATIC'],
    'RNAmod': ['m6A'],
    'Chromatin': ['chromatin'],
    'HypoxiaOxidative': ['hypoxia','oxidative_stress'],
    'Senescence': ['senescence'],
    'Exosome': ['exosome'],
    'ImmuneTME': ['immune','CAF'],
    'Lactylation': ['lactylation'],
}

def classify(model):
    for grp, kws in THEME_GROUPS.items():
        for kw in kws:
            if kw.lower() in model.lower():
                return grp
    return 'Other'

main['theme_group'] = main['model'].apply(classify)
main['sig'] = main['significant_fdr05'].astype(int)

groups = main['theme_group'].values
sig = main['sig'].values
n_total = len(main)
n_sig_total = sig.sum()

def chi_stat(sig_vec, groups_vec):
    stat = 0.0
    for g in np.unique(groups_vec):
        mask = groups_vec == g
        n_g = mask.sum()
        obs = sig_vec[mask].sum()
        exp = n_g * n_sig_total / n_total
        if exp > 0:
            stat += (obs - exp) ** 2 / exp
    return stat

observed_stat = chi_stat(sig, groups)
print("Observed heterogeneity statistic (chi-square-like, by theme group):", observed_stat)

rng = np.random.default_rng(20260906)
n_perm = 20000
perm_stats = np.zeros(n_perm)
for i in range(n_perm):
    perm_sig = rng.permutation(sig)
    perm_stats[i] = chi_stat(perm_sig, groups)

p_value = (perm_stats >= observed_stat).mean()
print(f"Permutation p-value (H0: significance independent of theme group): {p_value:.4f}")
print(f"Permutation null distribution: mean={perm_stats.mean():.2f}, 95th pct={np.percentile(perm_stats,95):.2f}, 99th pct={np.percentile(perm_stats,99):.2f}")

print("\nObserved per-group significant rate vs overall:")
overall_rate = n_sig_total / n_total
summary = main.groupby('theme_group')['sig'].agg(['sum','count'])
summary['rate'] = summary['sum'] / summary['count']
summary['overall_rate'] = overall_rate
print(summary.sort_values('rate', ascending=False))

"""Generate all documentation figures from sample data."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pybaseflow
from pybaseflow.separation import (
    lh, lh_multi, chapman, chapman_maxwell, boughton, eckhardt,
    ewma, furey, willems, ihacres, what,
    fixed, slide, local, ukih, part,
    strict_baseflow, _recursive_digital_filter,
)
from pybaseflow.estimate import recession_coefficient, bflow, bflow_recession_analysis
from pybaseflow.tracer import cmb

FIGDIR = "docs/assets/figures"
data = pybaseflow.load_sample_data()
Q = data['Q']
dates = data['dates']

# Estimate recession coefficient for methods that need it
strict = strict_baseflow(Q)
a = recession_coefficient(Q, strict)

# Pre-compute LH baseflow (needed by local, ukih)
b_LH = lh(Q)

# Shared plot style
plt.rcParams.update({
    'figure.figsize': (12, 4),
    'figure.dpi': 150,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'font.size': 10,
})

def save(fig, name):
    fig.savefig(f"{FIGDIR}/{name}.png", bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  {name}.png")


def plot_separation(dates, Q, b, method_name, color='#2196F3', fname=None):
    """Standard baseflow separation plot."""
    fig, ax = plt.subplots()
    ax.plot(dates, Q, color='#333333', linewidth=0.6, label='Streamflow')
    ax.fill_between(dates, 0, b, alpha=0.5, color=color, label=f'{method_name} baseflow')
    ax.fill_between(dates, b, Q, alpha=0.15, color='#FF5722', label='Quick flow')
    ax.set_ylabel('Discharge (ft$^3$/s)')
    ax.set_xlabel('')
    ax.legend(loc='upper right', framealpha=0.9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    fig.autofmt_xdate()
    bfi = np.sum(b) / np.sum(Q)
    ax.set_title(f'{method_name} (BFI = {bfi:.3f})')
    save(fig, fname or method_name.lower().replace(' ', '_').replace('-', '_'))


# =========================================================================
# 1. Overview: all methods comparison
# =========================================================================
print("Generating overview figures...")

fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(dates, Q, color='#333333', linewidth=0.8, label='Streamflow', zorder=10)

methods = {
    'Eckhardt': eckhardt(Q, a, BFImax=0.8),
    'Lyne-Hollick': b_LH,
    'Chapman-Maxwell': chapman_maxwell(Q, a),
    'EWMA': ewma(Q, e=0.01),
    'PART': part(Q, area=1600),
    'Fixed interval': fixed(Q, area=1600),
}
colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0', '#F44336', '#795548']
for (name, b), c in zip(methods.items(), colors):
    bfi = np.sum(b) / np.sum(Q)
    ax.plot(dates, b, linewidth=0.9, color=c, label=f'{name} ({bfi:.2f})', alpha=0.85)

ax.set_ylabel('Discharge (ft$^3$/s)')
ax.legend(loc='upper right', fontsize=8, framealpha=0.9, ncol=2)
ax.set_title('Comparison of baseflow separation methods — Fish River, Maine (2019–2020)')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
fig.autofmt_xdate()
save(fig, 'overview_comparison')

# =========================================================================
# 2. Individual method plots — digital filters
# =========================================================================
print("Generating digital filter figures...")

plot_separation(dates, Q, eckhardt(Q, a, BFImax=0.8), 'Eckhardt')
plot_separation(dates, Q, chapman_maxwell(Q, a), 'Chapman-Maxwell', '#4CAF50')
plot_separation(dates, Q, chapman(Q, a), 'Chapman', '#009688')
plot_separation(dates, Q, boughton(Q, a, C=0.05), 'Boughton', '#FF9800')
plot_separation(dates, Q, ewma(Q, e=0.01), 'EWMA', '#9C27B0')
plot_separation(dates, Q, furey(Q, a, A=0.5), 'Furey-Gupta', '#E91E63', 'furey')
plot_separation(dates, Q, willems(Q, a, w=0.5), 'Willems', '#00BCD4')
plot_separation(dates, Q, ihacres(Q, a, C=0.3, alpha_s=-0.5), 'IHACRES', '#3F51B5')

# Lyne-Hollick with pass comparison
print("Generating Lyne-Hollick multi-pass figure...")
fig, ax = plt.subplots()
ax.plot(dates, Q, color='#333333', linewidth=0.6, label='Streamflow')
for n, c, ls in [(1, '#F44336', '--'), (2, '#2196F3', '-'), (3, '#4CAF50', '-.')]:
    b = lh_multi(Q, num_pass=n)
    bfi = np.sum(b) / np.sum(Q)
    ax.plot(dates, b, color=c, linewidth=0.9, linestyle=ls,
            label=f'{n}-pass (BFI = {bfi:.3f})')
ax.set_ylabel('Discharge (ft$^3$/s)')
ax.legend(loc='upper right', framealpha=0.9)
ax.set_title('Lyne-Hollick filter — effect of number of passes')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
fig.autofmt_xdate()
save(fig, 'lyne_hollick_passes')

# =========================================================================
# 3. Filter family comparison (gamma=0 vs gamma=1)
# =========================================================================
print("Generating filter family figure...")

# Zoom into a 90-day window around a peak for clarity
peak_idx = np.argmax(Q)
start = max(0, peak_idx - 30)
end = min(len(Q), peak_idx + 60)
d_zoom = dates[start:end]
Q_zoom = Q[start:end]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.5), sharey=True)

# gamma=0 family
ax1.plot(d_zoom, Q_zoom, color='#333333', linewidth=0.8, label='Streamflow')
for name, b_full, c in [
    ('Chapman-Maxwell', chapman_maxwell(Q, a), '#2196F3'),
    ('Eckhardt', eckhardt(Q, a, BFImax=0.8), '#4CAF50'),
    ('Boughton', boughton(Q, a, C=0.05), '#FF9800'),
    ('EWMA', ewma(Q, e=0.01), '#9C27B0'),
]:
    ax1.plot(d_zoom, b_full[start:end], linewidth=1.0, color=c, label=name)
ax1.set_title(r'$\gamma = 0$ family (linear reservoir)')
ax1.set_ylabel('Discharge (ft$^3$/s)')
ax1.legend(fontsize=8, framealpha=0.9)
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

# gamma=1 family
ax2.plot(d_zoom, Q_zoom, color='#333333', linewidth=0.8, label='Streamflow')
for name, b_full, c in [
    ('Lyne-Hollick', b_LH, '#2196F3'),
    ('Chapman (1991)', chapman(Q, a), '#4CAF50'),
    ('Willems', willems(Q, a, w=0.5), '#FF9800'),
]:
    ax2.plot(d_zoom, b_full[start:end], linewidth=1.0, color=c, label=name)
ax2.set_title(r'$\gamma = 1$ family (signal processing)')
ax2.legend(fontsize=8, framealpha=0.9)
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

fig.suptitle('Recursive digital filter families — peak event detail', fontsize=12, y=1.02)
fig.tight_layout()
save(fig, 'filter_families')

# =========================================================================
# 4. Graphical methods
# =========================================================================
print("Generating graphical method figures...")

plot_separation(dates, Q, fixed(Q, area=1600), 'Fixed Interval', '#795548')
plot_separation(dates, Q, slide(Q, area=1600), 'Sliding Interval', '#607D8B')
plot_separation(dates, Q, local(Q, b_LH, area=1600), 'Local Minimum', '#FF5722')
plot_separation(dates, Q, ukih(Q, b_LH), 'UKIH', '#009688')

# =========================================================================
# 5. PART method
# =========================================================================
print("Generating PART figure...")
plot_separation(dates, Q, part(Q, area=1600), 'PART', '#F44336')

# =========================================================================
# 6. BFlow multi-pass
# =========================================================================
print("Generating BFlow figure...")
result = bflow(Q)
fig, ax = plt.subplots()
ax.plot(dates, Q, color='#333333', linewidth=0.6, label='Streamflow')

b1 = lh_multi(Q, num_pass=1)
b2 = lh_multi(Q, num_pass=2)
b3 = result['baseflow']
ax.fill_between(dates, 0, b3, alpha=0.4, color='#2196F3', label=f'Pass 3 (BFI = {result["BFI"]:.3f})')
ax.plot(dates, b1, color='#F44336', linewidth=0.7, linestyle='--',
        label=f'Pass 1 (BFI = {result["BFI_pass1"]:.3f})')
ax.plot(dates, b2, color='#FF9800', linewidth=0.7, linestyle='-.',
        label=f'Pass 2 (BFI = {result["BFI_pass2"]:.3f})')

ax.set_ylabel('Discharge (ft$^3$/s)')
ax.legend(loc='upper right', framealpha=0.9)
ax.set_title(f'BFlow separation — alpha factor = {result["alpha_factor"]:.4f}, '
             f'baseflow days = {result["baseflow_days"]:.0f}')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
fig.autofmt_xdate()
save(fig, 'bflow')

# =========================================================================
# 7. CMB (synthetic specific conductance)
# =========================================================================
print("Generating CMB figure...")
np.random.seed(42)
bf_frac = b_LH / np.maximum(Q, 0.01)
bf_frac = np.clip(bf_frac, 0, 1)
SC_synth = bf_frac * 300 + (1 - bf_frac) * 50 + np.random.normal(0, 5, len(Q))
SC_synth = np.maximum(SC_synth, 10)

b_cmb = cmb(Q, SC_synth, SC_BF=300, SC_RO=50)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                                gridspec_kw={'height_ratios': [2, 1]})

ax1.plot(dates, Q, color='#333333', linewidth=0.6, label='Streamflow')
ax1.fill_between(dates, 0, b_cmb, alpha=0.5, color='#2196F3', label='CMB baseflow')
ax1.fill_between(dates, b_cmb, Q, alpha=0.15, color='#FF5722', label='Quick flow')
bfi = np.sum(b_cmb) / np.sum(Q)
ax1.set_title(f'Conductivity Mass Balance separation (BFI = {bfi:.3f})')
ax1.set_ylabel('Discharge (ft$^3$/s)')
ax1.legend(loc='upper right', framealpha=0.9)

ax2.plot(dates, SC_synth, color='#FF9800', linewidth=0.5)
ax2.axhline(300, color='#2196F3', linestyle='--', linewidth=0.8, label='SC$_{BF}$ = 300')
ax2.axhline(50, color='#F44336', linestyle='--', linewidth=0.8, label='SC$_{RO}$ = 50')
ax2.set_ylabel('Specific conductance\n($\\mu$S/cm)')
ax2.legend(loc='upper right', framealpha=0.9, fontsize=8)
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
fig.autofmt_xdate()
fig.tight_layout()
save(fig, 'cmb')

# =========================================================================
# 8. IHACRES alpha_s sensitivity
# =========================================================================
print("Generating IHACRES sensitivity figure...")

fig, ax = plt.subplots()
ax.plot(d_zoom, Q_zoom, color='#333333', linewidth=0.8, label='Streamflow')
for alpha_s, c, ls in [
    (0.0, '#FF9800', '--'),
    (-0.3, '#2196F3', '-'),
    (-0.6, '#4CAF50', '-.'),
    (-0.9, '#9C27B0', ':'),
]:
    b = ihacres(Q, a=0.98, C=0.3, alpha_s=alpha_s)
    bfi = np.sum(b) / np.sum(Q)
    label = f'$\\alpha_s$ = {alpha_s}' + (' (= Boughton)' if alpha_s == 0 else '')
    ax.plot(d_zoom, b[start:end], color=c, linewidth=1.0, linestyle=ls,
            label=f'{label} (BFI = {bfi:.3f})')
ax.set_ylabel('Discharge (ft$^3$/s)')
ax.set_title(r'IHACRES filter — effect of $\alpha_s$ parameter')
ax.legend(fontsize=8, framealpha=0.9)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
save(fig, 'ihacres_sensitivity')

# =========================================================================
# 9. Recession coefficient estimation
# =========================================================================
print("Generating recession coefficient figure...")

fig, ax = plt.subplots()
ax.plot(dates, Q, color='#333333', linewidth=0.6, label='Streamflow')
strict_mask = strict_baseflow(Q)
# Plot strict baseflow days as markers
strict_dates = dates[strict_mask]
strict_Q = Q[strict_mask]
ax.scatter(strict_dates, strict_Q, color='#F44336', s=4, zorder=5,
           label=f'Strict baseflow days ({strict_mask.sum()} of {len(Q)})')
ax.set_ylabel('Discharge (ft$^3$/s)')
ax.set_title(f'Strict baseflow identification — recession coefficient a = {a:.4f}')
ax.legend(loc='upper right', framealpha=0.9)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
fig.autofmt_xdate()
save(fig, 'recession_coefficient')

# =========================================================================
# 10. Brutsaert-Nieber (bn77) drought flow points
# =========================================================================
print("Generating bn77 figure...")

from pybaseflow.separation import bn77

drought_idx = bn77(Q, L_min=5, snow_freeze_period=(400, 450),
                   observational_precision=0.1, quantile=0.9)

fig, ax = plt.subplots()
ax.plot(dates, Q, color='#333333', linewidth=0.6, label='Streamflow')
if len(drought_idx) > 0:
    drought_dates = dates[drought_idx.astype(int)]
    drought_Q = Q[drought_idx.astype(int)]
    ax.scatter(drought_dates, drought_Q, color='#2196F3', s=12, zorder=5,
               edgecolors='#0D47A1', linewidths=0.3,
               label=f'Drought flow points ({len(drought_idx)})')
ax.set_ylabel('Discharge (ft$^3$/s)')
ax.set_title('Brutsaert-Nieber drought flow identification (bn77)')
ax.legend(loc='upper right', framealpha=0.9)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
fig.autofmt_xdate()
save(fig, 'bn77')

print("\nAll figures generated.")

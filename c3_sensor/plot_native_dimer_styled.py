import RNA
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import os
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# =====================================================================
# SEQUENCE & STRUCTURE PREDICTION
# =====================================================================
target_frag = 'UUCCCAGGCAAAAAACUAGUGCUGUCCAGUGAGAAGACUGUGCUGACCCCUGCCACCAACCACAUGGGCAACGUCACCUUCACGAUCCCAGCCAACAGGGAGUUCAAGUCAGAAAAGGGGCGCAACAAGUUCGUGACCGUGCAGGCCACCUUCGGGACCCAAGUGGUGGAGAAGGUGGUGCUGGUCAGCCUGCAGAGCGGGUA'
sensor = 'AUCGGGUCUGCGGGCUGAGCAGCACCGCCUUCGUCACCGCCUGGGUACGGAGGGAGACCUGCACGGUCACGACCCUGAUGCCCCCCUUUUGUGACGCGAUAGCGCUGUGUGCUGGGACCGAGAAGGCGACGGUGCCGAUGUGGUUGGUCGCUGGCCUCAGCAUCGCCUUCUCCCUGGACCGCUCACCUUUUUUUUCUGGGA'

complex_seq = f'{target_frag}&{sensor}'
struct, mfe = RNA.cofold(complex_seq)

# Get original coordinates from ViennaRNA layout engine
coords = RNA.get_xy_coordinates(struct)

xs = []
ys = []
for i in range(len(complex_seq)):
    if complex_seq[i] != '&':
        try:
            pt = coords.get(len(xs))
            xs.append(pt.X)
            ys.append(pt.Y)
        except:
            break

# Identify base pairs
pairs = []
stack = []
for i, char in enumerate(struct):
    if char == '(':
        stack.append(i)
    elif char == ')':
        if stack:
            j = stack.pop()
            pairs.append((j, i))

# Index mapping (skipping '&' to match coordinate array indices)
seq_no_amp = complex_seq.replace('&', '')
coord_map = {}
idx = 0
for i, char in enumerate(struct):
    if char != '&':
        coord_map[i] = idx
        idx += 1

target_len = len(target_frag)

# =====================================================================
# INTERPOLATION ALGORITHM TO DOUBLE THE DISTANCE BETWEEN STRANDS
# =====================================================================
expansion_factor = 2.0

S_x = [0.0] * len(xs)
S_y = [0.0] * len(ys)
is_paired = [False] * len(xs)

for (i, j) in pairs:
    c_i = coord_map[i]
    c_j = coord_map[j]
    dx = xs[c_i] - xs[c_j]
    dy = ys[c_i] - ys[c_j]
    S_x[c_i] = 0.5 * dx
    S_y[c_i] = 0.5 * dy
    is_paired[c_i] = True
    S_x[c_j] = -0.5 * dx
    S_y[c_j] = -0.5 * dy
    is_paired[c_j] = True

for n in range(target_len):
    if not is_paired[n]:
        L, R = -1, -1
        for s in range(n-1, -1, -1):
            if is_paired[s]: L = s; break
        for s in range(n+1, target_len):
            if is_paired[s]: R = s; break
        if L != -1 and R != -1:
            dL = n - L; dR = R - n
            S_x[n] = (S_x[L]*dR + S_x[R]*dL)/(dL+dR)
            S_y[n] = (S_y[L]*dR + S_y[R]*dL)/(dL+dR)
        elif L != -1:
            S_x[n] = S_x[L]; S_y[n] = S_y[L]
        elif R != -1:
            S_x[n] = S_x[R]; S_y[n] = S_y[R]

for n in range(target_len, len(xs)):
    if not is_paired[n]:
        L, R = -1, -1
        for s in range(n-1, target_len-1, -1):
            if is_paired[s]: L = s; break
        for s in range(n+1, len(xs)):
            if is_paired[s]: R = s; break
        if L != -1 and R != -1:
            dL = n - L; dR = R - n
            S_x[n] = (S_x[L]*dR + S_x[R]*dL)/(dL+dR)
            S_y[n] = (S_y[L]*dR + S_y[R]*dL)/(dL+dR)
        elif L != -1:
            S_x[n] = S_x[L]; S_y[n] = S_y[L]
        elif R != -1:
            S_x[n] = S_x[R]; S_y[n] = S_y[R]

for n in range(len(xs)):
    xs[n] += S_x[n] * expansion_factor
    ys[n] += S_y[n] * expansion_factor

# =====================================================================
# REFINED PLOTTING — Publication-Quality Style
# =====================================================================

# --- Color palette ---
COLOR_BG        = '#0f1923'      # Deep dark navy background
COLOR_TARGET    = '#56CCF2'      # Cyan-blue for target mRNA backbone
COLOR_SENSOR    = '#F2994A'      # Warm amber for sensor RNA backbone
COLOR_HBOND     = '#ffffff'      # White H-bonds (subtle)
COLOR_HBOND_ALT = '#a0a0a0'      # Alternate grey for thin bonds

# Base colors — soft pastel fills
BASE_COLORS = {
    'A': '#FF6B6B',   # Coral red
    'U': '#4ECDC4',   # Teal
    'G': '#45B7D1',   # Sky blue
    'C': '#FFA07A',   # Light salmon
}

# Edge glow colors per strand
EDGE_TARGET = '#56CCF2'
EDGE_SENSOR = '#F2994A'

# --- Figure setup ---
fig, ax = plt.subplots(figsize=(18, 14), facecolor=COLOR_BG)
ax.set_facecolor(COLOR_BG)

# --- Title ---
title_text = ax.set_title(
    f"dsRNA Dimer · Expanded Structure\n"
    f"Target mRNA Fragment  ×  Sensor RNA   |   MFE: {mfe:.2f} kcal/mol",
    fontsize=18, fontweight='bold', color='#e0e0e0',
    fontfamily='sans-serif', pad=24, linespacing=1.6
)
title_text.set_path_effects([
    pe.withStroke(linewidth=0, foreground=COLOR_BG)
])

# --- 1. Hydrogen Bonds (draw first, lowest layer) ---
for (i, j) in pairs:
    c_i = coord_map[i]
    c_j = coord_map[j]
    ax.plot(
        [xs[c_i], xs[c_j]], [ys[c_i], ys[c_j]],
        color=COLOR_HBOND, linestyle='-', linewidth=0.6,
        alpha=0.25, zorder=1
    )

# --- 2. Backbone glow effect + solid line ---
# Glow layer (thick, blurred)
ax.plot(xs[:target_len], ys[:target_len],
        color=COLOR_TARGET, linewidth=6, alpha=0.12, zorder=2, solid_capstyle='round')
ax.plot(xs[target_len:], ys[target_len:],
        color=COLOR_SENSOR, linewidth=6, alpha=0.12, zorder=2, solid_capstyle='round')

# Core backbone
ax.plot(xs[:target_len], ys[:target_len],
        color=COLOR_TARGET, linewidth=1.8, alpha=0.85, zorder=3,
        solid_capstyle='round', label='Target mRNA')
ax.plot(xs[target_len:], ys[target_len:],
        color=COLOR_SENSOR, linewidth=1.8, alpha=0.85, zorder=3,
        solid_capstyle='round', label='Sensor RNA')

# --- 3. Bases (nucleotide beads) ---
bead_size = 30
font_size = 1.6

for i in range(len(xs)):
    base = seq_no_amp[i]
    fill_color = BASE_COLORS.get(base, '#888888')
    edge_color = EDGE_TARGET if i < target_len else EDGE_SENSOR

    # Outer glow
    ax.scatter(xs[i], ys[i], s=bead_size * 3.5, c=edge_color,
               alpha=0.08, edgecolors='none', zorder=4)

    # Bead
    ax.scatter(xs[i], ys[i], s=bead_size, c=fill_color,
               edgecolors=edge_color, linewidth=0.6, alpha=0.92, zorder=5)

    # Base letter
    ax.text(xs[i], ys[i], base, fontsize=font_size,
            ha='center', va='center', fontweight='bold',
            color='#1a1a2e', zorder=6,
            fontfamily='monospace')

# --- 4. Legend with refined styling ---
legend = ax.legend(
    loc='upper right', fontsize=11, frameon=True,
    facecolor='#1a2332', edgecolor='#2a3a4a',
    labelcolor='#d0d0d0', borderpad=1.0,
    handlelength=2.5, handleheight=1.2
)
legend.get_frame().set_alpha(0.85)

# --- 5. Base-color legend (small annotation) ---
legend_x = 0.02
legend_y = 0.04
base_labels = [('A', BASE_COLORS['A']), ('U', BASE_COLORS['U']),
               ('G', BASE_COLORS['G']), ('C', BASE_COLORS['C'])]
for idx_b, (bl, bc) in enumerate(base_labels):
    ax.annotate(
        f'  {bl}', xy=(legend_x + idx_b * 0.045, legend_y),
        xycoords='axes fraction', fontsize=10, fontweight='bold',
        color=bc, fontfamily='monospace',
        bbox=dict(boxstyle='round,pad=0.25', facecolor=COLOR_BG,
                  edgecolor=bc, alpha=0.6, linewidth=1.2)
    )

# --- Axes cleanup ---
ax.set_aspect('equal')
ax.axis('off')

# --- Subtle border ---
for spine in ax.spines.values():
    spine.set_visible(False)

# --- Save ---
output_file = os.path.join(os.getcwd(), 'True_dsRNA_Dimer_styled.png')
fig.savefig(output_file, dpi=300, bbox_inches='tight',
            facecolor=COLOR_BG, edgecolor='none', pad_inches=0.4)
plt.close(fig)
print(f"[SUCCESS] Styled dimer plotted and saved to: {output_file}")

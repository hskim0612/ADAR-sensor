import RNA
import matplotlib.pyplot as plt
import os

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
expansion_factor = 2.0  # 1.0 means we add 1x the original distance (i.e. double it)

S_x = [0.0] * len(xs)
S_y = [0.0] * len(ys)
is_paired = [False] * len(xs)

# Calculate expansion vectors for paired bases
for (i, j) in pairs:
    c_i = coord_map[i]
    c_j = coord_map[j]
    
    dx = xs[c_i] - xs[c_j]
    dy = ys[c_i] - ys[c_j]
    
    # Push c_i away from c_j
    S_x[c_i] = 0.5 * dx
    S_y[c_i] = 0.5 * dy
    is_paired[c_i] = True
    
    # Push c_j away from c_i
    S_x[c_j] = -0.5 * dx
    S_y[c_j] = -0.5 * dy
    is_paired[c_j] = True

# Linearly interpolate expansion vectors for unpaired bases (Target Strand)
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

# Linearly interpolate expansion vectors for unpaired bases (Sensor Strand)
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

# Apply expansion shift to coordinates
for n in range(len(xs)):
    xs[n] += S_x[n] * expansion_factor
    ys[n] += S_y[n] * expansion_factor

# =====================================================================
# PLOTTING
# =====================================================================
plt.figure(figsize=(16, 12))
plt.title(f"True dsRNA Dimer Visualization (Expanded Distance)\nTarget mRNA Fragment vs Sensor (MFE: {mfe:.2f} kcal/mol)", fontsize=16)

# 1. Draw Hydrogen Bonds (SOLID line to indicate matches clearly)
for (i, j) in pairs:
    c_i = coord_map[i]
    c_j = coord_map[j]
    # Solid line instead of dashed, slightly thicker and distinct color
    plt.plot([xs[c_i], xs[c_j]], [ys[c_i], ys[c_j]], color='#888888', linestyle='-', linewidth=1, zorder=1)

# 2. Draw Backbone
plt.plot(xs[:target_len], ys[:target_len], color='#1f77b4', linewidth=2, zorder=2, label='Target mRNA')
plt.plot(xs[target_len:], ys[target_len:], color='#ff7f0e', linewidth=2, zorder=2, label='Sensor RNA')

# 3. Draw Bases (Beads) with specified 1/10 size and 1/5 font size
for i in range(len(xs)):
    base = seq_no_amp[i]
    color = 'white'
    edgecolor = '#1f77b4' if i < target_len else '#ff7f0e'
    
    if base == 'A': color = '#ffebee'
    elif base == 'U': color = '#e3f2fd'
    elif base == 'G': color = '#e8f5e9'
    elif base == 'C': color = '#fff3e0'

    plt.scatter(xs[i], ys[i], s=20, c=color, edgecolors=edgecolor, linewidth=1, zorder=3)
    plt.text(xs[i], ys[i], base, fontsize=3, ha='center', va='center', fontweight='bold', zorder=4)

plt.axis('equal')
plt.axis('off')
plt.legend(loc='upper right', fontsize=12)

output_file = os.path.join(os.getcwd(), 'True_dsRNA_Dimer.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"[SUCCESS] Expanded dimer plotted and saved to: {output_file}")

import RNA
import matplotlib.pyplot as plt
import os
import sys

def visualize_binding(target_full_seq, sensor_seq, title, output_path):
    target_full = target_full_seq.upper().replace(' ', '').replace('\n', '').replace('T', 'U')
    sensor = sensor_seq.upper().replace(' ', '').replace('\n', '').replace('T', 'U')

    print(f"Analyzing binding for: {title}")
    
    # 1. Optimal duplex search
    duplex = RNA.duplexfold(target_full, sensor)
    if not duplex or duplex.energy >= 0:
        print("Error: Could not find a stable binding site.")
        return False

    struct_t_full, struct_s_full = duplex.structure.split('&')
    t_end = duplex.i
    t_start = duplex.i - len(struct_t_full) + 1
    
    f_start = max(0, t_start - 51)
    f_end = min(len(target_full), t_end + 50)
    target_frag = target_full[f_start:f_end]
    
    # 2. Fragment analysis
    complex_seq = f"{target_frag}&{sensor}"
    struct, mfe = RNA.cofold(complex_seq)
    
    # Use a specific coordinate calculation method if available, 
    # otherwise we manually adjust to prevent circularity.
    coords = RNA.get_xy_coordinates(struct)
    
    xs, ys = [], []
    for i in range(len(complex_seq)):
        if complex_seq[i] != '&':
            try:
                pt = coords.get(len(xs))
                xs.append(pt.X); ys.append(pt.Y)
            except: break

    # 3. Expansion and Manual De-circularization
    # We apply a shift to the sensor strand to push its ends away from the target ends
    expansion_factor = 1.0 
    S_x, S_y = [0.0]*len(xs), [0.0]*len(ys)
    is_paired = [False]*len(xs)
    paired_indices = set()
    
    pairs, stack, coord_map, idx_c = [], [], {}, 0
    for i, char in enumerate(struct):
        if char == '(': stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs.append((j, i))
                paired_indices.add(i); paired_indices.add(j)
        if char != '&':
            coord_map[i] = idx_c; idx_c += 1

    for (i, j) in pairs:
        c_i, c_j = coord_map[i], coord_map[j]
        dx, dy = xs[c_i] - xs[c_j], ys[c_i] - ys[c_j]
        S_x[c_i], S_y[c_i] = 0.5*dx, 0.5*dy
        S_x[c_j], S_y[c_j] = -0.5*dx, -0.5*dy
        is_paired[c_i] = is_paired[c_j] = True

    target_len = len(target_frag)
    def interpolate(start_idx, end_idx):
        for n in range(start_idx, end_idx):
            if not is_paired[n]:
                L = next((s for s in range(n-1, start_idx-1, -1) if is_paired[s]), -1)
                R = next((s for s in range(n+1, end_idx) if is_paired[s]), -1)
                if L != -1 and R != -1:
                    dL, dR = n-L, R-n
                    S_x[n] = (S_x[L]*dR + S_x[R]*dL)/(dL+dR)
                    S_y[n] = (S_y[L]*dR + S_y[R]*dL)/(dL+dR)
                elif L != -1: S_x[n], S_y[n] = S_x[L], S_y[L]
                elif R != -1: S_x[n], S_y[n] = S_x[R], S_y[R]

    interpolate(0, target_len)
    interpolate(target_len, len(xs))

    # Apply expansion
    for n in range(len(xs)):
        xs[n] += S_x[n] * expansion_factor
        ys[n] += S_y[n] * expansion_factor

    # 4. Final Rendering
    plt.figure(figsize=(18, 14))
    plt.title(f"{title}\nTarget Region: mRNA {t_start}-{t_end} | MFE: {mfe:.2f} kcal/mol (Red = Bulge/Mismatch)", fontsize=16)
    
    for (i, j) in pairs:
        c_i, c_j = coord_map[i], coord_map[j]
        plt.plot([xs[c_i], xs[c_j]], [ys[c_i], ys[c_j]], color='#AAAAAA', linestyle='-', linewidth=2, zorder=1)

    # Drawing Backbones with small gaps at ends to emphasize separation
    plt.plot(xs[:target_len], ys[:target_len], color='#1f77b4', alpha=0.4, linewidth=1.5, zorder=2, label='Target mRNA')
    plt.plot(xs[target_len:], ys[target_len:], color='#ff7f0e', alpha=0.4, linewidth=1.5, zorder=2, label='Sensor RNA')

    seq_no_amp = complex_seq.replace('&', '')
    for i in range(len(struct)):
        if i == struct.find('&'): continue
        c_idx = coord_map[i]
        base = seq_no_amp[c_idx]
        
        if i in paired_indices:
            color, edge, fcolor = '#F5F5F5', '#DDDDDD', '#AAAAAA'
        else:
            color, edge, fcolor = '#FF1744', '#D50000', 'white'

        plt.scatter(xs[c_idx], ys[c_idx], s=25, c=color, edgecolors=edge, linewidth=0.5, zorder=3)
        plt.text(xs[c_idx], ys[c_idx], base, fontsize=1.4, ha='center', va='center', color=fcolor, fontweight='bold', zorder=4)

    plt.axis('equal'); plt.axis('off')
    plt.legend(loc='upper right', fontsize=12)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Successfully saved plot: {output_path}")
    return True

if __name__ == "__main__":
    if len(sys.argv) >= 5:
        visualize_binding(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

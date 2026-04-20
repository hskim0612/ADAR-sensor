import pandas as pd
import RNA
import matplotlib.pyplot as plt
import os
import collections

# Python 3.10+ Compatibility Fix for collections.Mapping
if not hasattr(collections, 'Mapping'):
    import collections.abc
    collections.Mapping = collections.abc.Mapping
    collections.MutableMapping = collections.abc.MutableMapping
    collections.Sequence = collections.abc.Sequence

# Try to import forgi
try:
    import forgi.visual.mplotlib as fvm
    import forgi.graph.bulge_graph as fgb
    FORGI_AVAILABLE = True
except ImportError:
    FORGI_AVAILABLE = False
    print("⚠️ Forgi not available. Falling back to custom matplotlib plot.")

def visualize_c3_elite():
    csv_path = "c3_300bp_alubert_v2_refined.csv"
    if not os.path.exists(csv_path):
        print(f"❌ CSV not found: {csv_path}")
        return

    df = pd.read_csv(csv_path)
    elite_sites = [3436, 2457, 4155]
    
    if not os.path.exists("outcome"):
        os.makedirs("outcome")

    for site in elite_sites:
        try:
            # Match CACCA_Start
            target_row = df[df['CACCA_Start'] == site]
            if len(target_row) == 0:
                print(f"⚠️ Site {site} not found in {csv_path}. Skipping.")
                continue
            
            row = target_row.iloc[0]
            target = row['Target_Context'][:300]
            sensor = row['Sensor_300bp_Refined']
            
            print(f"🎨 Visualizing Site {site}...")
            
            # ViennaRNA cofold
            complex_seq = f"{target}&{sensor}"
            structure, energy = RNA.cofold(complex_seq)
            
            title = f"C3 300bp AI-Alu Sensor (Site {site})"
            output_path = f"outcome/C3_300bp_AI_Alu_Site_{site}_Forgi.png"
            
            if FORGI_AVAILABLE:
                # Use Forgi
                forgi_seq = complex_seq.replace('&', 'AAAA')
                forgi_struct = structure.replace('&', '....')
                
                bg = fgb.BulgeGraph()
                bg.from_dotbracket(forgi_struct, forgi_seq)
                
                plt.figure(figsize=(15, 15))
                fvm.plot_rna(bg) # Removed text_kwargs
                plt.title(f"{title}\nMFE: {energy:.2f} kcal/mol", fontsize=15)
                plt.savefig(output_path, dpi=300)
                plt.close()
                print(f"   ✅ Saved Forgi plot: {output_path}")
            else:
                # Fallback to a simpler plot or just print energy
                print(f"   ℹ️ MFE: {energy:.2f} kcal/mol (No Forgi plot generated)")
                
        except Exception as e:
            print(f"❌ Error at Site {site}: {e}")

if __name__ == "__main__":
    visualize_c3_elite()

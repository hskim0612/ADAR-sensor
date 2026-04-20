import sys
import collections
import os

# Essential monkey patches for Python 3.12 compatibility with forgi
if not hasattr(collections, 'Mapping'):
    import collections.abc
    collections.Mapping = collections.abc.Mapping
if not hasattr(collections, 'MutableMapping'):
    import collections.abc
    collections.MutableMapping = collections.abc.MutableMapping
if not hasattr(collections, 'Sequence'):
    import collections.abc
    collections.Sequence = collections.abc.Sequence

import matplotlib.pyplot as plt
import RNA
import forgi.graph.bulge_graph as fgb
import forgi.visual.mplotlib as fvm

# The exact target fragment and sensor we found
target_frag = 'UUCCCAGGCAAAAAACUAGUGCUGUCCAGUGAGAAGACUGUGCUGACCCCUGCCACCAACCACAUGGGCAACGUCACCUUCACGAUCCCAGCCAACAGGGAGUUCAAGUCAGAAAAGGGGCGCAACAAGUUCGUGACCGUGCAGGCCACCUUCGGGACCCAAGUGGUGGAGAAGGUGGUGCUGGUCAGCCUGCAGAGCGGGUA'
sensor = 'AUCGGGUCUGCGGGCUGAGCAGCACCGCCUUCGUCACCGCCUGGGUACGGAGGGAGACCUGCACGGUCACGACCCUGAUGCCCCCCUUUUGUGACGCGAUAGCGCUGUGUGCUGGGACCGAGAAGGCGACGGUGCCGAUGUGGUUGGUCGCUGGCCUCAGCAUCGCCUUCUCCCUGGACCGCUCACCUUUUUUUUCUGGGA'

# Cofold calculation
complex_seq = f"{target_frag}&{sensor}"
structure, energy = RNA.cofold(complex_seq)

# Connect with a loop for visualization (Forgi requirement)
# This simulates the two strands as one long strand folded on itself
full_seq = f"{target_frag}GAAA{sensor}"
full_struct = structure.replace('&', '(...)')

output_path = os.path.join(os.getcwd(), 'dsRNA_Binding_Structure.png')

try:
    bg = fgb.BulgeGraph()
    bg.from_dotbracket(full_struct, full_seq)
    
    plt.figure(figsize=(15, 10))
    fvm.plot_rna(bg, lighten=0.7, text_kwargs={'fontsize': 6})
    plt.title(f"Target mRNA (Fragment) vs Sensor dsRNA
MFE: {energy:.2f} kcal/mol", fontsize=14)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Saved beautiful dsRNA structure plot: {output_path}")
except Exception as e:
    print(f"Error plotting: {e}")

import pandas as pd
from ete4 import Tree

# Loading tree 
t = Tree(open("../dados/16S_final.tree.treefile"))

# Set outgroup
tips = pd.read_table("../dados/tip_meta.tsv")
tipsdict = tips.set_index('label').to_dict(orient='index')

# Loading node annotations
for node in t.traverse():
    if node.is_leaf:
        node.add_props(**tipsdict[node.name])
        node.name = tipsdict[node.name]['species']

# Reroot
exclude = [
    "Rubripirellula lacrimiformis",
    "Paludisphaera sp001610895",
    "Coprothermobacter proteolyticus",
    "JAADFF01 sp013153035",
]
outgroup = tips.query(f'''phylum != "Bacillota" and species not in {exclude}''')['species'].tolist()
t.set_outgroup(t.common_ancestor(outgroup))

t.write("16S_final_itol.nwk", parser=1)

t.explore()

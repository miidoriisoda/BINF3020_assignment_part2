import json
import os

import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO, Entrez, SeqIO
from matplotlib.patches import Patch

from helper_functions import (
    neighbor_joining,
    reroot_tree,
    sort_children_by_leaves,
    plot_nj_tree,
)

Entrez.email = "r.nayager@unsw.edu.au"

def hamming_matrix(seqs):
    """Calculate pairwise Hamming distances between aligned sequences."""
    arr = np.array([list(str(seq).upper()) for seq in seqs])

    n = len(arr)
    D = np.zeros((n, n), dtype=float)

    for i in range(n):
        D[i] = np.sum(arr != arr[i], axis=1)

    return D


def fetch_groups(accessions, cache="groups.json"):
    """Map accession numbers to coronavirus groups using NCBI taxonomy."""

    groups = {}

    # Load any previously cached results.
    if os.path.exists(cache):
        with open(cache) as f:
            groups = json.load(f)

    # Fetch only accessions that are missing from the cache
    missing = [acc for acc in accessions if acc not in groups]

    if missing:
        handle = Entrez.efetch(
            db="nuccore",
            id=",".join(missing),
            rettype="gb",
            retmode="text",
        )

        try:
            for rec in SeqIO.parse(handle, "genbank"):
                accession = rec.id.split(".")[0]
                lineage = rec.annotations.get("taxonomy", [])

                group = next(
                    (
                        taxon
                        for taxon in lineage
                        if taxon in {
                            "Alphacoronavirus",
                            "Betacoronavirus",
                            "Gammacoronavirus",
                            "Deltacoronavirus",
                        }
                    ),
                    "Other",
                )

                groups[accession] = group
        finally:
            handle.close()

        with open(cache, "w") as f:
            json.dump(groups, f, indent=2)

    return groups


def find_node(node, name):
    """Find a node by its name."""
    if node is None:
        return None

    if node.name == name:
        return node

    found = find_node(node.left, name)
    if found is not None:
        return found

    return find_node(node.right, name)


# --------------------------------------------------
# Read the multiple sequence alignment using AlignIO
# --------------------------------------------------
alignment = AlignIO.read(
    "data/clustal_MSA_coronavirus.fasta",
    "fasta",
)

records = list(alignment)

# Currently in format of name|accession
labels = [record.id.split("|")[0].strip() for record in records]

accessions = [
    record.id.split("|")[1].strip().split(".")[0]
    for record in records
]

seqs = [str(record.seq) for record in records]

assert len(set(labels)) == len(labels), "Duplicate sequence names"
assert len(set(accessions)) == len(accessions), "Duplicate accessions"

# Double-check 31 sequences were detected
print(f"Loaded {len(records)} aligned sequences")
print(f"Alignment length: {alignment.get_alignment_length()}")

# --------------------------------------------------
# (a) Calculate distances and construct NJ tree
# --------------------------------------------------

D = hamming_matrix(seqs)
tree = neighbor_joining(D, labels)

# --------------------------------------------------
# (b) Reroot using Breda virus
# --------------------------------------------------

breda = next(
    (label for label in labels if "Breda" in label),
    None,
)

if breda is None:
    raise ValueError("Breda virus was not found in the alignment.")

outgroup = find_node(tree, breda)

if outgroup is None:
    raise ValueError(f"Could not find Breda virus node: {breda}")

if tree.left is outgroup or tree.right is outgroup:
    rooted = tree
else:
    rooted = reroot_tree(tree, outgroup)

# --------------------------------------------------
# (c) Sort children by the number of descendant leaves
# --------------------------------------------------

sort_children_by_leaves(rooted)

# --------------------------------------------------
# (d) Fetch taxonomy and assign group colours
# --------------------------------------------------

groups = fetch_groups(accessions)

palette = {
    "Alphacoronavirus": "#FFCBE1",
    "Betacoronavirus": "#D6E5BD",
    "Gammacoronavirus": "#BCD8EC",
    "Deltacoronavirus": "#DCCCEC",
}

leaf_colours = {
    label: palette.get(groups.get(accession, "Other"), "#F9E1A8")
    for label, accession in zip(labels, accessions)
}

# Plot tree
fig, ax = plt.subplots(figsize=(11, 11))

plot_nj_tree(
    rooted,
    ax=ax,
    leaf_colours=leaf_colours,
    legend=palette,
    show_branch_lengths=False
)

# Add a legend for the four major coronavirus groups.
legend_handles = [
    Patch(facecolor=colour, label=group)
    for group, colour in palette.items()
]

legend_handles.append(
    Patch(facecolor="#F9E1A8", label="Other")
)

ax.legend(handles=legend_handles, loc="best", fontsize=9)

fig.tight_layout()
plt.savefig("problem3.svg", format="svg", bbox_inches="tight")
plt.show()

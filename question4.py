import os
import json
import numpy as np
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.Align import substitution_matrices

Entrez.email = "r.nayager@student.unsw.edu.au"

accession_codes = {
    "Human-SARS-CoV-2": "NC_045512",
    "Bat-CoV RaTG13": "MN996532",
    "Pangolin-CoV MP789": "MT121216",
}

blosum90 = substitution_matrices.load("BLOSUM90")
ALIGNMENT_CACHE = "spike_alignments.json"

def save_alignments(results, filename=ALIGNMENT_CACHE):
    """Save aligned sequences and scores to a JSON file."""
    data = {}

    for name, result in results.items():
        data[name] = {
            "human": result["human"],
            "other": result["other"],
            "score": float(result["score"]),
        }

    with open(filename, "w") as f:
        json.dump(data, f)

def load_alignments(filename=ALIGNMENT_CACHE):
    """Load previously calculated alignments."""
    with open(filename, "r") as f:
        data = json.load(f)

    return data

def blosum90_matrix(x, y, default=-8):
    try:
        result = blosum90[x][y]
        return result
    except KeyError:
        return default
# Global alignment from previous task
def global_alignment(seq1, seq2, scoring_function):
    # Gap penalty
    gap = -1

    n = len(seq1)
    m = len(seq2)

    # ---------------------------------------------------------
    # 1. Create and initialise the dynamic programming matrix
    # ---------------------------------------------------------

    # score[i][j] = best alignment score for
    # seq1[:i] and seq2[:j]
    score = [[0] * (m + 1) for _ in range(n + 1)]

    # Aligning a sequence with an empty sequence requires
    # inserting a gap for every character.
    for i in range(1, n + 1):
        score[i][0] = i * gap

    for j in range(1, m + 1):
        score[0][j] = j * gap

    # ---------------------------------------------------------
    # 2. Fill the dynamic programming matrix
    # ---------------------------------------------------------

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # Option 1: align the two characters
            diagonal = (
                score[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            # Option 2: align seq1[i-1] with a gap
            up = score[i - 1][j] + gap

            # Option 3: align seq2[j-1] with a gap
            left = score[i][j - 1] + gap

            # Keep the best possible score
            score[i][j] = max(diagonal, up, left)

    # ---------------------------------------------------------
    # 3. Traceback
    # ---------------------------------------------------------

    aligned_seq1 = []
    aligned_seq2 = []

    i = n
    j = m

    while i > 0 or j > 0:

        # Diagonal: seq1[i-1] aligned with seq2[j-1]
        if (
            i > 0
            and j > 0
            and score[i][j]
            == score[i - 1][j - 1]
            + scoring_function(seq1[i - 1], seq2[j - 1])
        ):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])

            i -= 1
            j -= 1

        # Up: seq1[i-1] aligned with a gap
        elif (
            i > 0
            and score[i][j] == score[i - 1][j] + gap
        ):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")

            i -= 1

        # Left: seq2[j-1] aligned with a gap
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])

            j -= 1

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, float(score[n][m])

def fetch_record(accession):
    """Retrieve one GenBank record from NCBI."""
    with Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gb",
        retmode="text"
    ) as handle:
        return SeqIO.read(handle, "genbank")

def extract_spike(record):
    """Extract and translate the spike coding sequence."""
    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene = feature.qualifiers.get("gene", [""])[0].lower()
        product = feature.qualifiers.get("product", [""])[0].lower()

        if gene == "s" or "surface glycoprotein" in product or "spike protein" in product:
            nucleotide_sequence = feature.extract(record.seq)
            return str(nucleotide_sequence.translate())

    raise ValueError(f"Spike CDS not found in {record.id}")

def window_similarity(aligned_a, aligned_b, window=100):
    """Fraction of identical, non-gap amino acids per alignment window."""
    if len(aligned_a) != len(aligned_b):
        raise ValueError("Aligned sequences must have equal lengths.")

    if len(aligned_a) < window:
        raise ValueError("Alignment is shorter than the selected window.")

    matches = np.array([
        (a == b) and (a != "-") and (b != "-")
        for a, b in zip(aligned_a, aligned_b)
    ], dtype=float)

    return np.convolve(
        matches,
        np.ones(window) / window,
        mode="valid"
    )

# Retrieve and extract spike sequences
records = {
    name: fetch_record(accession)
    for name, accession in accession_codes.items()
}

spike_sequences = {
    name: extract_spike(record)
    for name, record in records.items()
}

window = 100
results = {}

if os.path.exists(ALIGNMENT_CACHE):
    print("Loading saved alignments...")
    cached = load_alignments()

    for name, result in cached.items():
        similarity = window_similarity(
            result["human"],
            result["other"],
            window
        )

        human_positions = np.cumsum([
            aa != "-" for aa in result["human"]
        ])

        x = human_positions[
            window // 2 : window // 2 + len(similarity)
        ]

        results[name] = {
            **result,
            "similarity": similarity,
            "x": x,
        }

else:
    print("Calculating global alignments...")
    for other in ("Bat-CoV RaTG13", "Pangolin-CoV MP789"):
        aligned_human, aligned_other, score = global_alignment(
            spike_sequences["Human-SARS-CoV-2"],
            spike_sequences[other],
            blosum90_matrix
        )

        results[other] = {
            "human": aligned_human,
            "other": aligned_other,
            "score": score,
        }

    save_alignments(results)

    # Calculate plotting data after saving the expensive alignments
    for name, result in results.items():
        similarity = window_similarity(
            result["human"],
            result["other"],
            window
        )

        human_positions = np.cumsum([
            aa != "-" for aa in result["human"]
        ])

        x = human_positions[
            window // 2 : window // 2 + len(similarity)
        ]

        result["similarity"] = similarity
        result["x"] = x

# Plot spike similarity
colours = {
    "Bat-CoV RaTG13": "#FFCBE1",
    "Pangolin-CoV MP789": "#BCD8EC",
}

fig, ax = plt.subplots(figsize=(14, 9))

for name, result in results.items():
    ax.plot(
        result["x"],
        result["similarity"],
        color=colours[name],
        label=f"SARS-CoV-2 vs {name}"
    )

# SARS-CoV-2 RBD boundaries
ax.axvspan(
    319, 541,
    color="grey",
    alpha=0.2,
    label="RBD"
)

ax.set(
    xlabel="SARS-CoV-2 spike amino-acid position",
    ylabel=f"Similarity",
    ylim=(0, 1.05),
    title="Spike protein similarity"
)
ax.legend(loc="lower left")
fig.tight_layout()
fig.savefig("spike_similarity.svg", bbox_inches="tight")
plt.show()

for name, result in results.items():
    print(f"{name}: alignment score = {result['score']}")
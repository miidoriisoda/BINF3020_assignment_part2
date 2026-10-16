import json
import os

import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO, Entrez, SeqIO
from matplotlib.patches import Patch

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

print(labels)
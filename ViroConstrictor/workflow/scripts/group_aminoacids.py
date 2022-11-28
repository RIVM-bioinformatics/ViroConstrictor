import ast
import sys

import pandas as pd
from Bio import SeqIO

_, input_files, output_files, space = sys.argv

input_files = input_files.split()
output_files = output_files.split()
space = pd.DataFrame.from_dict(ast.literal_eval(space))

## process the input files into individual amino acid records and store the information in the SeqRecords dataframe
SeqRecords = pd.DataFrame()
for f in input_files:
    for s in list(SeqIO.parse(f, "fasta")):
        sample = s.id.split(".")[0]
        aa_feat = s.id.split(".")[1]
        # select the rows that contain the sample name from snakemake.params.space
        data = space.loc[space["sample"] == sample]
        data = data.explode("AA_FEAT_NAMES").reset_index(drop=True)[
            ["sample", "Virus", "RefID", "AA_FEAT_NAMES"]
        ]
        data = data.loc[data["AA_FEAT_NAMES"] == aa_feat]
        data["AA_SEQ"] = str(s.seq)
        SeqRecords = pd.concat([SeqRecords, data], ignore_index=True)


for x in SeqRecords["Virus"].unique():
    for y in SeqRecords.loc[SeqRecords["Virus"] == x]["RefID"].unique():
        for z in SeqRecords.loc[
            (SeqRecords["Virus"] == x) & (SeqRecords["RefID"] == y)
        ]["AA_FEAT_NAMES"].unique():
            if f"results/Virus~{x}/RefID~{y}/aminoacids/{z}.faa" in output_files:
                data = SeqRecords.loc[
                    (SeqRecords["Virus"] == x)
                    & (SeqRecords["RefID"] == y)
                    & (SeqRecords["AA_FEAT_NAMES"] == z)
                ]
                with open(f"results/Virus~{x}/RefID~{y}/aminoacids/{z}.faa", "w") as f:
                    for index, row in data.iterrows():
                        f.write(
                            f">{row['sample']}.{row['AA_FEAT_NAMES']}\n{row['AA_SEQ']}\n"
                        )

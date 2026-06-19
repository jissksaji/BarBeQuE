#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import re
import sys
from pathlib import Path

"""parser= argparse.Argumentparser(description= "Script options")
parser.add_argumet("--taxids",help = "insert the taxids to analyse")
parser.add_argumet("--report ",help = "The consensus report to anylse")
parser.add_argument("--reference",help= "taxonomy lists from db")"""


workdir = Path("/home/saj/barbeque-alpha/BarBeQuE/work/e6/07126346f6859f00a2d38ba675e751")

teeliste = workdir / "teeliste.tsv"
consensus = workdir / "ITS2_collapsed_custom.cluster_consensus.tsv"
taxonomy = workdir / "ITS2_collapsed_custom.cluster_taxonomy.tsv"   

print("Workdir:", workdir)
print("Workdir exists:", workdir.exists())

print("Teeliste path:", teeliste)
print("Teeliste exists:", teeliste.exists())
print("Is file:", teeliste.is_file())

taxids_extracted = []
#selecting the taxids
with open(teeliste, "r", encoding="utf-8") as liste:
    for item in liste:
        columns = item.strip().split("\t")

        if len(columns) < 3:
            continue

        taxid = columns[2].strip()

        if taxid.isdigit():
            taxids_extracted.append(int(taxid))

print(taxids_extracted)

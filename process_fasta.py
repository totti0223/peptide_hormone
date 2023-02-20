import pandas as pd
from utils import *
import json

full_pos = read_fasta(path ="misc/Arabidopdis_peptide_full.fasta")


signaldf = gff3_to_df("misc/output.gff3")
# get the column 0 of signaldf, and split by "__", assign a new column named common_name
signaldf["common_name"] = signaldf[0].str.split("__", expand = True)[0]

l = []
for record in full_pos:
    annot = record.description
    # search for mature sequence positon
    query = annot.split(" ")[0]
    _df = signaldf[signaldf["common_name"].str.contains(query)]
    # iterate column and print value of _df
    if len(_df) != 0:
        start, end, score = _df[3].values[0], _df[4].values[0], _df[5].values[0]
        l.append([annot,query,start, end, score, record.seq])
    else:
        l.append([annot, query, 0, 0, 0, record.seq])

df = pd.DataFrame(l, columns = ["annot", "query", "SP6_start", "SP6_end", "SP6_score", "seq"])
print(df.dtypes)
df.to_csv("annotated_230220.csv")
#print(df.iloc[0,:])
# output.gff3がシグナル配列そのもののgff3ファイル
# region_output.gff3がシグナル配列の詳細を含むgff3ファイル

# column no 3 to 4 がシグナル配列の開始終了位置

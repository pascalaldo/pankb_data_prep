import pandas as pd

def family_summary_table(family, gtdb_meta_path, family_summary_path):
    df = pd.read_csv(gtdb_meta_path, low_memory = False, index_col=0).loc[:,["Family","gc_percentage", "genome_size"]]
    df['source'] = 'ncbi'
    df['gc_content'] = round((df['gc_percentage'])*0.01, 3)
    df['genome_len'] = df['genome_size'].astype(int)
    df = df.loc[df["Family"] == family, :]

    df.to_csv(family_summary_path)

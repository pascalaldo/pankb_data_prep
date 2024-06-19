import pandas as pd
import json

def full_summary_table(gtdb_meta_path, full_summary_path):
    df = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0).loc[
        :, ["Family", "gc_percentage", "genome_size"]
    ]
    df["source"] = "ncbi"
    df["gc_content"] = round((df["gc_percentage"]) * 0.01, 3)
    df["genome_len"] = df["genome_size"].astype(int)

    df.to_csv(full_summary_path)

def family_summary_table(family, gtdb_meta_path, family_summary_path):
    if not family.startswith("f__"):
        family = f"f__{family}"
    df = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0).loc[
        :, ["Family", "gc_percentage", "genome_size"]
    ]
    df["source"] = "ncbi"
    df["gc_content"] = round((df["gc_percentage"]) * 0.01, 3)
    df["genome_len"] = df["genome_size"].astype(int)
    df = df.loc[df["Family"] == family, :]

    df.to_csv(family_summary_path)


def species_pangenome_summary(
    analysis_name,
    gp_binary_path,
    summary_v2_path,
    gtdb_meta_path,
    species_summary_csv_path,
    species_summary_json_path,
):
    df_gp_binary = pd.read_csv(gp_binary_path, index_col="Gene", low_memory=False)
    df_pangene_summary = pd.read_csv(summary_v2_path, low_memory=False, index_col=0)
    df_gtdb_meta = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0)

    genomes = df_gp_binary.columns.tolist()
    n_genomes = len(genomes)
    n_genes = int(df_pangene_summary.shape[0])

    # TODO
    # print(set(df_gtdb_meta.loc[genomes, "Family"]))
    # print(set(df_gtdb_meta.loc[genomes, "Species"]))
    # assert len(set(df_gtdb_meta.loc[genomes, "Family"])) == 1 # Make sure that all the genomes in the analysis are from the same family
    # assert len(set(df_gtdb_meta.loc[genomes, "Species"])) == 1 # Make sure that all the genomes in the analysis are from the same species

    species = str(df_gtdb_meta.loc[genomes[0], "Organism"]).replace("s__", "")
    family = str(df_gtdb_meta.loc[genomes[0], "Family"]).replace("f__", "")

    core_len = int(
        df_pangene_summary.loc[
            df_pangene_summary["pangenome_class_2"] == "Core", "pangenome_class_2"
        ].count()
    )
    rare_len = int(
        df_pangene_summary.loc[
            df_pangene_summary["pangenome_class_2"] == "Rare", "pangenome_class_2"
        ].count()
    )
    accessory_len = int(
        df_pangene_summary.loc[
            df_pangene_summary["pangenome_class_2"] == "Accessory", "pangenome_class_2"
        ].count()
    )

    df = pd.DataFrame(
        {
            "Pangenome_analyses": [analysis_name],
            "Family": [family],
            "Species": [species],
            "N_of_genome": [n_genomes],
            "N_of_gene": [n_genes],
            "N_of_core": [core_len],
            "N_of_rare": [rare_len],
            "N_of_accessory": [accessory_len],
            "Openness": ["Open"],
        }
    )
    df.to_csv(species_summary_csv_path, index=False)

    json_data = {
        "Family": family,
        "Species": species,
        "Number_of_genome": n_genomes,
        "Gene_class": [core_len, accessory_len, rare_len],
        "Openness": "Open",
    }
    with open(species_summary_json_path, "w") as f:
        json.dump(json_data, f)


def pangenome_summary(species_summary_path, species_list_out_path, genome_count_out_path, gene_count_out_path):
    df = pd.read_csv(species_summary_path)

    json_data = {
        "Family": df["Family"].tolist(),
        "Species": df["Species"].tolist(),
        "Openness": df["Openness"].tolist(),
        "N_of_genome": df["N_of_genome"].astype(int).tolist(),
        "Pangenome_analyses": df["Pangenome_analyses"].tolist(),
        "Gene_class": list(
            map(
                list,
                zip(
                    *[
                        df[x].astype(int).tolist()
                        for x in ["N_of_core", "N_of_accessory", "N_of_rare"]
                    ]
                ),
            )
        ),
    }
    with open(species_list_out_path, "w") as f:
        json.dump(json_data, f)

    family_group = df.groupby("Family")
    family_group["N_of_genome"].sum().sort_values(ascending=False).to_json(genome_count_out_path, orient='index')

    family_group["N_of_gene"].sum().sort_values(ascending=False).to_json(gene_count_out_path, orient='index')

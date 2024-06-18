import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.preprocessing import MinMaxScaler
import numpy as np
import string

def kMeansRes(scaled_data, k, alpha_k=0.02):
    '''
    # Calculating clusters from https://medium.com/towards-data-science/an-approach-for-choosing-number-of-clusters-for-k-means-c28e614ecb2c
    Parameters 
    ----------
    scaled_data: matrix 
        scaled data. rows are samples and columns are features for clustering
    k: int
        current k for applying KMeans
    alpha_k: float
        manually tuned factor that gives penalty to the number of clusters
    Returns 
    -------
    scaled_inertia: float
        scaled inertia value for current k           
    '''
    
    inertia_o = np.square((scaled_data - scaled_data.mean(axis=0))).sum()
    # fit k-means
    kmeans = KMeans(n_clusters=k, random_state=0).fit(scaled_data)
    scaled_inertia = kmeans.inertia_ / inertia_o + alpha_k * k
    return scaled_inertia

def chooseBestKforKMeans(scaled_data, k_range):
    ans = []
    for k in k_range:
        scaled_inertia = kMeansRes(scaled_data, k)
        ans.append((k, scaled_inertia))
    results = pd.DataFrame(ans, columns = ['k','Scaled Inertia']).set_index('k')
    best_k = results.idxmin()[0]
    return best_k, results
        
def generate_mash_list(genomes, gtdb_meta_path, mash_file_path, mash_list_path):
    df_mash = pd.read_csv(mash_file_path, index_col=0)
    df_mash = df_mash.loc[genomes, :]
    df_gtdb = pd.read_csv(gtdb_meta_path, index_col='genome_id')
    df_gtdb = df_gtdb.loc[genomes, :]

    df_mash_corr = df_mash.corr()

    # choose features
    data_for_clustering = df_mash.copy()
    data_for_clustering.fillna(0, inplace=True)

    # create data matrix
    data_matrix = np.matrix(data_for_clustering).astype(float)
    # scale the data
    mms = MinMaxScaler()
    scaled_data = mms.fit_transform(np.asarray(data_matrix))

    # choose k range
    if len(df_mash) <= 21:
        max_range = len(df_mash) - 1
    else:
        max_range = 20

    k_range=range(2, max_range)
    # compute adjusted intertia
    best_k, results = chooseBestKforKMeans(scaled_data, k_range)
    n_clusters = best_k

    # create output folder
    Agg_hc = AgglomerativeClustering(n_clusters = n_clusters, metric = 'euclidean', linkage = 'ward')
    y_hc = Agg_hc.fit_predict(df_mash_corr)
    df_hclusts = pd.DataFrame(index=df_mash_corr.index, columns=['hcluster', 'color_code'])
    df_hclusts['hcluster'] = y_hc

    group_dict = {}

    for i in range(24):
        group_dict[i] = string.ascii_uppercase[i]

    df_hclusts['cluster'] = df_hclusts['hcluster'].apply(lambda c: group_dict[int(c)])
    df_hclusts['genome_id'] = df_hclusts.index
    mash_df = df_hclusts.sort_values(by = 'cluster')[['genome_id', 'cluster']]

    mash_df.to_csv(mash_list_path, index=False)
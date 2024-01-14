import pandas as pd

id_pathway_df = pd.read_csv('/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/ID_PATHWAY_TABLE.txt', sep='\t', header=None, names=['ID', 'PATHWAY'])
ids_coincidentes = set(map(str.strip, open('/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/IDs_coincidentes2.txt')))

data = {'all_pathways': [], 'differentially_expresed': [], 'PATHWAY': []}

for _, row in id_pathway_df.iterrows():
    ids_en_pathway = set(map(str.strip, row['ID'].split(',')))
    
    data['all_pathways'].append(len(ids_en_pathway))
    
    coincidentes = [gene for gene in ids_en_pathway if gene in ids_coincidentes]
    data['differentially_expresed'].append(len(coincidentes))
    
    data['PATHWAY'].append(row['PATHWAY'])

result_df = pd.DataFrame(data)

result_df.to_csv('/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/resultados_wildtype_pks14.csv', sep=',', index=False)


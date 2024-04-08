import pandas as pd
import decoupler as dc

# Gets the three Prior Knowledge sources from the OmniPath database via decoupleR and outputs them as tsv files. 
progeny = dc.get_progeny(organism='human')
# dorothea = dc.get_dorothea(organism='human')
msigdb = dc.get_resource('MSigDB', organism='human')
collectri = dc.get_collectri(organism='human', split_complexes=False)

msigdb_hallmarks = msigdb.loc[msigdb['collection'] == 'hallmark'].rename(columns={'genesymbol':'target', 'geneset':'source'}).drop_duplicates()

progeny.to_csv('progeny__source.tsv', sep='\t', index=False)
# dorothea.to_csv('dorothea__source.tsv', sep='\t', index=False)
collectri.to_csv('collectri__source.tsv', sep='\t', index=False)
msigdb_hallmarks.to_csv('msigdb_hallmarks__source.tsv', sep='\t', index=False)


# for benchmark
# markers = dc.get_resource('PanglaoDB')
# markers = markers[(markers['human']=='True')&(markers['canonical_marker']=='True')]
# markers = markers[~markers.duplicated(['cell_type', 'genesymbol'])].rename(columns={'genesymbol':'target', 'cell_type':'source'})
# markers.to_csv('biomarkers__source.tsv', sep='\t', index=False)

import pandas as pd
import decoupler as dc

# Gets the three Prior Knowledge sources from the OmniPath database via decoupleR and outputs them as tsv files. 
progeny = dc.get_progeny(organism='human')
dorothea = dc.get_dorothea(organism='human')
msigdb = dc.get_resource('MSigDB')

msigdb_hallmarks = msigdb.loc[msigdb['collection'] == 'hallmark'].rename(columns={'genesymbol':'target', 'geneset':'source'}).drop_duplicates()
progeny.to_csv('progeny__source.tsv', sep='\t', index=False)
dorothea.to_csv('dorothea__source.tsv', sep='\t', index=False)
msigdb_hallmarks.to_csv('msigdb_hallmarks__source.tsv', sep='\t', index=False)




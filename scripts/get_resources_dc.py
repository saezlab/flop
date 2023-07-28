import pandas as pd
import decoupler as dc

progeny = dc.get_progeny(organism='human')
dorothea = dc.get_dorothea(organism='human')
msigdb = dc.get_resource('MSigDB')
# random = dorothea
# random['target'] = dorothea['target'].sample(frac=1).values
# random['weight'] = dorothea['weight'].sample(frac=1).values
# random.drop_duplicates(inplace=True, subset=['source', 'target'])
msigdb_hallmarks = msigdb.loc[msigdb['collection'] == 'hallmark'].rename(columns={'genesymbol':'target', 'geneset':'source'}).drop_duplicates()
progeny.to_csv('progeny__source.tsv', sep='\t', index=False)
dorothea.to_csv('dorothea__source.tsv', sep='\t', index=False)
msigdb_hallmarks.to_csv('msigdb_hallmarks__source.tsv', sep='\t', index=False)
# random.to_csv('random__source.tsv', sep='\t', index=False)



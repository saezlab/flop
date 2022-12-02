import pandas as pd
import decoupler as dc
import os


progeny = dc.get_progeny(organism='human')
dorothea = dc.get_dorothea(organism='human')
msigdb = dc.get_resource('MSigDB')
msigdb_hallmarks = msigdb.loc[msigdb['collection'] == 'hallmark'].rename(columns={'genesymbol':'target', 'geneset':'source'}).drop_duplicates()

try:
    os.mkdir('./dc_resources')
except FileExistsError:
    pass

progeny.to_csv('./dc_resources/progeny__source.tsv', sep='\t', index=False)
dorothea.to_csv('./dc_resources/dorothea__source.tsv', sep='\t', index=False)
msigdb_hallmarks.to_csv('./dc_resources/msigdb_hallmarks__source.tsv', sep='\t', index=False)
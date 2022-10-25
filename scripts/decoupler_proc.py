import decoupler as dc
import pandas as pd
import sys
import numpy as np

input_file = sys.argv[1]
input_data = pd.read_csv(input_file, sep='\t').transpose()
input_data.fillna(0, inplace=True)
input_data.rename(columns = input_data.iloc[0], inplace=True)
input_data.drop(input_data.index[0], inplace=True)
mat = input_data.astype(float)
network = dc.get_progeny(organism = 'human')
res_decoupler = dc.run_consensus(mat, network, min_n=0)
#output_file = input_file.replace('__decouplerinput.tsv', '__decoupleroutput.tsv')
print(res_decoupler[0])
#res_decoupler.to_csv(output_file, sep='\t')









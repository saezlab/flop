import decoupler as dc
import pandas as pd
import sys

#Input management
input_file = sys.argv[1]
input_data = pd.read_csv(input_file, sep='\t').transpose()
input_data.fillna(0, inplace=True)
input_data.rename(columns = input_data.iloc[0], inplace=True)
input_data.drop(input_data.index[0], inplace=True)
mat = input_data.astype(float)

#Decoupler analysis
network = dc.get_progeny(organism = 'human')
dc_result = dc.decouple(mat, network, min_n=0)

#Output formatting
methods = ['mlm_estimate', 'ulm_estimate', 'wsum_norm', 'consensus_estimate']
suffixes = ['mlm', 'ulm', 'wsum', 'cons']
selected_results = list(map(dc_result.get, methods))
for i in range(len(methods)):
    newtags = [name + '__' + suffixes[i] for name in selected_results[i].index.values]
    selected_results[i].index = newtags
    dc_output = selected_results[i]
    output_file = input_file.replace('__decouplerinput.tsv', '__{}__decoupleroutput.tsv'.format(suffixes[i]))
    dc_output.to_csv(output_file, sep='\t')










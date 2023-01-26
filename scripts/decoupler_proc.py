import decoupler as dc
import pandas as pd
import sys

#Input management
input_file = sys.argv[1]
resource = sys.argv[2]
scriptdir = sys.argv[0].replace('decoupler_proc.py', '')
resource_file = '{}dc_resources/{}__source.tsv'.format(scriptdir,resource)
input_data = pd.read_csv(input_file, sep='\t').transpose()
input_data.fillna(0, inplace=True)
input_data.rename(columns = input_data.iloc[0], inplace=True)
input_data.drop(input_data.index[0], inplace=True)
mat = input_data.astype(float)

#Decoupler analysis
network = pd.read_csv(resource_file, sep='\t')
if 'weight' in network.columns:
    dc_result = dc.decouple(mat, network)
else:
    dc_result = dc.decouple(mat, network, weight=None)

#Output formatting
methods = ['consensus_estimate']
suffixes = ['cons']
selected_results = list(map(dc_result.get, methods))
for i in range(len(methods)):
    newtags = [name + '__' + suffixes[i] for name in selected_results[i].index.values]
    selected_results[i].index = newtags
    dc_output = selected_results[i]
    output_file = input_file.replace('__decouplerinput.tsv', '__{}__{}__decoupleroutput.tsv'.format(suffixes[i], resource))
    dc_output.to_csv(output_file, sep='\t')








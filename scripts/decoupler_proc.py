import decoupler as dc
import pandas as pd
import sys

#Input management
input_file = sys.argv[1]
resource = sys.argv[2]
scriptdir = sys.argv[0].replace('decoupler_proc.py', '')
resource_file = '{}dc_resources/{}__source.tsv'.format(scriptdir, resource)
input_data = pd.read_csv(input_file, sep='\t').transpose()
input_data.fillna(0, inplace=True)
input_data.rename(columns = input_data.iloc[0], inplace=True)
input_data.drop(input_data.index[0], inplace=True)
mat = input_data.astype(float)

# input_file = 'Sweet18__yes_v_no__stat__decouplerinput.tsv'
# resource = 'dorothea'
# scriptdir = './scripts/'


#Decoupler analysis
network = pd.read_csv(resource_file, sep='\t')
if 'weight' in network.columns:
    dc_result = dc.run_ulm(mat, network)
else:
    dc_result = dc.run_ulm(mat, network, weight=None)

#Output formatting
# methods = ['consensus_estimate']
# suffixes = ['cons']
# selected_results = list(map(dc_result.get, methods))
# for i in range(len(methods)):
#     newtags = [name + '__' + suffixes[i] for name in selected_results[i].index.values]
#     selected_results[i].index = newtags
#     dc_output = selected_results[i]
#     output_file = input_file.replace('__decouplerinput.tsv', '__{}__{}__decoupleroutput.tsv'.format(suffixes[i], resource))
#     dc_output.to_csv(output_file, sep='\t')

selected_results = dc_result[1]
newtags = [name + '__' + 'mean' for name in selected_results.index.values]
selected_results.index = newtags
dc_output = selected_results
output_file = input_file.replace('__decouplerinput.tsv', '__{}__{}__decoupleroutput.tsv'.format('mean', resource))
dc_output.to_csv(output_file, sep='\t')
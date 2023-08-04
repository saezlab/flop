import decoupler as dc
import pandas as pd
import sys

# Performs GSEA with decoupleR. The input is a wide-format table with pipelines as columns and genes as rows.
# The output is a table with the activity scores and p-values for each pipeline.

# Input management
input_file = sys.argv[1]
resource = sys.argv[2]
scriptdir = sys.argv[0].replace('decoupler_proc.py', '')
resource_file = '{}dc_resources/{}__source.tsv'.format(scriptdir, resource)
input_data = pd.read_csv(input_file, sep='\t').transpose()
input_data.fillna(0, inplace=True)
input_data.rename(columns = input_data.iloc[0], inplace=True)
input_data.drop(input_data.index[0], inplace=True)
mat = input_data.astype(float)

# Decoupler analysis
network = pd.read_csv(resource_file, sep='\t')
if 'weight' in network.columns:
    dc_result = dc.run_ulm(mat, network)
else:
    dc_result = dc.run_ulm(mat, network, weight=None)

# Output formatting
# Retrieves the activity scores and p-values from the decoupler result object

act_scores, pvalues = dc_result

indexes = pvalues.index.values
adj_pvalues = pd.DataFrame(index=indexes, columns=pvalues.columns.values)

for i in indexes:
    adj_pvalues.loc[i] = dc.p_adjust_fdr(pvalues.loc[i].values)

act_scores_newtags = [name + '__' + 'ulm' + '__' + 'act' for name in act_scores.index.values]
pval_newtags = [name + '__' + 'ulm' + '__' + 'padj' for name in adj_pvalues.index.values]
act_scores.index = act_scores_newtags
adj_pvalues.index = pval_newtags
dc_output = pd.concat([act_scores, adj_pvalues], axis=0)

output_file = input_file.replace('__decouplerinput.tsv', '__{}__{}__decoupleroutput.tsv'.format('mean', resource))
dc_output.to_csv(output_file, sep='\t')

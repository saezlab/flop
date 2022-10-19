#!/usr/bin/env nextflow

// Define input parameters
params.counts = "$baseDir/data/test_input/GSE103001_GeneLevel_Raw_data.tsv"
params.scripts_dir = "$baseDir/scripts"

// Normalize counts files using VSN
process normalize_vsn {
 
    input:
    path 'raw_counts.tsv'
    path 'scripts_dir'
 
    output:
    path '*__norm.tsv'
    
    script:
    """
    Rscript scripts_dir/vsn.R --counts raw_counts.tsv
    """

}

// Normalize counts files using TMMs
process normalize {
 
    input:
    path 'raw_counts.tsv'
    path 'scripts_dir'
 
    output:
    path '*__norm.tsv'
    
    script:
    """
    Rscript scripts_dir/tmm.R --counts raw_counts.tsv
    """

}

// Normalize counts files using log2 and quantile normalization
process normalize_log2quant {
 
    input:
    path 'raw_counts.tsv'
    path 'scripts_dir'
 
    output:
    path '*__norm.tsv'
    
    script:
    """
    Rscript scripts_dir/log2quant.R --counts raw_counts.tsv
    """

}

// Diff exp with DESeq2
process normalize_log2quant {
 
    input:
    path 'raw_counts.tsv'
    path 'scripts_dir'
 
    output:
    path '*__norm.tsv'
    
    script:
    """
    Rscript scripts_dir/log2quant.R --counts raw_counts.tsv
    """

}


workflow {

    // Apply normalization
    normalize_vsn(params.counts, params.scripts_dir)
    normalize_tmm(params.counts, params.scripts_dir)
    normalize_log2quant(params.counts, params.scripts_dir)

    // Apply diff expression analysis


}
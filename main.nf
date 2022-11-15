#!/usr/bin/env nextflow

// Normalize counts files using VSN
process normalize_vsn {
 
    input:
    path 'scripts_dir'
    path 'raw_counts.tsv'
 
    output:
    path '*__norm.tsv'
    
    script:
    """
    Rscript scripts_dir/vsn.R --counts raw_counts.tsv
    """

}

// Normalize counts files using TMMs
process normalize_tmm {
 
    input:
    path 'scripts_dir'
    path 'raw_counts.tsv'
 
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
    path 'scripts_dir'
    path 'raw_counts.tsv'
 
    output:
    path '*__norm.tsv'
    
    script:
    """
    Rscript scripts_dir/log2quant.R --counts raw_counts.tsv
    """

}

// Diff exp with DESeq2
process diffexp_deseq2 {
 
    input:
    path 'scripts_dir'
    path 'raw_counts.tsv'
    path 'metadata.tsv'
 
    output:
    path '*__deseq2__de.tsv'
    
    script:
    """
    Rscript scripts_dir/deseq2.R --counts raw_counts.tsv --meta metadata.tsv
    """

}


// Diff exp with limma
process diffexp_limma {
 
    input:
    path 'scripts_dir'
    path 'norm_expr.tsv'
    path 'metadata.tsv'
 
    output:
    path '*__limma__de.tsv'
    
    script:
    """
    Rscript scripts_dir/limma.R --norm norm_expr.tsv --meta metadata.tsv
    """

}

// Input parameters
params.counts = "$baseDir/data/test_input/GSE103001_GeneLevel_Raw_data.tsv"
params.metadata = "$baseDir/data/test_input/GSE103001_filtered_metadata.tsv"
params.scripts_dir = "$baseDir/scripts"

// Workflow definition
workflow {

    // Apply normalization
    normalize_vsn(params.scripts_dir, params.counts)
    normalize_tmm(params.scripts_dir, params.counts)
    normalize_log2quant(params.scripts_dir, params.counts)

    // Apply diff expression analysis from raw counts
    diffexp_deseq2(params.scripts_dir, params.counts, params.metadata)

    // Execute limma downstream of VSN normalization
    diffexp_limma(params.scripts_dir, normalize_vsn.out, params.metadata)



}
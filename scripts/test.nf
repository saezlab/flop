params.scripts_dir = projectDir

process normalize {
    publishDir "$params.scripts_dir/results/$datasetID/norm_output", mode: 'copy'
 
    input:
    path scripts_dir
    tuple val(datasetID), val(dataID), path(counts), path(meta)
    each norm_method
 
    output:
    tuple val(datasetID),val(dataID), path ('*__norm.tsv')
    
    script:
    """
    Rscript ${scripts_dir}/${norm_method}.R --counts ${counts}
    """
}

//Database selection
workflow {    
    Channel
        .fromFilePairs("$params.scripts_dir/data/*/*_{*_countdata,*_metadata}.tsv")
        .map{it -> tuple it[1][0].parent.baseName, it[0], it[1][0], it[1][1]}
        .view()
        .set{datasets}

    Channel
        .of('vsn', 'log2quant', 'tmm')
        .set{norm_methods}
    
    //normalization channels
    normalize(params.scripts_dir,datasets, norm_methods)
        .combine(datasets, by: [0,1])
        .view{"Normalization: $it"}
        .set{normalised_vals}

}
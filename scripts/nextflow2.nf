
params.scripts_dir = projectDir

//normalization
process normalize {
    publishDir "$params.scripts_dir/results/$datasetID/norm_output", mode: 'copy'
 
    input:
    path scripts_dir
    tuple val(datasetID), val(dataID), path(counts), path(meta)
    each norm_method
 
    output:
    tuple val(datasetID), val(dataID), path ('*__norm.tsv')
    
    script:
    """
    Rscript ${scripts_dir}/${norm_method}.R --counts ${counts}
    """
}

// Diff exp with limma
process diffexp_limma {
 
    input:
    path scripts_dir
    tuple val(datasetID), val(dataID), path(norm_file), path(counts), path(meta)

    output:
    tuple val(datasetID), val(dataID), path ('*__limma__de.tsv')
    
    script:

    """
    Rscript ${scripts_dir}/limma.R --norm ${norm_file} --meta ${meta}
    """

}

// Diff exp with DESeq2
process diffexp_deseq2 {
 
    input:
    path scripts_dir
    tuple val(datasetID), val(dataID), path(counts), path(meta)
 
    output:
    tuple val(datasetID), val(dataID), path ('*__deseq2__de.tsv')
    
    script:

    """
    Rscript ${scripts_dir}/deseq2.R --counts ${counts} --meta ${meta}
    """

}

// Diff exp with EdgeR
process diffexp_edger{
 
    input:
    path scripts_dir
    tuple val(datasetID), val(dataID), path(counts), path(meta)
 
    output:
    tuple val(datasetID), val(dataID), path ('*__edger__de.tsv')
    
    script:

    """
    Rscript ${scripts_dir}/edger.R --counts ${counts} --meta ${meta}
    """

}

//Differential expression output file merger
process merge_de{
    publishDir "$params.scripts_dir/results/$datasetID/dc_input", mode: 'copy'
 
    input:
    path scripts_dir
    tuple val(datasetID), val(dataID), path(diffexpr_files)
    each method
 
    output:
    tuple val(datasetID), path ('*__decouplerinput.tsv')
    
    script:

    """
    Rscript ${scripts_dir}/merge_de.R --dataset ${dataID} --param ${method}
    """

}

//Functional analysis 
process func_decoupler{
    publishDir "$params.scripts_dir/results/$datasetID/dc_output", mode: 'move'
 
    input:
    path scripts_dir
    tuple val(datasetID), path(decoupler_files)
    each resources
 
    output:
    path ('*__decoupleroutput.tsv')
    
    script:

    """
    python3 ${scripts_dir}/decoupler_proc.py ${decoupler_files} ${resources}
    """
}



workflow {
    Channel
        .fromFilePairs("$params.scripts_dir/data/*/*_{*_countdata,*_metadata}.tsv")
        .map{it -> tuple it[1][0].parent.baseName, it[0], it[1][0], it[1][1]}
        .set {datasets}
    
    Channel
        .of('logFC', 'stat')
        .set {diffexpr_methods}
    
    Channel
        .of('vsn', 'log2quant', 'tmm')
        .set{norm_methods}

    Channel
        .of('progeny', 'dorothea', 'msigdb_hallmarks')
        .set{resources}
    
    //normalization channels
    normalize(params.scripts_dir,datasets, norm_methods)
        .combine(datasets, by: [0,1])
        //.view{"Normalization: $it \n"}
        .set{normalised_vals}

    //differential analysis channels
    diffexp_deseq2(params.scripts_dir,datasets).set {deseq2}
    diffexp_edger(params.scripts_dir,datasets).set {edger}
    
    diffexp_limma(params.scripts_dir,normalised_vals)
        .set{limma}

    //functional analysis channels
    limma
        .mix(deseq2, edger)
        .groupTuple(by:[0,1])
        //.view{"Differential analysis: $it \n"}
        .set {diffexpr_files}
    
    merge_de(params.scripts_dir,diffexpr_files,diffexpr_methods)
        //.view{"Decoupler input: $it \n"}
        .set {mergede}
    
    func_decoupler(params.scripts_dir, mergede, resources)
        .set {decoupler}
}
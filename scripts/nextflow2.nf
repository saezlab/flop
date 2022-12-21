
params.scripts_dir = projectDir

process get_prsources{
    publishDir "$params.scripts_dir/dc_resources", mode: 'copy'

    input:
    path scripts_dir

    output:
    path ("*__source.tsv")

    script:

    """
    python3 ${scripts_dir}/get_resources_dc.py
    """
}

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
    publishDir "$params.scripts_dir/results/$datasetID/dc_output", mode: 'copy'
 
    input:
    path scripts_dir
    tuple val(datasetID), path(decoupler_files)
    each resources
 
    output:
    tuple val(datasetID), path ('*__decoupleroutput.tsv')
    
    script:

    """
    python3 ${scripts_dir}/decoupler_proc.py ${decoupler_files} ${resources}
    """
}


process decoupler_merger{
    publishDir "$params.scripts_dir/results/$datasetID", mode: 'copy'

    input:
    path scripts_dir
    tuple val(datasetID), path (decoupler_results)

    output:
    tuple val(datasetID), path ("*__result.tsv")

    script:

    """
    Rscript ${scripts_dir}/decoupler_merger.R --dataset ${datasetID} --files "${decoupler_results}"
    """
}

process rank_analysis{
    publishDir "$params.scripts_dir/results/$datasetID/plots", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), path (analysis_results)

    output:
    path ("*.png")

    script:

    """
    Rscript ${scripts_dir}/rank_analysis.R --file ${analysis_results}
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
  
    //Prior knowledge sources
    get_prsources(params.scripts_dir)
        .flatten()
        .map{it -> it.baseName.toString().replaceAll(/__source/, "")}
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
        .groupTuple(by:0)
        .map{it -> tuple(it[0],it[1].flatten())}
        //.view{"Decoupler out: $it"}
        .set {decoupler}

    decoupler_merger(params.scripts_dir, decoupler)
        .set {merged_results}

    rank_analysis(params.scripts_dir, merged_results)
        .set {rank}

}
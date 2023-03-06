
params.scripts_dir = projectDir

//Downloads and stores prior knowledge sources
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

//Performs filtering, normalisation and differential expression analysis
process diffexp_analysis{
    input:
    path scripts_dir
    tuple val(datasetID), val(biocontext), path(counts), path(meta)
    each pipelines
    each status
 
    output:
    tuple val(datasetID), val(biocontext), val(status), path("*__de.rds")
    
    script:

    """
    Rscript ${scripts_dir}/diffexp_analysis.R --dataset ${datasetID} --counts "${counts}" --meta "${meta}" --pipeline "${pipelines}" --status ${status} --bio ${biocontext}
    """
}

//Merges the results from different pipelines
process merge_de{
    publishDir "$params.scripts_dir/results/$datasetID/$status/", mode: 'copy'
    input:
    path scripts_dir
    tuple val(datasetID), val(biocontext), val(status), path(diffexpr_files)
    each method
 
    output:
    tuple val(datasetID), val(biocontext), val(status), path ('*__decouplerinput.tsv')

    //afterScript "rm -f ${diffexpr_files}"
    
    script:

    """
    Rscript ${scripts_dir}/merge_de.R --dataset ${datasetID} --bio ${biocontext} --param ${method} --files "${diffexpr_files}"
    """
}

//Functional analysis 
process func_decoupler{
    publishDir "$params.scripts_dir/results/$datasetID/$status/", mode: 'copy'

    input:
    path scripts_dir
    tuple val(datasetID), val(biocontext), val(status), path(decoupler_files)
    each resources
 
    output:
    tuple val(datasetID), val(status), path ('*__decoupleroutput.tsv')

    //afterScript "rm -f ${decoupler_files}"
    
    script:

    """
    python3 ${scripts_dir}/decoupler_proc.py ${decoupler_files} ${resources}
    """
}

//Merges the results from DecoupleR
process decoupler_merger{
    publishDir "$params.scripts_dir/results/$datasetID/$status", mode: 'copy'

    input:
    path scripts_dir
    tuple val(datasetID), val(status), path (decoupler_results)

    output:
    tuple val(datasetID), val(status), path ("*__result.tsv")

    script:

    """
    Rscript ${scripts_dir}/decoupler_merger.R --dataset ${datasetID} --files "${decoupler_results}"
    """
}

//Rank analysis
process rank_analysis{
    publishDir "$params.scripts_dir/results/$datasetID/$status/plots", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), val(status), path (analysis_results)

    output:
    path ("*.png")

    script:

    """
    Rscript ${scripts_dir}/rank_analysis.R --file ${analysis_results}
    """
}

//Rand index analysis
process rand_index_analysis{
    publishDir "$params.scripts_dir/results/$datasetID/$status/plots", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), val(status), path (analysis_results)

    output:
    path ("*.png")

    script:

    """
    Rscript ${scripts_dir}/rand_index_analysis.R --file ${analysis_results}
    """
}


workflow {
    Channel
        .fromFilePairs("$params.scripts_dir/data/*/*_{*_countdata,*_metadata}.tsv")
        .map{it -> tuple it[1][0].parent.baseName, it[0], it[1][0], it[1][1]}
        //.view()
        .set {datasets}
    
    Channel
        .of('vsn_norm limma_analysis', 'tmm_norm limma_analysis', 'log2quant_norm limma_analysis', 'edger_analysis', 'deseq2_analysis')
        .set {pipelines}
    
    Channel
        .of('filtered', 'unfiltered')
        .set {status}
  
    Channel
        .of('logFC', 'stat')
        .set {diffexpr_methods}
    
    get_prsources(params.scripts_dir)
        .flatten()
        .map{it -> it.baseName.toString().replaceAll(/__source/, "")}
        .set{resources}
    
    diffexp_analysis(params.scripts_dir, datasets, pipelines, status)
        .groupTuple(by:[0,1,2])
        //.view{"Differential analysis: $it \n"}
        .set {diffexpr_files}
    
    merge_de(params.scripts_dir,diffexpr_files,diffexpr_methods)
        //.view{"Decoupler input: $it \n"}
        .set {mergede}
    
    func_decoupler(params.scripts_dir, mergede, resources)
        .groupTuple(by:[0,1])
        .map{it -> tuple(it[0],it[1],it[2].flatten())}
        .view{"Decoupler out: $it"}
        .set {decoupler}

    decoupler_merger(params.scripts_dir, decoupler)
        .set {merged_results}

    rank_analysis(params.scripts_dir, merged_results)
        .set {rank}

    rand_index_analysis(params.scripts_dir, merged_results)
        .set {rank}

}


params.scripts_dir = projectDir

//Downloads and stores prior knowledge sources
process get_prsources{
    cpus 2
    memory '32 GB'
	time '10h'
    executor 'slurm'

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
    cpus 16
    memory '32 GB'
	time '10h'
    executor 'slurm'

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
    cpus 16
    memory '32 GB'
	time '10h'
    executor 'slurm'

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
    cpus 16
    memory '32 GB'
	time '10h'
    executor 'slurm'
 
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
    cpus 2
    memory '32 GB'
	time '10h'
    executor 'slurm'

    input:
    path scripts_dir
    tuple val(datasetID), val(status), path (decoupler_results)

    output:
    tuple val(datasetID), path ("*__result.tsv")

    //afterScript "rm -f ${decoupler_results}"

    script:

    """
    Rscript ${scripts_dir}/decoupler_merger.R --dataset ${datasetID} --file ${decoupler_results} --status ${status}
    """
}

process subset_merger{
    cpus 2
    memory '32 GB'
	time '10h'
    executor 'slurm'

    publishDir "$params.scripts_dir/results/fullmerged/", mode: 'copy'

    input:
    path scripts_dir
    tuple val(datasetID), path (subset_files)

    output:
    tuple val(datasetID), path("*fullmerge.tsv")

    script:

    """
    Rscript ${scripts_dir}/subset_merger.R --dataset ${datasetID} --files "${subset_files}"
    """
}

//Rank analysis
process rank_analysis{
    cpus 2
    memory '32 GB'
	time '10h'
    executor 'slurm'

    publishDir "$params.scripts_dir/results/rank", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), path (analysis_results)

    output:
    tuple val(datasetID), path ("*__rank.tsv")

    script:

    """
    Rscript ${scripts_dir}/rank_analysis.R --dataset ${datasetID} --file ${analysis_results}
    """
}

//Rand index analysis
process rand_index_analysis{
    cpus 2
    memory '32 GB'
	time '10h'
    executor 'slurm'

    publishDir "$params.scripts_dir/results/rand_index", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), path (analysis_results)

    output:
    tuple val(datasetID), path ("*__randindex.tsv")

    script:

    """
    Rscript ${scripts_dir}/rand_index_analysis.R --dataset ${datasetID} --file ${analysis_results}
    """
}


workflow {
    Channel
        .fromFilePairs("$params.scripts_dir/data/*/*_{*_countdata,*_metadata}.tsv")
        .map{it -> tuple it[1][0].parent.baseName, it[0], it[1][0], it[1][1]}
        //.view()
        .set {datasets}
    
    Channel
        .of('vsn_norm limma_analysis', 'voom_norm limma_analysis', 'tmm_norm limma_analysis', 'log2quant_norm limma_analysis', 'edger_analysis', 'deseq2_analysis')
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
        .collectFile() { it ->
        [ "${it[0]}__${it[1]}.txt", "${it[2].parent}/${it[2].name}\n"]
        }
        .map{it -> tuple it.baseName.toString().replaceAll(/.txt/, "").split("__")[0], it.baseName.toString().replaceAll(/.txt/, "").split("__")[1], it}
        //.view()
        .set {decoupler}

    decoupler_merger(params.scripts_dir, decoupler)
        .map{it -> tuple it[0].split("_")[0], it[1]}
        .groupTuple(by:0)
        .set {subset_results}

    subset_merger(params.scripts_dir, subset_results)
        .set {full_results}

    rank_analysis(params.scripts_dir, full_results)
        .set {rank}

    rand_index_analysis(params.scripts_dir, full_results)
        .set {rank}

}

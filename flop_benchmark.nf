params.scripts_dir = projectDir
params.data_folder = "$params.scripts_dir/data"
params.parent_folder = projectDir
params.perturbation = ""

//Downloads and stores prior knowledge sources
process get_prsources{
    publishDir "$params.scripts_dir/scripts/dc_resources", mode: 'copy'

    input:
    path scripts_dir

    output:
    path ("*__source.tsv")

    script:

    """
    python3 ${scripts_dir}/scripts/get_resources_dc.py
    """
}

//Performs filtering, normalisation and differential expression analysis
process contrast_creator{
    input:
    path scripts_dir
    tuple val(subsetID), path(subset_dir)
 
    output:
    path("*.rds")
    
    script:

    """
    Rscript ${scripts_dir}/scripts/contrast_creator.R --file_dir ${subset_dir}
    """
}

//Performs filtering, normalisation and differential expression analysis
process diffexp_analysis{
    input:
    path scripts_dir
    tuple val(subsetID), val(biocontext), path(counts), path(meta)
    each pipelines
    each status
 
    output:
    tuple val(subsetID), val(biocontext), val(status), path("*__de.rds")
    
    script:

    """
    Rscript ${scripts_dir}/scripts/diffexp_analysis.R --dataset ${subsetID} --counts "${counts}" --meta "${meta}" --pipeline "${pipelines}" --status ${status} --bio ${biocontext}
    """
}

//Merges the results from different pipelines
process merge_de{

    input:
    path scripts_dir
    tuple val(subsetID), val(biocontext), val(status), path(diffexpr_files)
    each method
 
    output:
    tuple val(subsetID), val(biocontext), val(status), path ('*__decouplerinput.tsv')

    //afterScript "rm -f ${diffexpr_files}"
    
    script:

    """
    Rscript ${scripts_dir}/scripts/merge_de.R --dataset ${subsetID} --bio ${biocontext} --param ${method} --files "${diffexpr_files}"
    """
}

//Functional analysis 
process func_decoupler{

    input:
    path scripts_dir
    tuple val(subsetID), val(biocontext), val(status), path(decoupler_files)
    each resources
 
    output:
    tuple val(subsetID), val(status), path ('*__decoupleroutput.tsv')

    //afterScript "rm -f ${decoupler_files}"
    
    script:

    """
    python3 ${scripts_dir}/scripts/decoupler_proc.py ${decoupler_files} ${resources}
    """
}

//Merges the results from DecoupleR
process decoupler_merger{

    input:
    path scripts_dir
    tuple val(subsetID), val(status), path (decoupler_results)

    output:
    tuple val(subsetID), path ("*__result.tsv")

    //afterScript "rm -f ${decoupler_results}"

    script:

    """
    Rscript ${scripts_dir}/scripts/decoupler_merger.R --dataset ${subsetID} --file ${decoupler_results} --status ${status}
    """
}

process subset_merger{
    publishDir "$params.parent_folder/flop_benchmark_results/fullmerged/", mode: 'copy'

    input:
    path scripts_dir
    tuple val(datasetID), path (subset_files)

    output:
    tuple val(datasetID), path("*fullmerge.tsv")

    script:

    """
    Rscript ${scripts_dir}/scripts/subset_merger.R --dataset ${datasetID} --files "${subset_files}"
    """
}

//Rank analysis
process rank_analysis{
    publishDir "$params.parent_folder/flop_benchmark_results/rank", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), path (analysis_results)

    output:
    tuple val(datasetID), path ("*__rank.tsv")

    script:

    """
    Rscript ${scripts_dir}/scripts/rank_analysis.R --dataset ${datasetID} --file ${analysis_results}
    """
}

//Rand index analysis
process rand_index_analysis{
    publishDir "$params.parent_folder/flop_benchmark_results/rand_index", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), path (analysis_results)
    each perturbation_dataset

    output:
    tuple val(datasetID), path ("*__randindex.tsv")

    when: perturbation_dataset =~ datasetID

    script:

    """
    Rscript ${scripts_dir}/scripts/rand_index_analysis.R --dataset ${datasetID} --file ${analysis_results}
    """
}

//Rand index analysis
process jaccard_analysis{
    publishDir "$params.parent_folder/flop_benchmark_results/jaccard", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), path (analysis_results)

    output:
    tuple val(datasetID), path ("*__jaccard.tsv")

    script:

    """
    Rscript ${scripts_dir}/scripts/jaccard_analysis.R --dataset ${datasetID} --file ${analysis_results}
    """
}


workflow {
    Channel
        .fromPath("$params.data_folder/*", type: 'dir')
        .map{it -> tuple it[-1], it}
        // .view()
        .set {datasets}
    
    Channel
        .of('vsn_norm limma_analysis', 'voom_norm limma_analysis', 'tmm_norm limma_analysis', 'log2quant_norm limma_analysis', 'edger_analysis', 'deseq2_analysis')
        .set {pipelines}
    
    Channel
        .of('filtered', 'unfiltered')
        .set {status}
  
    Channel
        .of('logFC', 'stat')
        // .view()
        .set {diffexpr_methods}
    
    Channel
        .of(params.perturbation)
        .map{it -> it.split(" ")}
        .flatten()
        // .view{"Perturbation datasets: ${it}"}
        .set{perturbation_datasets}
    
    get_prsources(params.scripts_dir)
        .flatten()
        .map{it -> it.baseName.toString().replaceAll(/__source/, "")}
        .set{resources}
    
    contrast_creator(params.scripts_dir, datasets)
        .flatten()
        // .view()
        .map{it -> tuple (
            it.name.toString().split("__")[0],
            it.name.toString().split("__")[1],
            it)}
        // .view()
        .groupTuple(by:[0,1])
        .map{it -> tuple it[0], it[1], it[2][0], it[2][1]}
        // .view()
        .set {contrasts}

    diffexp_analysis(params.scripts_dir, contrasts, pipelines, status)
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

    rand_index_analysis(params.scripts_dir, full_results, perturbation_datasets)
        .set {randindex}

    jaccard_analysis(params.scripts_dir, full_results)
        .set {jaccard}

}
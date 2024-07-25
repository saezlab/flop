params.scripts_dir = projectDir
params.data_folder = "$params.scripts_dir/data"
params.parent_folder = projectDir

def asciiBanner = 
"##########################################################################\n" +
" Welcome to\n" +
".------..------..------..------.\n" +
"|F.--. ||L.--. ||O.--. ||P.--. |\n" +
"| :(): || :/\\: || :/\\: || :/\\: |\n" +
"| ()() || (__) || :\\/: || (__) |\n" +
"|  -- F||  -- L||  -- O||  -- P|\n" +
"'------''------''------''------'\n" +
"\n" +
" The FunctionaL Omics Preprocessing platform is\n" +
" a workflow meant to evaluate the impact of different\n" +
" normalization and differential expression tools on the\n" +
" resulting functional space, in the context of bulk RNA-seq data.\n" +
"\n" +
" ##########################################################################\n"

println asciiBanner

//Downloads and stores prior knowledge sources
process get_prsources{
    storeDir "$params.scripts_dir/scripts/dc_resources"

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
    path("*.qs")
    
    script:

    """
    Rscript ${scripts_dir}/scripts/contrast_creator.R --file_dir ${subset_dir}
    """

}

//Performs filtering, normalisation and differential expression analysis
process diffexp_analysis {

    input:
    path scripts_dir
    tuple val(subsetID), val(biocontext), path(counts), path(meta)
    each pipelines
    each status

    output:
    tuple val(subsetID), val(biocontext), val(status), path("*__de.qs")

    script:
    """
    Rscript ${scripts_dir}/scripts/diffexp_analysis.R \
    --dataset ${subsetID} \
    --counts "${counts}" \
    --meta "${meta}" \
    --pipeline "${pipelines}" \
    --status ${status} \
    --bio ${biocontext} \
    --filterbyexpr-libsize ${params.filterbyexpr_libsize} \
    --filterbyexpr-mincount ${params.filterbyexpr_mincount} \
    --filterbyexpr-mintotalcount ${params.filterbyexpr_mintotalcount} \
    --filterbyexpr-largen ${params.filterbyexpr_largen} \
    --filterbyexpr-minprop ${params.filterbyexpr_minprop} \
    --diffexp-limma-ndups ${params.diffexp_limma_ndups} \
    --diffexp-limma-spacing ${params.diffexp_limma_spacing} \
    --diffexp-limma-block "${params.diffexp_limma_block}" \
    --diffexp-limma-weights ${params.diffexp_limma_weights} \
    --diffexp-limma-method ${params.diffexp_limma_method} \
    --diffexp-deseq2-test ${params.diffexp_deseq2_test} \
    --diffexp-deseq2-fitType ${params.diffexp_deseq2_fitType} \
    --diffexp-deseq2-quiet ${params.diffexp_deseq2_quiet} \
    --diffexp-deseq2-minReplicatesForReplace ${params.diffexp_deseq2_minReplicatesForReplace} \
    --diffexp-deseq2-parallel ${params.diffexp_deseq2_parallel} \
    --diffexp-deseq2-betaprior ${params.diffexp_deseq2_betaprior} \
    --diffexp-edger-calcnormfactors-method ${params.diffexp_edger_calcnormfactors_method} \
    --diffexp-edger-calcnormfactors-refColumn ${params.diffexp_edger_calcnormfactors_refColumn} \
    --diffexp-edger-calcnormfactors-logratiotrim ${params.diffexp_edger_calcnormfactors_logratiotrim} \
    --diffexp-edger-calcnormfactors-sumtrim ${params.diffexp_edger_calcnormfactors_sumtrim} \
    --diffexp-edger-calcnormfactors-doweighting ${params.diffexp_edger_calcnormfactors_doweighting} \
    --diffexp-edger-calcnormfactors-acutoff ${params.diffexp_edger_calcnormfactors_acutoff} \
    --diffexp-edger-calcnormfactors-p ${params.diffexp_edger_calcnormfactors_p} \
    --diffexp-edger-estimatedisp-priordf ${params.diffexp_edger_estimatedisp_priordf} \
    --diffexp-edger-estimatedisp-trendmethod ${params.diffexp_edger_estimatedisp_trendmethod} \
    --diffexp-edger-estimatedisp-tagwise ${params.diffexp_edger_estimatedisp_tagwise} \
    --diffexp-edger-estimatedisp-mixeddf ${params.diffexp_edger_estimatedisp_mixeddf} \
    --diffexp-edger-estimatedisp-span ${params.diffexp_edger_estimatedisp_span} \
    --diffexp-edger-estimatedisp-minrowsum ${params.diffexp_edger_estimatedisp_minrowsum} \
    --diffexp-edger-estimatedisp-gridlength ${params.diffexp_edger_estimatedisp_gridlength} \
    --diffexp-edger-estimatedisp-gridrange "${params.diffexp_edger_estimatedisp_gridrange}" \
    --diffexp-edger-estimatedisp-robust ${params.diffexp_edger_estimatedisp_robust} \
    --diffexp-edger-estimatedisp-winsortailp "${params.diffexp_edger_estimatedisp_winsortailp}" \
    --diffexp-edger-estimatedisp-tol ${params.diffexp_edger_estimatedisp_tol} \
    --diffexp-edger-glmqlfit-dispersion ${params.diffexp_edger_glmqlfit_dispersion} \
    --diffexp-edger-glmqlfit-libsize ${params.diffexp_edger_glmqlfit_libsize} \
    --diffexp-edger-glmqlfit-offset ${params.diffexp_edger_glmqlfit_offset} \
    --diffexp-edger-glmqlfit-weights ${params.diffexp_edger_glmqlfit_weights} \
    --diffexp-edger-glmqlfit-abundancetrend ${params.diffexp_edger_glmqlfit_abundancetrend} \
    --diffexp-edger-glmqlfit-avelogcpm ${params.diffexp_edger_glmqlfit_avelogcpm} \
    --diffexp-edger-glmqlfit-robust ${params.diffexp_edger_glmqlfit_robust}
    """
}


//Merges the output of differential expression analysis files
process output_merge_de{

    publishDir "$params.parent_folder/flop_results/diffexp", mode: 'copy'

    input:
    path scripts_dir
    tuple val(datasetID), path(diffexpr_files)


    output:
    tuple val(datasetID), path ('*__deresults.tsv')

    //afterScript "rm -f ${diffexpr_files}"
    
    script:

    """
    Rscript ${scripts_dir}/scripts/diffexp_merger.R --dataset ${datasetID} --files "${diffexpr_files}"
    """

}

//Merges the results from different pipelines
process downstream_merge_de{

    input:
    path scripts_dir
    tuple val(subsetID), val(biocontext), val(status), path(diffexpr_files)
    each method
    val(flop_ngenes_threshold)
 
    output:
    tuple val(subsetID), val(biocontext), val(status), path ('*__decouplerinput.tsv'), optional: true

    //afterScript "rm -f ${diffexpr_files}"
    
    script:

    """
    Rscript ${scripts_dir}/scripts/merge_de.R --dataset ${subsetID} --bio ${biocontext} --param ${method} --files "${diffexpr_files}" --threshold ${flop_ngenes_threshold}
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
    python3 ${scripts_dir}/scripts/decoupler_proc.py \
        ${decoupler_files} \
        ${resources} \
        ${params.decoupler_runulm_batch_size} \
        ${params.decoupler_runulm_minn} \
        ${params.decoupler_runulm_verbose} \
        ${params.decoupler_runulm_useraw}
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

    publishDir "$params.parent_folder/flop_results/funcomics/fullmerged/", mode: 'copy'

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

    publishDir "$params.parent_folder/flop_results/funcomics/rank", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), path (func_results), path (de_results)

    output:
    tuple val(datasetID), path ("*__rank.tsv")

    script:

    """
    Rscript ${scripts_dir}/scripts/rank_analysis.R --dataset ${datasetID} --func_file ${func_results} --de_file ${de_results}
    """

}

//Top/bottom features overlap analysis
process top_bottom_overlap_analysis{

    publishDir "$params.parent_folder/flop_results/funcomics/overlap", mode: 'move'

    input:
    path scripts_dir
    tuple val(datasetID), path (func_results), path (de_results)
    val(flop_pval_threshold)

    output:
    tuple val(datasetID), path ("*__overlap.tsv")

    script:

    """
    Rscript ${scripts_dir}/scripts/top_bottom_overlap_analysis.R --dataset ${datasetID} --func_file ${func_results} --de_file ${de_results} --pval_thresh ${flop_pval_threshold}
    """

}


workflow {
    Channel
        .fromPath("$params.data_folder/*", type: 'dir')
        .map{it -> tuple it[-1], it}
        .set {datasets}
    
    Channel
        .of(params.diffexp_pipelines.split(","))
        .set {pipelines}
    
    Channel
        .of(params.filtering.split(","))
        .set {status}
  
    Channel
        .of('logFC', 'stat')
        .set {diffexpr_metrics}
    
    get_prsources(params.scripts_dir)
        .flatten()
        .map{it -> it.baseName.toString().replaceAll(/__source/, "")}
        .set{resources}
    
    contrast_creator(params.scripts_dir, datasets)
        .flatten()
        .map{it -> tuple (
            it.name.toString().split("__")[0],
            it.name.toString().split("__")[1],
            it)}
        .groupTuple(by:[0,1], size: 2)
        .map{it -> tuple it[0], it[1], it[2][0], it[2][1]}
        .set {contrasts}

    diffexp_analysis(params.scripts_dir, contrasts, pipelines, status)
        .multiMap { it -> downstream: output: it}
        .set {diffexpr}
    
    diffexpr.downstream
        .groupTuple(by:[0,1,2], size: params.diffexp_pipelines.length())
        .set {diffexpr_files}
    
    diffexpr.output
        .collectFile() { it ->
        [ "${it[0]}__defiles.txt", "${it[3]}\n"]
        }
        .map{it -> tuple it.baseName.toString().replaceAll(/.txt/, "").split("__")[0].split("_")[0], it}
        .groupTuple(by:0)
        .set {diffexpr_out}
    
    output_merge_de(params.scripts_dir,diffexpr_out)
        .set {diffexpr_merged}

    if (params.functional_analysis == true){
        func_decoupler(params.scripts_dir, diffexpr_merged, resources)
            .collectFile() { it ->
            [ "${it[0]}__${it[1]}.txt", "${it[2].parent}/${it[2].name}\n"]
            }
            .map{it -> tuple it.baseName.toString().replaceAll(/.txt/, "").split("__")[0], it.baseName.toString().replaceAll(/.txt/, "").split("__")[1], it}
            .set {decoupler}
    
        decoupler_merger(params.scripts_dir, decoupler)
            .map{it -> tuple it[0].split("_")[0], it[1]}
            .groupTuple(by:0)
            .map{it -> tuple it[0], it[1]}
            .set {subset_results}
    
        subset_merger(params.scripts_dir, subset_results)
            .concat(diffexpr_merged)
            .groupTuple(by:0, size:2)
            .map{it -> tuple it[0], it[1][0], it[1][1]}
            .set {full_results}
    
        rank_analysis(params.scripts_dir, full_results)
            .set {rank}
    
        top_bottom_overlap_analysis(params.scripts_dir, full_results, params.flop_pval_threshold)
            .set {jaccard}
    }
}

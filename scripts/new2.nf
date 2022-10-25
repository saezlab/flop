
params.scripts_dir = "/mnt/c/Users/victo/Onedrive - Universidad Politécnica de Madrid/Documentos/1º Master/Internship/flop_benchmark/scripts"

//normalization
process normalize {
 
    input:
    path scripts_dir
    tuple val(dataID), path(datafiles)
    each norm_method
 
    output:
    tuple val(dataID),path ('*__norm.tsv')
    
    script:
    def (counts, meta) = datafiles

    """
    Rscript ${scripts_dir}/${norm_method}.R --counts ${counts}
    """
}

// Diff exp with limma
process diffexp_limma {
 
    input:
    path scripts_dir
    tuple val(dataID), path(norm_file),path(datafiles)

    output:
    tuple val(dataID), path ('*__limma__de.tsv')
    
    script:
    def (counts, meta) = datafiles

    """
    Rscript ${scripts_dir}/limma.R --norm ${norm_file} --meta ${meta}
    """

}

// Diff exp with DESeq2
process diffexp_deseq2 {
 
    input:
    path scripts_dir
    tuple val(dataID), path(datafiles)
 
    output:
    tuple val(dataID), path ('*__deseq2__de.tsv')
    
    script:
    def (counts, meta) = datafiles

    """
    Rscript ${scripts_dir}/deseq2.R --counts ${counts} --meta ${meta}
    """

}

// Diff exp with EdgeR
process diffexp_edger{
 
    input:
    path scripts_dir
    tuple val(dataID), path(datafiles)
 
    output:
    tuple val(dataID), path ('*__edger__de.tsv')
    
    script:
    def (counts, meta) = datafiles

    """
    Rscript ${scripts_dir}/edger.R --counts ${counts} --meta ${meta}
    """

}

//Differential expression output file merger
process merge_de{
 
    input:
    path scripts_dir
    tuple val(dataID), path(diffexpr_files)
    each method
 
    output:
    path ('*__decouplerinput.tsv')
    
    script:

    """
    Rscript ${scripts_dir}/merge_de.R --dataset ${dataID} --param ${method}
    """

}

//Functional analysis 
process func_decoupler{
 
    input:
    path scripts_dir
    path(decoupler_files)
 
    output:
    path ('*__decoupleroutput.tsv')
    
    script:

    """
    Rscript ${scripts_dir}/decoupler.R --input ${decoupler_files}
    """

}


workflow {
    Channel
        .fromFilePairs('/mnt/c/Users/victo/Onedrive - Universidad Politécnica de Madrid/Documentos/1º Master/Internship/flop_benchmark/scripts/data/*_{GeneLevel_Raw_data,*_metadata}.tsv')
        .set {datasets}
    
    Channel
        .of('padj', 'logFC', 'stat')
        .set {diffexpr_methods}
    
    Channel
        .of('vsn', 'log2quant', 'tmm')
        .set{norm_methods}
    
    //normalization channels
    normalize(params.scripts_dir,datasets, norm_methods)
        .combine(datasets, by: 0)
        .set{normalised_vals}
        

    


    //differential analysis channels
    diffexp_deseq2(params.scripts_dir,datasets).set {deseq2}
    diffexp_edger(params.scripts_dir,datasets).set {edger}
    
    diffexp_limma(params.scripts_dir,normalised_vals)
        .set{limma}

    //functional analysis channels
    limma
        .mix(deseq2, edger)
        .groupTuple()
        .set {diffexpr_files}
    
    merge_de(params.scripts_dir,diffexpr_files,diffexpr_methods)
        .view()
        .set {mergede}
    //func_decoupler(params.scripts_dir, mergede).set {decoupler}
}
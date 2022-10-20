
params.scripts_dir = "/mnt/c/Users/victo/Onedrive - Universidad Politécnica de Madrid/Documentos/1º Master/Internship/flop_benchmark/scripts"

process normalize_vsn {
 
    input:
    path scripts_dir
    tuple val(dataID), path(datafiles)
 
    output:
    tuple val(dataID), path ('*__vsn__norm.tsv')
    
    script:
    def (counts, meta) = datafiles

    """
    Rscript ${scripts_dir}/vsn.R --counts ${counts}
    """
}

process normalize_tmm {
 
    input:
    path scripts_dir
    tuple val(dataID), path(datafiles)
 
    output:
    tuple val(dataID), path ('*__tmm__norm.tsv')
    
    script:
    def (counts, meta) = datafiles

    """
    Rscript ${scripts_dir}/tmm.R --counts ${counts}
    """
}

process normalize_log2quant {
 
    input:
    path scripts_dir
    tuple val(dataID), path(datafiles)
 
    output:
    tuple val(dataID), path ('*__log2quant__norm.tsv')
    
    script:
    def (counts, meta) = datafiles

    """
    Rscript ${scripts_dir}/log2quant.R --counts ${counts}
    """
}

// Diff exp with limma
process diffexp_limma {
 
    input:
    path scripts_dir
    tuple val(dataID_norm), path(norm_files)
    tuple val(dataID_data), path(datafiles)
    
    when:
    dataID_norm == dataID_data

    output:
    tuple val(dataID_norm), path ('*__limma__de.tsv')
    
    script:
    def (counts, meta) = datafiles

    """
    Rscript ${scripts_dir}/limma.R --norm ${norm_files} --meta ${meta}
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

process merge_de{
 
    input:
    path scripts_dir
    tuple val(dataID), path(diffexpr_files)
 
    output:
    path ('*__decouplerinput.tsv')
    
    script:

    """
    Rscript ${scripts_dir}/merge_de.R --dataset ${dataID} --param padj
    """

}


workflow {
    Channel
        .fromFilePairs('/mnt/c/Users/victo/Onedrive - Universidad Politécnica de Madrid/Documentos/1º Master/Internship/flop_benchmark/scripts/data/*_{GeneLevel_Raw_data,*_metadata}.tsv')
        .set {files}
    
    normalize_vsn(params.scripts_dir,files).set {vsn}
    normalize_tmm(params.scripts_dir,files).set {tmm}
    normalize_log2quant(params.scripts_dir,files).set {log2quant}

    vsn
        .mix(tmm, log2quant)
        .groupTuple()
        .set {norm_files}
    
    diffexp_limma(params.scripts_dir,norm_files,files).set {limma}
    diffexp_deseq2(params.scripts_dir,files).set {deseq2}
    diffexp_edger(params.scripts_dir,files).set {edger}

    limma
        .mix(deseq2, edger)
        .groupTuple()
        .set {diffexpr_files}
    diffexpr_files.view()
    merge_de(params.scripts_dir,diffexpr_files).set {mergede}
    mergede.view()
}
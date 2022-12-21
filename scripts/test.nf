params.scripts_dir = projectDir

process get_prsources{
    publishDir "$params.scripts_dir/dc_resources", mode: 'move'

    input:
    path scripts_dir

    output:
    path ("*__source.tsv")

    script:

    """
    python3 ${scripts_dir}/get_resources_dc.py
    """
}

//Database selection
workflow {    
    get_prsources(params.scripts_dir)
        .map{it -> it.baseName.toString()}
        .map{it -> it.replaceAll(/__source/, "")}
        .view()
        .set {pr_sources}

}
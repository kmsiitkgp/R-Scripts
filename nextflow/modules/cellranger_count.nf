process CELLRANGER_COUNT {

    tag "Running CellRanger count on ${sample_id}"
    label 'process_high'

    input:
    val(sample_id)
    path(raw_fastq_dir)
    path(cellranger_index_dir)
    val(cellranger_args)

    output:
    // Matrices
    path("${sample_id}/outs/raw_feature_bc_matrix"),         emit: raw_matrix_dir
    path("${sample_id}/outs/filtered_feature_bc_matrix"),     emit: filt_matrix_dir

    // Reports & Metrics
    path("${sample_id}/outs/web_summary.html"),                emit: web_summary
    path("${sample_id}/outs/metrics_summary.csv"),            emit: metric_summary
    path("${sample_id}/outs/*.json"),                         emit: run_info

    // Data files
    path("${sample_id}/outs/*.h5"),                         emit: h5_files
    path("${sample_id}/outs/*.cloupe"),                     emit: cloupe
    path("${sample_id}/outs/*.bam*"),                         emit: bam_files,             optional: true

    // Logs and Info
    path("${sample_id}/outs/analysis"),                     emit: analysis_dir,            optional: true
    path("${sample_id}.CELLRANGER_COUNT.error.log"),        emit: cellranger_count_error_log  // Process log

    script:

    def LOG = "${sample_id}.CELLRANGER_COUNT.error.log"

    """

    cellranger count \
        ${cellranger_args} \
        --transcriptome "${cellranger_index_dir}" \
        --fastqs "${raw_fastq_dir}" \
        --sample "${sample_id}" \
        --id "${sample_id}" \
        --localcores "${task.cpus}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: CellRanger count failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: CellRanger count completed for ${sample_id}" >> "${LOG}"
    """
}


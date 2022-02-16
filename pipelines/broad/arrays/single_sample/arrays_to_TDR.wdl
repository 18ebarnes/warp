version 1.0

import "Arrays.wdl" as ArraysPipeline

workflow arrays_outputs_to_TDR {
    meta {
        description: "Push outputs of Arrays.wdl to TDR dataset table ArraysOutputsTable."
    }

    input {
        # inputs to wrapper task 
        String workspace_name
        String workspace_project
        String workspace_project
        String gcp_project
        String tdr_dataset_id
        String tdr_target_table_name
        # String workflow_name # Arrays in this example

        # required inputs to Arrays.wdl
        String chip_well_barcode
        String sample_alias
        String sample_lsid
        String reported_gender
        File red_idat_cloud_path
        File green_idat_cloud_path
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File dbSNP_vcf
        File dbSNP_vcf_index
        File haplotype_database_file
        File variant_rsids_file
        Int disk_size
        Int preemptible_tries
        String environment
        File vault_token_path
    }

    call ArraysPipeline.Arrays {
        input:
            chip_well_barcode          = chip_well_barcode,
            sample_alias               = sample_alias,
            sample_lsid                = sample_lsid,
            reported_gender            = reported_gender,
            red_idat_cloud_path        = red_idat_cloud_path,
            green_idat_cloud_path      = green_idat_cloud_path, 
            ref_fasta                  = ref_fasta, 
            ref_fasta_index            = ref_fasta_index, 
            ref_dict                   = ref_dict, 
            dbSNP_vcf                  = dbSNP_vcf, 
            dbSNP_vcf_index            = dbSNP_vcf_index, 
            haplotype_database_file    = haplotype_database_file,
            variant_rsids_file         = variant_rsids_file, 
            disk_size                  = disk_size, 
            preemptible_tries          = preemptible_tries, 
            environment                = environment, 
            vault_token_path           = vault_token_path
    }

    call format_arrays_outputs {
        input:
            chip_well_barcode_output                            = Arrays.chip_well_barcode_output,
            analysis_version_number_output                      = Arrays.analysis_version_number_output,
            gtc_file                                            = Arrays.gtc_file,
            # red_idat_md5_cloud_path                             = Arrays.red_idat_md5_cloud_path,
            # green_idat_md5_cloud_path                           = Arrays.green_idat_md5_cloud_path,
            # output_vcf_md5_cloud_path                           = Arrays.output_vcf_md5_cloud_path,
            output_vcf                                          = Arrays.output_vcf,
            output_vcf_index                                    = Arrays.output_vcf_index,
            baf_regress_metrics_file                            = Arrays.baf_regress_metrics_file,
            # contamination_metrics_file                          = Arrays.contamination_metrics_file,
            # reference_fingerprint_vcf                           = Arrays.reference_fingerprint_vcf,
            # reference_fingerprint_vcf_index                     = Arrays.reference_fingerprint_vcf_index,
            # output_fingerprint_vcf                              = Arrays.output_fingerprint_vcf,
            # output_fingerprint_vcf_index                        = Arrays.output_fingerprint_vcf_index,
            # output_fingerprint_json_file                        = Arrays.output_fingerprint_json_file,
            arrays_variant_calling_detail_metrics_file          = Arrays.arrays_variant_calling_detail_metrics_file,
            arrays_variant_calling_summary_metrics_file         = Arrays.arrays_variant_calling_summary_metrics_file,
            arrays_variant_calling_control_metrics_file         = Arrays.arrays_variant_calling_control_metrics_file,
            # arrays_subset_variant_calling_detail_metrics_file   = Arrays.arrays_subset_variant_calling_detail_metrics_file,
            # arrays_subset_variant_calling_summary_metrics_file  = Arrays.arrays_subset_variant_calling_summary_metrics_file,
            # arrays_subset_variant_calling_control_metrics_file  = Arrays.arrays_subset_variant_calling_control_metrics_file,
            fingerprint_detail_metrics_file                     = Arrays.fingerprint_detail_metrics_file,
            fingerprint_summary_metrics_file                    = Arrays.fingerprint_summary_metrics_file,
            genotype_concordance_summary_metrics_file           = Arrays.genotype_concordance_summary_metrics_file,
            genotype_concordance_detail_metrics_file            = Arrays.genotype_concordance_detail_metrics_file,
            genotype_concordance_contingency_metrics_file       = Arrays.genotype_concordance_contingency_metrics_file
            # chip_well_barcode_params_file                       = Arrays.chip_well_barcode_params_file
    }

    call ingest_outputs_to_tdr {
        input:
            workspace_name          = workspace_name,
            workspace_project       = workspace_project,
            workspace_bucket        = workspace_bucket,
            gcp_project             = gcp_project,
            tdr_dataset_id          = tdr_dataset_id,
            tdr_target_table_name   = tdr_target_table_name,
            ingest_json             = format_arrays_outputs.ingest_outputs_json
    }


}

task format_arrays_outputs {
    input {
        String  chip_well_barcode_output
        Int     analysis_version_number_output
        String? baf_regress_metrics_file
        String  gtc_file

        String? output_vcf
        String? output_vcf_index

        String  arrays_variant_calling_detail_metrics_file
        String? arrays_variant_calling_summary_metrics_file
        String? arrays_variant_calling_control_metrics_file

        String? fingerprint_detail_metrics_file
        String? fingerprint_summary_metrics_file

        String? genotype_concordance_summary_metrics_file
        String? genotype_concordance_detail_metrics_file
        String? genotype_concordance_contingency_metrics_file

        # this value is added during the python script that recodes json
        #String last_modified_date
    }

    command <<<
        echo -e "chip_well_barcode_output\tanalysis_version_number_output\tbaf_regress_metrics_file\tgtc_file\t\
                 output_vcf\toutput_vcf_index\t\
                 arrays_variant_calling_detail_metrics_file\tarrays_variant_calling_summary_metrics_file\tarrays_variant_calling_control_metrics_file\t\
                 fingerprint_detail_metrics_file\tfingerprint_summary_metrics_file\t\
                 genotype_concordance_summary_metrics_file\tgenotype_concordance_detail_metrics_file\tgenotype_concordance_contingency_metrics_file" \
                 > ingestDataset_arrays_outputs.tsv
        
        echo -e "~{chip_well_barcode_output}\t~{analysis_version_number_output}\t~{baf_regress_metrics_file}\t~{gtc_file}\t\
                 ~{output_vcf}\t~{output_vcf_index}\t\
                 ~{arrays_variant_calling_detail_metrics_file}\t~{arrays_variant_calling_summary_metrics_file}\t~{arrays_variant_calling_control_metrics_file}\t\
                 ~{fingerprint_detail_metrics_file}\t~{fingerprint_summary_metrics_file}\t\
                 ~{genotype_concordance_summary_metrics_file}\t~{genotype_concordance_detail_metrics_file}\t~{genotype_concordance_contingency_metrics_file}" \
                 >> ingestDataset_arrays_outputs.tsv

        python3 << CODE
            import pandas as pd

            tsv_df = pd.read_csv("ingestDataset_arrays_outputs.tsv", sep="\t")
            tsv_df = tsv_df.dropna(axis=1, how="all")  # drop columns if no value (optional outputs etc)

            outputs = tsv_df.to_json("ingestDataset_arrays_outputs.json", orient="records")  # write json file

        CODE
    >>>

    runtime {

    }

    output {
        File ingest_outputs_tsv = "ingestDataset_arrays_outputs.tsv"
        File ingest_outputs_json = "ingestDataset_arrays_outputs.json"
    }
}

task ingest_outputs_to_tdr {
    input {
        String workspace_name
        String workspace_project
        String workspace_bucket
        String gcp_project
        # String workflow_name
        String tdr_dataset_id
        String tdr_target_table_name

        File   outputs_json
    }

    command {

        python3 WDL_write_arrays_wdl_outputs_to_TDR_ArraysOutputsTable.py -s submission_id \
                                                                          -p ~{workspace_project} \
                                                                          -w ~(workspace_name) \
                                                                          -g ~{gcp_project} \
                                                                          -b ~{workspace_bucket} \
                                                                          -d ~{tdr_dataset_id} \
                                                                          -t ~{tdr_target_table_name} \
                                                                          -j ~{outputs_json}
    }

    runtime {

    }

    output {

    }
}
version 1.0

workflow ExtractOptimizeSingleSample { 
    input {
        String monitoring_script="gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
        File input_vcf
        File input_vcf_index
        String base_file_name
        String sample_name_calls
        File gtr_vcf
        File gtr_vcf_index
        File gtr_highconf_intervals
        String sample_name_gtr
        String gatk_docker
        String jukebox_vc_docker
        String flow_order
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_fasta_sdf
        File runs_file
        Array[File] annotation_intervals
        String score_key
    }

    call AnnotateSampleVCF {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            docker = gatk_docker,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            output_basename = base_file_name
    }
    
    call ExtractSample { 
        input: 
            sample_name = sample_name_calls,
            monitoring_script = monitoring_script,
            input_vcf = AnnotateSampleVCF.output_vcf_file,
            input_vcf_index = AnnotateSampleVCF.output_vcf_index,
            docker=jukebox_vc_docker,
            no_address=true,

    }

    call FilterSampleVCF{
        input:
            monitoring_script = monitoring_script,
            input_vcf = ExtractSample.output_vcf_file,
            no_address=true,
            docker=jukebox_vc_docker,
    }

    call FilterSymbolicAlleles{
        input:
            monitoring_script = monitoring_script,
            base_file_name = base_file_name,
            docker=jukebox_vc_docker,
            input_vcf = FilterSampleVCF.output_vcf_file,
            input_vcf_index = FilterSampleVCF.output_vcf_index,
            no_address=true
         }

    call CompareToGroundTruth{
        input:
            gtr_vcf = gtr_vcf,
            gtr_vcf_index = gtr_vcf_index, 
            gtr_highconf = gtr_highconf_intervals,
            left_sample_name = sample_name_calls,
            right_sample_name = sample_name_gtr,
            input_vcf = FilterSymbolicAlleles.output_vcf,
            input_vcf_index = FilterSymbolicAlleles.output_vcf_index,
            input_vcf_name = base_file_name,
            
            ref_fasta = ref_fasta,
            ref_dict = ref_dict,
            ref_fasta_sdf = ref_fasta_sdf,
            runs_file = runs_file,
            annotation_intervals = annotation_intervals,

            disk_size = 2*(ceil(size(input_vcf, "GB") + size(gtr_vcf, "GB") + size(ref_fasta, "GB") + size(ref_fasta_sdf, "GB")))+10,
            docker = jukebox_vc_docker,
            monitoring_script = monitoring_script,
            no_address = true,
            flow_order = flow_order

    }

    call EvaluateResults {
        input:
          h5_input = CompareToGroundTruth.compare_h5,
          input_vcf_name = base_file_name,
          disk_size = 10,
          docker = jukebox_vc_docker,
          monitoring_script = monitoring_script,
          no_address = true
    }

    call HardThresholdVCF { 
        input: 
          thresholds = EvaluateResults.thresholds_report,
          input_vcf = input_vcf,
          input_vcf_index = input_vcf_index,
          score_key = score_key,
          output_basename = base_file_name,
          disk_size = 3*ceil(size(input_vcf, "GB")) + 14,
          docker = gatk_docker,
          monitoring_script = monitoring_script, 
          no_address = false
    }

    output {
        File output_vcf = HardThresholdVCF.output_vcf
        File output_vcf_index = HardThresholdVCF.output_vcf_index
        File eval_report_h5 = EvaluateResults.eval_report_h5
        File thresholds_report = EvaluateResults.thresholds_report
    }
}

task ExtractSample {
    input {
        File monitoring_script
        String sample_name
        File input_vcf
        File input_vcf_index
        String docker
        Boolean no_address
    }
    String output_vcf = basename(input_vcf, ".vcf.gz") + ".~{sample_name}.vcf.gz"
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} > monitoring.log &
        source ~/.bashrc
        conda activate genomics.py3
        bcftools view -s ~{sample_name} ~{input_vcf} -o  tmp.vcf.gz -O z 
        bcftools index -t tmp.vcf.gz
        bcftools view -h tmp.vcf.gz | sed 's/COMBINED_TREE_SCORE,Number=A,/COMBINED_TREE_SCORE,Number=\.,/g' > hdr.fixed.txt
        bcftools reheader -h hdr.fixed.txt tmp.vcf.gz | bcftools view -Oz -o ~{output_vcf} - 
        bcftools index -t ~{output_vcf}
    >>>
    output {
        File monitoring_log = "monitoring.log"
        File output_vcf_file = "~{output_vcf}"
        File output_vcf_index = "~{output_vcf}.tbi"
    }
    runtime {
        memory: "8GB"
        disks: "local-disk " + (ceil(size(input_vcf, "GB")) * 3 + 10) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 1
    }
}

task FilterSampleVCF{
    input{
        File monitoring_script
        File input_vcf
        String docker
        Boolean no_address
    }

    String output_vcf = basename(input_vcf, ".vcf.gz") + ".filter_unused.vcf.gz"

    command <<<
        bash ~{monitoring_script} > monitoring.log &
        source ~/.bashrc
        conda activate genomics.py3
        bcftools view -a -i 'GT[*]="alt"' ~{input_vcf} -o ~{output_vcf} -O z 
        bcftools index -t ~{output_vcf}
    >>>

    output {
        File monitoring_log = "monitoring.log"
        File output_vcf_file = "~{output_vcf}"
        File output_vcf_index = "~{output_vcf}.tbi"
    }

    runtime {
        memory: "8GB"
        disks: "local-disk " + (ceil(size(input_vcf, "GB")) *2 + 10) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 1
    }
}

task FilterSymbolicAlleles {
    input {
        File monitoring_script
        String base_file_name
        File input_vcf
        File input_vcf_index
        Boolean no_address
        String docker
    }

    String output_vcf_name = "~{base_file_name}" + ".vcf.gz"
    command <<<
        bash ~{monitoring_script} > monitoring.log &
        source ~/.bashrc
        conda activate genomics.py3
        gatk --java-options "-Xmx10g"  SelectVariants \
            -V ~{input_vcf} \
            -O ~{output_vcf_name}.tmp.vcf.gz \
            --remove-unused-alternates
        gatk --java-options "-Xmx10g" SelectVariants \
            -V ~{output_vcf_name}.tmp.vcf.gz \
            -O ~{output_vcf_name} \
            --exclude-non-variants \
            --select-type-to-exclude SYMBOLIC
        >>>
    runtime {
        memory: "12 GB"
        cpu: 1
        disks: "local-disk " + (ceil(size(input_vcf, "GB")) *4 +10) + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File output_vcf = "~{output_vcf_name}"
        File output_vcf_index = "~{output_vcf_name}.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task CompareToGroundTruth {
  input {
    File monitoring_script
    String left_sample_name
    String right_sample_name
    File input_vcf
    File input_vcf_index

    String input_vcf_name

    File gtr_vcf
    File gtr_vcf_index
    File gtr_highconf

    File? interval_list
    File runs_file
    File ref_fasta
    File ref_dict
    File ref_fasta_sdf
    String flow_order
    Array[File] annotation_intervals
    Int disk_size
    String docker
    Boolean no_address
  }

  String used_flow_order = (if flow_order=="" then "TACG" else flow_order)

  command <<<
    bash ~{monitoring_script} > monitoring.log &

    source ~/.bashrc
    conda activate genomics.py3

    python -m tarfile -e ~{ref_fasta_sdf} ~{ref_fasta}.sdf

    run_comparison_pipeline.py \
            --n_parts 0 \
            --hpol_filter_length_dist 12 10 \
            --input_prefix $(echo "~{input_vcf}" | sed 's/\(.vcf.gz\|.vcf\)$//') \
            --output_file ~{input_vcf_name}.comp.h5 \
            --gtr_vcf ~{gtr_vcf} \
            --highconf_intervals ~{gtr_highconf} \
            --runs_intervals ~{runs_file} \
            --reference ~{ref_fasta} \
            --reference_dict ~{ref_dict} \
            --call_sample_name ~{left_sample_name} \
            --ignore_filter_status \
            --flow_order ~{used_flow_order} \
            --truth_sample_name ~{right_sample_name}\
            --annotate_intervals ~{sep=" --annotate_intervals " annotation_intervals} \
            --n_jobs 8 \
            --output_suffix '' \
            --output_interval interval_~{input_vcf_name}.comp.bed
    >>>
  runtime {
    memory: "32 GB"
    cpu: 16
    disks: "local-disk " + disk_size + " SSD"
    docker: docker
    noAddress: no_address
  }
  output {
    File compare_h5 = "~{input_vcf_name}.comp.h5"
    File monitoring_log = "monitoring.log"
    Array[File] comparison_beds = glob("~{input_vcf_name}.*.bed")
    File output_interval = "interval_~{input_vcf_name}.comp.bed"
  }
}

task EvaluateResults {
  input {
    File monitoring_script
    File h5_input
    String input_vcf_name
    Int disk_size
    String docker
    Boolean no_address
  }
  command <<<
    bash ~{monitoring_script} > monitoring.log &

    source ~/.bashrc
    conda activate genomics.py3

    evaluate_concordance.py \
            --input_file ~{h5_input} \
            --output_prefix ~{input_vcf_name}.report \
            --use_for_group_testing variant_type

    >>>
  runtime {
    memory: "32 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
  }
  output {
    File eval_report_h5 = "~{input_vcf_name}.report.h5"
    File thresholds_report = "~{input_vcf_name}.report.thresholds.txt"
    File monitoring_log = "monitoring.log"
  }
}

task HardThresholdVCF {
  input { 
    File thresholds
    File input_vcf
    File input_vcf_index
    String output_basename
    String score_key
    Int disk_size
    String docker
    File monitoring_script
    Boolean no_address
  }

  command <<<
    set -eo pipefail
    bash ~{monitoring_script} > monitoring.log &    
    
    tail -n +2 ~{thresholds} > ~{thresholds}.body
    cat ~{thresholds}.body | awk 'BEGIN { FS=","} \
    {print "-filter \"VARIANT_TYPE == \047" $1 "\047 && ~{score_key} \
     < "$2 "\" --filter-name LOW_SCORE_"$1 " "}' | tr '\n' ' ' > tmpfile
    
    echo "Contents of tmpfile"
    cat tmpfile

    params=$( cat tmpfile )
    echo "Params"
    echo $params

    echo gatk --java-options "-Xmx20g"  VariantFiltration \
    -V ~{input_vcf} \
    -O ~{output_basename}.vcf.gz \
    "$params" >tmpfile
    . ./tmpfile
  >>>
  runtime {
    memory: "32 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"    
    File monitoring_log = "monitoring.log"
  }
}

task AnnotateSampleVCF {
    input {
        File input_vcf
        File input_vcf_index
        String output_basename
        Int disk_size = ceil(size(input_vcf, "GB") * 2) + 50
        String docker
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String flow_order = "TGCA"
    }

    command <<<
        gatk --java-options "-Xmx15g" \
            VariantAnnotator \
            -O ~{output_basename}.vcf.gz \
            -V ~{input_vcf} \
            -R ~{ref_fasta} \
            -G StandardFlowBasedAnnotation \
            --flow-order-for-annotations ~{flow_order}
    >>>
    runtime {
        memory: "16 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }

    output {
        File output_vcf_file = "~{output_basename}.vcf.gz"
        File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
    }
}




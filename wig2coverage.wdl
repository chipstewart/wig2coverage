task wig2coverage_task_1 {
    File interval_list
    File wig_coverage_file
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions
    String stub=basename(wig_coverage_file,".txt.gz")
    command {
        set -euo pipefail

        echo "${interval_list}"
        echo "${wig_coverage_file}"
        echo "${stub}"

        ls -lath
        
        python /opt/src/wig2cov.py -i ${interval_list} -w ${wig_coverage_file}

        ls -lath

    }

    output {
        File annotated_interval_file="${stub}.covered_bases.tsv"
        String total_targeted_covered_bases=read_string("${stub}.total_targeted_covered_bases.txt")
        String total_targeted_bases=read_string("${stub}.total_targeted_bases.txt")
        String total_non_targeted_covered_bases=read_string("${stub}.total_non_targeted_covered_bases.txt")
        File covered_bases_summary_file="${stub}.covered_bases_summary.txt"
    }

    runtime {
        docker : "docker.io/chipstewart/wig2coverage_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '3'}"
    }

    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }
}

workflow wig2coverage {

    File interval_list
    File wig_coverage_file

    call wig2coverage_task_1 {
        input: 
        interval_list=interval_list,
        wig_coverage_file=wig_coverage_file
    }

}
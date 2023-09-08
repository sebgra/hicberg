# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_2:
    input:
        genome = GENOME,
        r1 = R1,
        r2 = R2,
        forward_sorted_bam = temp("/home/sardine/Bureau/hic_test/1.sorted.bam"),
        reverse_sorted_bam = temp("/home/sardine/Bureau/hic_test/2.sorted.bam")
    
    output:
        forward_mapped_bam = "/home/sardine/Bureau/hic_test/group1.1.bam",


    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUTPUT} -r {RATE} \
        -t {THREADS} -m {MODE}  -e {ENZYMES_0} -e {ENZYMES_1} -s {SENSITIVITY} -n {BANK_NAME} -R {ROUNDS} \
        -M {MAGNITUDE}  --start-stage groups  --exit-stage build
        """
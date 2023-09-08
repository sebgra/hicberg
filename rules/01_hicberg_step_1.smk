# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_1:
    input:
        genome = GENOME,
        r1 = R1,
        r2 = R2,
        size_fragments = "/home/sardine/Bureau/hic_test/fragments_fixed_sizes.txt",

    params:
        main_directory = directory("/home/sardine/Bureau/hic_test/"),

    output: 
        forward_sorted_bam = temp("/home/sardine/Bureau/hic_test/1.sorted.bam"),
        reverse_sorted_bam = temp("/home/sardine/Bureau/hic_test/2.sorted.bam")

    shell: 
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUTPUT} -r {RATE} \
        -t {THREADS} -m {MODE}  -e {ENZYMES_0} -e {ENZYMES_1} -s {SENSITIVITY} -n {BANK_NAME} -R {ROUNDS} \
        -M {MAGNITUDE}   --start-stage bam  --exit-stage groups
        """
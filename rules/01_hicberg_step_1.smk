# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_1:
    input:
        genome = GENOME,
        r1 = R1,
        r2 = R2,
        for_sam = "/home/sardine/Bureau/hic_test/1.sam",
        rev_sam = "/home/sardine/Bureau/hic_test/2.sam",

    params:
        main_directory = directory("/home/sardine/Bureau/hic_test/"),

    output: 
        unmapped_for_bam = "/home/sardine/Bureau/hic_test/group0.1.bam",
        unmapped_rev_bam = "/home/sardine/Bureau/hic_test/group0.2.bam",
        mapped_for_bam = "/home/sardine/Bureau/hic_test/group1.1.bam",
        mapped_rev_bam = "/home/sardine/Bureau/hic_test/group1.2.bam",
        multi_mapped_for_bam = "/home/sardine/Bureau/hic_test/group2.1.bam",
        multi_mapped_rev_bam = "/home/sardine/Bureau/hic_test/group2.2.bam"

    shell: 
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUTPUT} -r {RATE} \
        -t {THREADS} -m {MODE}  -e {ENZYMES_0} -e {ENZYMES_1} -s {SENSITIVITY} -n {BANK_NAME} -R {ROUNDS} \
        -M {MAGNITUDE} -f --start-stage bam  --exit-stage stats
        """
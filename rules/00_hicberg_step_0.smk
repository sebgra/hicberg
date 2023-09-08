# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_0:
    input:
        genome = GENOME,
        r1 = R1,
        r2 = R2,
    output:
        index_folder = directory("/home/sardine/Bureau/hic_test/index/"),
        chrom_sizes = "/home/sardine/Bureau/hic_test/fragments_fixed_sizes.txt",


    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUTPUT} -r {RATE}  -t {THREADS} -m {MODE}  -e {ENZYMES_0} \
        -e {ENZYMES_1} -s {SENSITIVITY} -n {BANK_NAME} -R {ROUNDS} -M {MAGNITUDE} -f --start-stage fastq  --exit-stage bam
        """
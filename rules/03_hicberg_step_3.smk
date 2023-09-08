# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_3:
    input:
        genome = GENOME,
        r1 = R1,
        r2 = R2,
        mapped_for_bam = "/home/sardine/Bureau/hic_test/group1.1.bam",
        
        
    
    output:
        unrescued_matrix = "/home/sardine/Bureau/hic_test/unrescued_map.cool",

    
    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUTPUT} -r {RATE} \
        -t {THREADS} -m {MODE}  -e {ENZYMES_0} -e {ENZYMES_1} -s {SENSITIVITY} -n {BANK_NAME} -R {ROUNDS} \
        -M {MAGNITUDE} -f --start-stage build  --exit-stage stats
        """
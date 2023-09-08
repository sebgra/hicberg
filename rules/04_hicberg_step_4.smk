# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_4:
    input:
        genome = GENOME,
        r1 = R1,
        r2 = R2,
        unrescued_matrix = "/home/sardine/Bureau/hic_test/unrescued_map.cool",
        
        
    
    output:
        restriction_map = "/home/sardine/Bureau/hic_test/restriction_map.npy",

    
    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUTPUT} -r {RATE} \
        -t {THREADS} -m {MODE}  -e {ENZYMES_0} -e {ENZYMES_1} -s {SENSITIVITY} -n {BANK_NAME} -R {ROUNDS} \
        -M {MAGNITUDE} -f --start-stage stats  --exit-stage rescue
        """
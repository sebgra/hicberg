# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_4:
    input:
        genome = lambda w: join(REF_DIR, config[w.libraries]['ref']),
        r1 = lambda w: join(config[w.libraries]['R1']),
        r2 = lambda w: join(config[w.libraries]['R2']),
        unrescued_matrix = join(OUT_DIR, '{libraries}', "unrescued_map.cool"),
        
    params:
    #     main_directory = directory("/home/sardine/Bureau/hic_test/"),
        name = '{libraries}',
    
    output:
        # restriction_map = temp("/home/sardine/Bureau/hic_test/"),
        restriction_map = join(OUT_DIR, '{libraries}', "restriction_map.npy")

    
    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUT_DIR} -r {RATE} \
        -t {THREADS} -m {MODE}  -e {ENZYMES_0} -e {ENZYMES_1} -s {SENSITIVITY} -n {params.name} -R {ROUNDS} \
        -M {MAGNITUDE}  --start-stage stats  --exit-stage rescue
        """
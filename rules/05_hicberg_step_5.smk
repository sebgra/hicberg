# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_5:
    input:
        genome = lambda w: join(REF_DIR, config[w.libraries]['ref']),
        r1 = lambda w: join(config[w.libraries]['R1']),
        r2 = lambda w: join(config[w.libraries]['R2']),
        restriction_map = join(OUT_DIR, '{libraries}', "restriction_map.npy")
        
    params:
        name = '{libraries}',
        
    
    output:
        # rescued_matrix = temp("/home/sardine/Bureau/hic_test/contacts/matrices/rescued_map.cool")
        rescued_matrix = join(OUT_DIR, '{libraries}', "contacts", "matrices", "rescued_map.cool")
    
    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUT_DIR} -r {RATE} \
        -t {THREADS} -m {MODE}  -e {ENZYMES_0} -e {ENZYMES_1} -s {SENSITIVITY} -n {params.name} -R {ROUNDS} \
        -M {MAGNITUDE} --start-stage rescue  --exit-stage final
        """
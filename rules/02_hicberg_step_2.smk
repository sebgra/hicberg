# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_2:
    input:
        genome = lambda w: join(REF_DIR, config[w.libraries]['ref']),
        r1 = lambda w: join(config[w.libraries]['R1']),
        r2 = lambda w: join(config[w.libraries]['R2']),
        forward_sorted_bam = temp(join(OUT_DIR, '{libraries}', "1.sorted.bam")),
        reverse_sorted_bam = temp(join(OUT_DIR, '{libraries}', "2.sorted.bam"))
    params:
    #     main_directory = directory("/home/sardine/Bureau/hic_test/"),
        name = '{libraries}',
    
    output:
        forward_mapped_bam = join(OUT_DIR, '{libraries}', "group1.1.bam")


    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUT_DIR} -r {RATE} \
        -t {THREADS} -m {MODE}  -e {ENZYMES_0} -e {ENZYMES_1} -s {SENSITIVITY} -n {params.name} -R {ROUNDS} \
        -M {MAGNITUDE}  --start-stage groups  --exit-stage build
        """
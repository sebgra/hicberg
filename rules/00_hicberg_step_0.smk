# Set parameters.
shell.prefix("set -euo pipefail;")


rule hicberg_step_0:
    input:
        genome = lambda w: join(REF_DIR, config[w.libraries]['ref']),
        r1 = lambda w: join(config[w.libraries]['R1']),
        r2 = lambda w: join(config[w.libraries]['R2']),
        # name = species
    output:
        # index_folder = directory(join(OUT_DIR, '{libraries}', "index/")),
        chrom_sizes = join(OUT_DIR, '{libraries}', "fragments_fixed_sizes.txt"),

    params:
        name = '{libraries}',



    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUT_DIR} -r {RATE}  -t {THREADS} -m {MODE}  -e {ENZYMES_0} \
        -e {ENZYMES_1} -s {SENSITIVITY} -n {params.name} -R {ROUNDS} -M {MAGNITUDE} -f --start-stage fastq  --exit-stage bam
        """
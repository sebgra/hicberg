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
        sampling_rate = lambda w: samples.sampling_rates[w.libraries],
        enzyme = lambda w: samples.enzymes[w.libraries],
        mode = lambda w: samples.modes[w.libraries],
        rounds = lambda w: samples.rounds[w.libraries],
        magnitude = lambda w: samples.magnitudes[w.libraries],
        max_reports = lambda w: samples.max_reports[w.libraries],
        
    
    output:
        rescued_matrix = join(OUT_DIR, '{libraries}', "contacts", "matrices", "rescued_map.cool")
    
    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUT_DIR} -r {params.sampling_rate}  -t {THREADS} \
        -m {params.mode}  -e {params.enzyme} -s very-sensitive -n {params.name} -R {params.rounds} -M {params.magnitude} -k {params.max_reports} \
        --start-stage rescue  --exit-stage final -f
        """
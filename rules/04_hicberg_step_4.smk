# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_4:
    input:
        genome = lambda w: join(REF_DIR, config[w.libraries]['ref']),
        r1 = lambda w: join(config[w.libraries]['R1']),
        r2 = lambda w: join(config[w.libraries]['R2']),
        unrescued_matrix = join(OUT_DIR, '{libraries}', "unrescued_map.cool"),
        
    params:
        name = '{libraries}',
        sampling_rate = lambda w: samples.sampling_rates[w.libraries],
        enzyme = lambda w: samples.enzymes[w.libraries],
        mode = lambda w: samples.modes[w.libraries],
        kernel_size = lambda w: samples.kernel_sizes[w.libraries],
        deviation = lambda w: samples.deviations[w.libraries],
        max_reports = lambda w: samples.max_reports[w.libraries],
        circularities = lambda w: samples.circularity[w.libraries],
    
    output:
        restriction_map = temp(join(OUT_DIR, '{libraries}', "restriction_map.npy"))

    threads: 16

    
    shell:
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUT_DIR} -r {params.sampling_rate}  -t {threads} \
        -m {params.mode}  -e {params.enzyme} -s very-sensitive -n {params.name} -K {params.kernel_size} -d {params.deviation} -k {params.max_reports} \
        -c {params.circularities} --start-stage stats  --exit-stage rescue -f
        """
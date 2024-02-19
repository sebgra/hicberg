# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_3:
    input:
        genome = lambda w: join(REF_DIR, config[w.libraries]['ref']),
        r1 = lambda w: join(config[w.libraries]['R1']),
        r2 = lambda w: join(config[w.libraries]['R2']),
        forward_mapped_bam = join(OUT_DIR, '{libraries}', "group1.1.bam")
        
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
        forward_mapped_bam = temp(join(OUT_DIR, '{libraries}', "unrescued_map.cool")),

    threads: 16
    
    shell:
        """
        hicberg pipeline -o {OUT_DIR} -r {params.sampling_rate}  -t {threads} \
        -m {params.mode}  -e {params.enzyme} -s very-sensitive -n {params.name} -K {params.kernel_size} -d {params.deviation} -k {params.max_reports} \
        -c {params.circularities} --start-stage build  --exit-stage stats -f {input.genome} {input.r1} {input.r2}
        """
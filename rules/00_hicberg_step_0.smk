# Set parameters.
shell.prefix("set -euo pipefail;")


rule hicberg_step_0:
    input:
        genome = lambda w: join(REF_DIR, config[w.libraries]['ref']),
        r1 = lambda w: join(config[w.libraries]['R1']),
        r2 = lambda w: join(config[w.libraries]['R2']),
    output:
        chrom_sizes = join(OUT_DIR, '{libraries}', "fragments_fixed_sizes.txt"),
    log: 'hicberg_step_0-{libraries}.log'
    params: # extracted from samples.csv
        name = '{libraries}',
        sampling_rate = lambda w: samples.sampling_rates[w.libraries],
        enzyme = lambda w: samples.enzymes[w.libraries],
        mode = lambda w: samples.modes[w.libraries],
        kernel_size = lambda w: samples.kernel_sizes[w.libraries],
        deviation = lambda w: samples.deviations[w.libraries],
        max_reports = lambda w: samples.max_reports[w.libraries],
        circularities = lambda w: samples.circularity[w.libraries],

    threads: 2



    shell:
        """
        hicberg pipeline  -o {OUT_DIR} -r {params.sampling_rate}  -t {threads} \
        -m {params.mode}  -e {params.enzyme} -s very-sensitive -n {params.name} -K {params.kernel_size} -d {params.deviation} -k {params.max_reports} -f \
        -c {params.circularities} --start-stage fastq  --exit-stage bam -f {input.genome} {input.r1} {input.r2}
        """
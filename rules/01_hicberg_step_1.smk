# Set parameters.
shell.prefix("set -euo pipefail;")

rule hicberg_step_1:
    input:
        genome = lambda w: join(REF_DIR, config[w.libraries]['ref']),
        r1 = lambda w: join(config[w.libraries]['R1']),
        r2 = lambda w: join(config[w.libraries]['R2']),
        size_fragments = join(OUT_DIR, '{libraries}', "fragments_fixed_sizes.txt"),

    params:
    #     main_directory = directory("/home/sardine/Bureau/hic_test/"),
        name = '{libraries}',
        sampling_rate = lambda w: samples.sampling_rates[w.libraries],
        enzyme = lambda w: samples.enzymes[w.libraries],
        mode = lambda w: samples.modes[w.libraries],
        rounds = lambda w: samples.rounds[w.libraries],
        magnitude = lambda w: samples.magnitudes[w.libraries],
        max_reports = lambda w: samples.max_reports[w.libraries],
        circularities = lambda w: samples.circularity[w.libraries],


    output: 
        # forward_sorted_bam = temp("/home/sardine/Bureau/hic_test/1.sorted.bam"),
        # reverse_sorted_bam = temp("/home/sardine/Bureau/hic_test/2.sorted.bam"),
        forward_sorted_bam = temp(join(OUT_DIR, '{libraries}', "1.sorted.bam")),
        reverse_sorted_bam = temp(join(OUT_DIR, '{libraries}', "2.sorted.bam"))

    threads: 16

    shell: 
        """
        hicberg pipeline -g {input.genome} --fq-for {input.r1} --fq-rev {input.r2} -o {OUT_DIR} -r {params.sampling_rate} -t {threads} \
        -m {params.mode}  -e {params.enzyme} -s very-sensitive -n {params.name} -R {params.rounds} -M {params.magnitude} -k {params.max_reports} \
        -c {params.circularities} --start-stage bam  --exit-stage groups -f
        """
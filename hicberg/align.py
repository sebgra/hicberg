import sys
import shutil
from os import getcwd
from os.path import join
from pathlib import Path
import subprocess as sp
import uuid
from hicberg import logger
from hicberg.utils import _run_command


def hic_build_index(
    genome: str,
    output_dir: str = None,
    cpus: int = 1,
    verbose: bool = False,
    aligner: str = "bowtie2",
) -> None:
    """
    Building of bowtie2 index (.bt2l files) for read alignment.

    Parameters
    ----------
    genome : str
        Path to the genome file along which reads are going to be aligned.
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the Bowtie2 index files should be stored, by default None
    verbose : bool, optional
        Set wether or not the shell command should be printed, by default False
    aligner : str, optional
        Aligner algorithm to use for alignment, by default bowtie2
        Supported aligners: 'bowtie2', 'bwa', 'minimap2'.
    """

    logger.info("Start building index for alignment")

    supported_aligners = {
        "bowtie2": "bowtie2-build",
        "bwa": "bwa index",  # bwa index is used for building the index
        "minimap2": "minimap2",  # minimap2 does not have a separate index build command like bowtie2 or bwa
    }

    if aligner not in supported_aligners:
        raise ValueError(
            f"Aligner '{aligner}' is not supported. Supported aligners are: {', '.join(supported_aligners.keys())}"
        )

    aligner_command = supported_aligners[aligner].split()[0]  # Get the base command

    print(f"aligner_command: {aligner_command}")

    # Map aligner names to their primary executable names
    aligner_executables = {
        "bowtie2": "bowtie2-build",  # For index building
        "bwa": "bwa",  # For index building (via 'bwa index')
        "minimap2": "minimap2",  # No separate index building step
    }

    if aligner not in aligner_executables:
        raise ValueError(
            f"Aligner '{aligner}' is not supported. Supported aligners are: {', '.join(aligner_executables.keys())}"
        )

    base_command = aligner_executables[aligner]

    # --- Robustly check for aligner executable existence in PATH ---
    if shutil.which(base_command) is None:
        raise RuntimeError(
            f"Aligner '{base_command}' not found in your system's PATH. "
            f"Please ensure {aligner} is installed and its executable is accessible.\n"
            f"You can often install it with: conda install {aligner}"
        )
    logger.info(f"Aligner '{aligner}' ({base_command}) found in PATH.")

    genome_path = Path(genome)

    if not genome_path.is_file():

        raise ValueError(f"Genome file {genome} not found")

    if output_dir is None:
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(
            f"Output path {output_path} does not exist. Please provide existing ouput path."
        )

    sample = Path(genome).stem
    index_path = Path(output_dir, sample)

    match aligner:

        case "bowtie2":

            cmd_index = f"bowtie2-build -q -f --threads {cpus} --large-index {genome} {index_path}"
            _run_command([cmd_index], description="Indexing genome using Bowtie2")

            return index_path

        case "bwa":

            cmd_index = f"bwa index -p {join(index_path)} {genome}"

            _run_command([cmd_index], description="Indexing genome using BWA")

            return index_path

        case "minimap2":

            logger.info(
                f"Minimap2 does not require a separate index build step like Bowtie2 or BWA."
            )
            logger.info(
                f"When using minimap2, the index is often built on-the-fly during alignment."
            )

        case "_":  # Default case

            pass
        
# --- Helper functions for aligner-specific parameters ---

def _get_bowtie2_params(sensitivity: str, max_alignment: int) -> list[str]:
    """Get Bowtie2 parameters based on sensitivity and max_alignment."""
    params = [f"--{sensitivity}"]
    if max_alignment is not None and max_alignment != -1:
        params.extend(["-k", str(max_alignment)])
    else:
        params.append("-a") # -a means report all alignments for Bowtie2.
    return params

def _get_bwa_params(
    bwa_subcommand: str,
    sensitivity: str,
    max_alignment: int = None,
    cpus: int = None # cpus is used by 'mem' and 'aln', but not 'samse' directly
) -> list[str]:
    """
    Get BWA parameters for a specific subcommand based on sensitivity and other settings.

    Parameters
    ----------
    bwa_subcommand : str
        The specific BWA subcommand ('mem', 'aln', 'samse').
    sensitivity : str
        Sensitivity setting ('very-fast', 'fast', 'sensitive', 'very-sensitive').
        Note: The interpretation of 'sensitivity' varies significantly between subcommands.
    max_alignment : Optional[int]
        Maximum number of alignments to be returned. Used by 'mem' (conceptually via -a) and 'samse' (-n).
    cpus : Optional[int]
        Number of threads. Used by 'mem' (-t) and 'aln' (-t).

    Returns
    -------
    list[str]
        A list of BWA command-line parameters for the specified subcommand.
    """
    params = []

    match bwa_subcommand:
        case "mem":
            if cpus is not None:
                params.extend(["-t", str(cpus)])

            # BWA-MEM specific sensitivity parameters
            if sensitivity == "very-fast":
                params.extend(["-k", "25", "-B", "6", "-O", "8,8", "-E", "2,2", "-T", "50"])
            elif sensitivity == "fast":
                params.extend(["-k", "22", "-B", "5", "-O", "7,7", "-E", "1,1", "-T", "40"])
            elif sensitivity == "sensitive":
                # BWA-MEM default parameters are often considered 'sensitive'
                pass
            elif sensitivity == "very-sensitive":
                params.extend(["-k", "15", "-B", "2", "-O", "4,4", "-E", "0,0", "-T", "10"])
            else:
                logger.warning(f"Unknown BWA-MEM sensitivity '{sensitivity}'. Using default parameters.")
            
            # For Hi-C, '-a' (output all found alignments) is common for bwa mem.
            params.append("-a")

        case "aln":
            if cpus is not None:
                params.extend(["-t", str(cpus)])

            # BWA-ALN specific sensitivity parameters
            if sensitivity == "very-fast":
                params.extend(["-n", "10", "-o", "2", "-q", "10"]) # Looser alignment parameters
            elif sensitivity == "fast":
                params.extend(["-n", "5", "-o", "1", "-q", "5"])
            elif sensitivity == "sensitive":
                params.extend(["-n", "2", "-o", "1", "-q", "0"]) # Common for Hi-C 'aln'
            elif sensitivity == "very-sensitive":
                params.extend(["-n", "1", "-o", "0", "-q", "0"]) # Very stringent
            else:
                logger.warning(f"Unknown BWA-ALN sensitivity '{sensitivity}'. Using default parameters (-n 2 -o 1 -q 0).")
                params.extend(["-n", "2", "-o", "1", "-q", "0"]) # Fallback to 'sensitive'

        case "samse":
            # BWA-SAMSE specific parameters
            # sensitivity parameter is largely irrelevant for samse, as it processes aln output
            if max_alignment is not None and max_alignment != -1:
                params.extend(["-n", str(max_alignment)])
            # No other common parameters needed for 'sensitivity' or 'cpus' for samse itself.

        case _:
            raise ValueError(f"Unsupported BWA subcommand: '{bwa_subcommand}'. Must be 'mem', 'aln', or 'samse'.")

    return params

def _get_minimap2_params(sensitivity: str, max_alignment: int) -> list[str]:
    """Get Minimap2 parameters based on sensitivity and max_alignment."""
    params = []
    # Mapping generic sensitivity to Minimap2's presets or manual tweaks
    if sensitivity == "very-fast":
        params.extend(["-x", "sr", "-k", "28", "-w", "20"])
    elif sensitivity == "fast":
        params.extend(["-x", "sr", "-k", "24", "-w", "15"])
    elif sensitivity == "sensitive":
        params.extend(["-x", "sr"])
    elif sensitivity == "very-sensitive":
        params.extend(["-x", "sr", "-k", "15", "-w", "5", "-B", "2", "-O", "2,10", "-E", "0,0", "-m", "20", "--max-chain-skip", "50"])
    else:
        logger.warning(f"Unknown Minimap2 sensitivity '{sensitivity}'. Using '-x sr' preset.")
        params.extend(["-x", "sr"])

    if max_alignment is not None and max_alignment != -1:
        params.extend(["-N", str(max_alignment)])
        params.append("--secondary=yes")
    else:
        params.append("--secondary=yes") # Output all secondary alignments by default

    return params


def hic_align(
    index: str,
    genome: str,
    fq_for: str,
    fq_rev: str,
    sensitivity: str = "very-sensitive",
    max_alignment: int = None,
    cpus: int = 1,
    output_dir: str = None,
    verbose: bool = False,
    aligner: str = "bowtie2",
    read_type: str = "short"
) -> None:
    """
    Alignment of reads from HiC experiments along an indexed genome.

    Parameters
    ----------
    index : str
        Path to the index of the genome along which reads are going to be aligned (path to .bt2l files). Default to None, index files are searched to sample_name/data/index/sample_name.
    genome : str
        Path to the genome to perform alignment on.
    fq_for : str
        Path to .fasta containing set of reads to align (forward mate).
    fq_rev : str
        Path to .fasta containing set of reads to align (forward mate).
    sensitivity : str, optional
        Sensitivity of the alignment., by default 'very_sensitive'
    max_alignment : int, optional
        Maximum number of alignments to be returned, by default None
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the alignment files (.sam) should be stored, by default None
    verbose : bool, optional
        Set wether or not the shell command should be printed, by default False
    aligner : str, optional
        Aligner algorithm to use for alignment, by default bowtie2
        Supported aligners: 'bowtie2', 'bwa', 'minimap2'.
    """

    logger.info("Start aligning reads")

    fq_for_path, fq_rev_path = Path(fq_for), Path(fq_rev)

    if not fq_for_path.is_file() or not fq_rev_path.is_file():

        raise IOError(
            f"Wrong path to fastq files : {fq_for_path} or {fq_rev_path} given. \
                    Pease provide existing files."
        )

    if output_dir is None:
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(
            f"Output path {output_path} does not exist. Please provide existing output path."
        )

    index_path = Path(output_path / index) if index is not None else None

    match aligner:

        case "bowtie2":
            
            bt2_params = _get_bowtie2_params(sensitivity=sensitivity, max_alignment=max_alignment)
   
            cmd_alignment_for = [
            "bowtie2",
            *bt2_params,
            "-p", str(cpus),
            "-x", str(index_path),
            "-S", str(output_path / "1.sam"),
            str(fq_for_path),
        ]
            
            cmd_alignment_rev = [
            "bowtie2",
            *bt2_params,
            "-p", str(cpus),
            "-x", str(index_path),
            "-S", str(output_path / "2.sam"),
            str(fq_rev_path),
        ]

            _run_command(
                cmd_alignment_for,
                description="Forward reads alignment with Bowtie2",
                verbose=True,
            )
            _run_command(
                cmd_alignment_rev,
                description="Reverse reads alignment with Bowtie2",
                verbose=True,
            )

        case "bwa":
            if max_alignment is None or max_alignment == -1:
                logger.info("Using BWA-MEM workflow.")
                # Call _get_bwa_params for 'mem'
                bwa_mem_params = _get_bwa_params("mem", sensitivity, cpus=cpus)
                
                # Construct cmd_alignment_for and cmd_alignment_rev using bwa_mem_params
                # ... (rest of BWA-MEM block as before, removing '-t' if added inside _get_bwa_params)
                cmd_alignment_for = [
                    "bwa", "mem",
                    *bwa_mem_params, # Now includes '-t' from _get_bwa_params
                    str(genome),
                    str(fq_for_path),
                    "-o", str(output_path / "1.sam")
                ]
                
                cmd_alignment_rev = [
                    "bwa", "mem",
                    *bwa_mem_params,
                    str(genome),
                    str(fq_rev_path),
                    "-o", str(output_path / "2.sam")
                ]
                

                _run_command(
                    cmd_alignment_for,
                    description="Forward reads alignment with BWA",
                    verbose=True,
                )
                _run_command(
                    cmd_alignment_rev,
                    description="Reverse reads alignment with BWA",
                    verbose=True,
                )

            elif max_alignment is not None:

                logger.info("Using BWA-MEM workflow.")
                # Call _get_bwa_params for 'mem'
                bwa_aln_params = _get_bwa_params("aln", sensitivity, cpus=cpus)
                bwa_samse_params = _get_bwa_params("samse", sensitivity, cpus=cpus)
                
                # Construct cmd_alignment_for and cmd_alignment_rev using bwa_mem_params
                # ... (rest of BWA-MEM block as before, removing '-t' if added inside _get_bwa_params)
                cmd_pre_alignment_for = [
                    "bwa", "aln",
                    *bwa_aln_params, # Now includes '-t' from _get_bwa_params
                    "-f", str(output_path / "1.sai"),
                    str(genome),
                    str(fq_for_path),
                ]
                
                
                cmd_alignment_for = [
                    "bwa", "samse",
                    *bwa_samse_params, # Now includes '-t' from _get_bwa_params
                    "-f", str(output_path / "1.sam"),
                    str(genome),
                    str(output_path / "1.sai"),
                    str(fq_for_path),
                ]
                
                
                cmd_pre_alignment_rev = [
                    "bwa", "aln",
                    *bwa_aln_params,
                    "-f", str(output_path / "2.sai"),
                    str(genome),
                    str(fq_rev_path),
                ]
                
                cmd_alignment_rev = [
                    "bwa", "samse",
                    *bwa_samse_params, # Now includes '-t' from _get_bwa_params
                    "-f", str(output_path / "2.sam"),
                    str(genome),
                    str(output_path / "2.sai"),
                    str(fq_for_path),
                ]
                

                _run_command(
                    cmd_pre_alignment_for,
                    description="Forward reads alignment with BWA",
                    verbose=True,
                )
                
                _run_command(
                    cmd_alignment_for,
                    description="Forward reads alignment with BWA",
                    verbose=True,
                )
                
                _run_command(
                    cmd_pre_alignment_rev,
                    description="Reverse reads alignment with BWA",
                    verbose=True,
                )
                
                _run_command(
                    cmd_alignment_rev,
                    description="Reverse reads alignment with BWA",
                    verbose=True,
                )

        case "minimap2":

            if max_alignment is None or max_alignment == -1:

                cmd_alignment_rev = f"minimap2 -a --secondary yes -t {cpus} -o {output_path / '1.sam'} {genome} {fq_for}"
                cmd_alignment_for = f"minimap2 -a --secondary yes -t {cpus} -o {output_path / '2.sam'} {genome} {fq_rev}"

            elif max_alignment is not None:

                cmd_alignment_for = f"minimap2 -a --secondary yes -N {max_alignment} -t {cpus} -o {output_path / '1.sam'} {genome} {fq_for}"
                cmd_alignment_rev = f"minimap2 -a --secondary yes -N {max_alignment} -t {cpus} -o {output_path / '2.sam'} {genome} {fq_rev}"

            _run_command(
                [cmd_alignment_for],
                description="Forward reads alignment with Minimap2",
                verbose=True,
            )
            _run_command(
                [cmd_alignment_rev],
                description="Reverse reads alignment with Minimap2",
                verbose=True,
            )
        # if verbose:

        #     logger.info(cmd_alignment_for)
        #     logger.info(cmd_alignment_rev)

        # p_for = sp.Popen(
        #     [cmd_alignment_for], shell=True, stdout=sp.PIPE, stderr=sp.PIPE
        # )
        # stdout_for, stderr_for = p_for.communicate()
        # p_rev = sp.Popen(
        #     [cmd_alignment_rev], shell=True, stdout=sp.PIPE, stderr=sp.PIPE
        # )
        # stdout_rev, stderr_rev = p_rev.communicate()

        # if stdout_for:
        #     logger.info(stdout_for.decode("ascii"))
        # if stderr_for:
        #     logger.info(stderr_for.decode("ascii"))

        # if stdout_rev:
        #     logger.info(stdout_rev.decode("ascii"))
        # if stderr_rev:
        #     logger.info(stderr_rev.decode("ascii"))

        case "_":

            pass

    logger.info(f"Alignment saved at {output_path}")


def hic_view(
    sam_for: str = "1.sam",
    sam_rev: str = "2.sam",
    cpus: int = 1,
    output_dir: str = None,
    verbose: bool = False,
) -> None:
    """
    Conversion of .sam alignment files to .bam alignment format (using samtools).

    Parameters
    ----------
    sam_for : str, optional
        Path to forward .sam alignment file, by default "1.sam"
    sam_rev : str, optional
        Path to reverse .sam alignment file, by default "2.sam"
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the alignment files (.bam) should be stored, by default None
    verbose : bool, optional
        Set wether or not the shell command should be printed, by default False
    """

    logger.info("Start converting .sam to .bam")

    try:

        sp.check_output(["samtools", "--help"])

    except OSError:

        raise RuntimeError(
            "Samtools not found; check if it is installed and in $PATH\n install Samtools with : conda install samtools"
        )

    if output_dir is None:
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(
            f"Output path {output_path} does not exist. Please provide existing output path."
        )

    cmd_view_for = f"samtools view -h  -b {output_path / sam_for} -o {output_path / '1.bam'} --threads {cpus}"
    cmd_view_rev = f"samtools view -h  -b {output_path / sam_rev} -o {output_path / '2.bam'} --threads {cpus}"

    _run_command(
        [cmd_view_for], description=".sam to .bam conversion of forward alignments"
    )
    _run_command(
        [cmd_view_rev], description=".sam to .bam conversion of reverse alignments"
    )

    # Delete .sam files after .bam conversion
    (output_path / sam_for).unlink()
    (output_path / sam_rev).unlink()

    logger.info(f"Compressed  alignment done at {output_path}")


def hic_sort(
    bam_for: str = "1.bam",
    bam_rev: str = "2.bam",
    cpus: int = 1,
    output_dir: str = None,
    verbose: bool = False,
) -> None:
    """
    Sort .bam alignment files by read_name  (using samtools).

    Parameters
    ----------
    bam_for : str, optional
        Forward alignment file to be sorted, by default "1.bam"
    bam_rev : str, optional
        Reverse alignment file to be sorted, by default "2.bam"
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the alignment files (.bam) should be stored, by default None
    verbose : bool, optional
        Set wether or not the shell command should be printed, by default False

    """
    logger.info("Start sorting .bam alignment files")

    try:

        sp.check_output(["samtools", "--help"])

    except OSError:

        raise RuntimeError(
            "Samtools not found; check if it is installed and in $PATH\n install Samtools with : conda install samtools"
        )

    if output_dir is None:
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(
            f"Output path {output_path} does not exist. Please provide existing ouput path."
        )

    id_for = uuid.uuid4()
    id_rev = uuid.uuid4()

    cmd_sort_for = f"samtools sort -n -T {id_for} {output_path / '1.bam'} -o {output_path / '1.sorted.bam'} --threads {cpus}"
    cmd_sort_rev = f"samtools sort -n -T {id_rev} {output_path / '2.bam'} -o {output_path / '2.sorted.bam'} --threads {cpus}"

    _run_command(
        [cmd_sort_for],
        description="Sorting of forward alignments considering read names",
    )
    _run_command(
        [cmd_sort_rev],
        description="Sorting of reverse alignments considering read names",
    )

    (output_path / "1.bam").unlink()
    (output_path / "2.bam").unlink()

    logger.info(f"Sorted alignment done at {output_path}")


def hic_index(
    bam_for: str = "1.sorted.bam",
    bam_rev: str = "2.sorted.bam",
    cpus: int = 1,
    output_dir: str = None,
    verbose: bool = False,
) -> None:
    """
    Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast random access.

    Parameters
    ----------
    bam_for : str, optional
        Forward alignment file to be indexed, by default "1.sorted.bam"
    bam_rev : str, optional
        Reverse alignment file to be indexed,, by default "2.sorted.bam"
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the alignment files (.bam) should be stored, by default None
    verbose : bool, optional
        Set wether or not the shell command should be printed, by default False

    """

    try:

        sp.check_output(["samtools", "--help"])

    except OSError:

        raise RuntimeError(
            "Samtools not found; check if it is installed and in $PATH\n install Samtools with : conda install samtools"
        )

    if output_dir is None:
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(
            f"Output path {output_path} does not exist. Please provide existing output path."
        )

    cmd_index_for = f"samtools index -b {bam_for} -@ {cpus}"
    cmd_index_rev = f"samtools index -b {bam_rev} -@ {cpus}"

    if verbose:

        logger.info(cmd_index_for)
        logger.info(cmd_index_rev)

    sp.run([cmd_index_for], shell=True)
    sp.run([cmd_index_rev], shell=True)

    logger.info(f"Indexed alignment done at {output_path}")

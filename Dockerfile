FROM  condaforge/mambaforge

LABEL Name=hicberg Version=0.0.1

RUN apt-get update && apt-get install -y gcc
RUN mamba config --add channels bioconda
RUN mamba install -c conda-forge -y pip \
  && mamba clean -afy

RUN apt update && \
    apt install -y sudo && \
    addgroup --gid 1000 nonroot && \
    adduser --uid 1000 --gid 1000 --disabled-password --gecos "" nonroot && \
    echo 'nonroot ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers

RUN mamba install bioconda::bowtie2
RUN mamba install bioconda::samtools


COPY --chown=nonroot:nonroot * ./ /app/
WORKDIR /app
RUN chmod 777 /app/*


# RUN pip install -Ur requirements.txt
RUN pip install -e .

USER nonroot


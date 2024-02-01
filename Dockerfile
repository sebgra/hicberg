FROM  condaforge/mambaforge

LABEL Name=hicberg Version=0.0.1

COPY * ./ /app/
WORKDIR /app

# RUN apt-get update && apt-get install -y wget bzip2 
# RUN wget -qO-  https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba 
# RUN touch /root/.bashrc 
# RUN ./bin/micromamba shell init -s bash -p /opt/conda  
# RUN grep -v '[ -z "\$PS1" ] && return' /root/.bashrc  > /opt/conda/bashrc   # this line has been modified 
# RUN apt-get clean autoremove --yes 
# RUN rm -rf /var/lib/{apt,dpkg,cache,log}

# SHELL ["bash", "-l" ,"-c"]

RUN apt-get update && apt-get install -y gcc
RUN mamba config --add channels bioconda
RUN mamba install -c conda-forge -y pip \
  && mamba clean -afy

RUN mamba install bioconda::bowtie2
RUN mamba install bioconda::samtools
# RUN pip install -Ur requirements.txt
RUN pip install -e .

ENTRYPOINT [ "hicberg" ]
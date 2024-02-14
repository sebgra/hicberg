# Get image
FROM  condaforge/mambaforge

# Add args for building
ARG UID
ARG GID

# Metadata
LABEL Name=hicberg Version=0.0.1

# Update base components and set-up mamba
RUN apt-get update && apt-get install -y gcc
RUN mamba config --add channels bioconda
RUN mamba install -c conda-forge -y pip \
  && mamba clean -afy

# Add non root user
RUN apt update && \
    apt install -y sudo && \
    addgroup --gid $GID nonroot && \
    adduser --uid $UID --gid $GID --disabled-password --gecos "" nonroot && \
    echo 'nonroot ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers

#Install bowtie2 and samtools
RUN mamba install bioconda::bowtie2
RUN mamba install bioconda::samtools

# Copy necessary files
COPY --chown=nonroot:nonroot * ./ /app/

# Set working directory
WORKDIR /app
RUN chmod 777 /app/*


# Install HiC-BERG
RUN pip install -e .

# Switch to non root user
USER nonroot


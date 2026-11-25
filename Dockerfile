# mobilink with the programs that `mobilink run` calls (DIAMOND, minimap2, Kraken2).
# Build:  docker build -t mobilink .
# Use:    docker run --rm -v "$PWD":/data mobilink mobilink db --out /data/mobilink_db
#         docker run --rm -v "$PWD":/data mobilink mobilink run /data/reads.fastq.gz --db /data/mobilink_db --out /data/results
# The databases and the Kraken2 index are not in the image: they are large and have their own licences.
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/mobilink
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 numpy scipy pandas click python-edlib "diamond>=2.1" minimap2 kraken2 pip \
    && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache-dir --no-deps /opt/mobilink && mobilink --version && diamond version && minimap2 --version
WORKDIR /data
CMD ["mobilink", "--help"]

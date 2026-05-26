FROM mambaorg/micromamba:1.5.8

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        fastp=0.23.4 \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH

CMD ["fastp", "--version"]

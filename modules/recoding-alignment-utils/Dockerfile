FROM mambaorg/micromamba:1.5.8

USER root

# Bio CLI tools + Python; pip-installable libs come next.
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        python=3.12 \
        samtools=1.21 \
        bcftools=1.21 \
        pip \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH

# Python deps for the in-house scripts in scripts/.
# Kept narrow on purpose — extend if a new script imports something else.
RUN pip install --no-cache-dir \
        pysam==0.23.3 \
        anndata>=0.10 \
        pandas>=2.2,<3.0 \
        numpy \
        matplotlib \
        seaborn \
        biopython

CMD ["python", "--version"]

FROM mambaorg/micromamba:1.5.8

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        bedtools=2.31.1 \
        ucsc-bedgraphtobigwig=445 \
    && micromamba clean --all --yes
 
ENV PATH=/opt/conda/bin:$PATH

CMD ["bedGraphToBigWig"]



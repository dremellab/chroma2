FROM continuumio/miniconda3:24.5.0-0
RUN conda install -y -c bioconda -c conda-forge subread=2.0.6 \
    && conda clean -ay
ENV PATH="/opt/conda/bin:$PATH"

FROM condaforge/miniforge3:latest

ENV CONDA_DEFAULT_ENV=psps

RUN set -ex \
    && conda config --add channels bioconda \
    && conda config --set channel_priority stric \
    && conda create --name ${CONDA_DEFAULT_ENV} \
        python=3.8 bcftools tabix bgzip --yes \
    && echo "conda activate ${CONDA_DEFAULT_ENV}" >> ~/.bashrc

SHELL ["conda", "run", "-n", "${CONDA_DEFAULT_ENV}", "/bin/bash", "-c"]
RUN set -ex \
    && conda install --channels bioconda --yes \
        cyvcf2 gffutils liftover numpy pathlib2 \
        pandarallel pandas pysam pybedtools tqdm \
    && conda clean --all --yes

CMD [ "/bin/bash" ]
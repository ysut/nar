FROM ubuntu:16.04 AS builder

ARG MANTA_VERSION=1.6.0
ARG MANTA_DIR=/opt/manta-${MANTA_VERSION}

WORKDIR /opt

RUN set -ex \
    && apt-get update -qq \
    && apt-get install -yqq --no-install-recommends \
        bzip2 gcc g++ make python wget zlib1g-dev

RUN set -ex \
    && mkdir -p ${MANTA_DIR} \
    && wget -c --no-check-certificate \
        https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2 \
    && tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2 \
    && ./manta-${MANTA_VERSION}.release_src/configure --jobs=6 --prefix=${MANTA_DIR} \
    && make -j6 install


# === 2nd stage ===
FROM ubuntu:16.04 AS manta

ENV PATH=/opt/manta-1.6.0/bin:$PATH

RUN set -ex \
    && apt-get update -qq \
    && apt-get install -yqq --no-install-recommends \
        bzip2 \
        python \
        zlib1g-dev

COPY --from=builder /opt/manta-1.6.0 /opt/manta-1.6.0

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu \
    && chsh -s /bin/bash \
    && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu

# RUN mkdir -p $OPT
# COPY --from=builder $OPT $OPT
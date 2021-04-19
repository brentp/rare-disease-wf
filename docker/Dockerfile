FROM ubuntu:18.04
MAINTAINER brentp

ARG slivar_version=v0.2.1
ARG htslib_version=1.12
ARG glnexus_version=v1.3.1

RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     curl wget ca-certificates git-core less netbase \
     openjdk-8-jre-headless unzip \
     g++ cmake autoconf make file valgrind \
     libjemalloc-dev libzip-dev libsnappy-dev libbz2-dev zlib1g-dev liblzma-dev libzstd-dev \
     python3-pyvcf bcftools pv

ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.1

RUN  git clone --depth 1 --branch $glnexus_version https://github.com/dnanexus-rnd/GLnexus && \
     cd GLnexus && \
     cmake -DCMAKE_BUILD_TYPE=Release . && make -j8 &&  \
     mv /GLnexus/glnexus_cli /usr/local/bin && \
     rm -rf /GLnexus

RUN \
    git clone --depth 1 https://github.com/ebiggers/libdeflate.git && \
    cd libdeflate && make -j 2 CFLAGS='-fPIC -O3' libdeflate.a && \
    cp libdeflate.a /usr/local/lib && cp libdeflate.h /usr/local/include && \
    cd .. && rm -rf libdeflate && \
    git clone --recursive --depth 1 --branch $htslib_version https://github.com/samtools/htslib && \
    cd htslib && autoheader && autoconf && \
    ./configure --with-libdeflate && \
    cd .. && make -j4 CFLAGS="-fPIC -O3" -C htslib install && \
    git clone --depth 1 --branch $htslib_version https://github.com/samtools/samtools && \
    cd samtools && autoheader && autoconf && \
    ./configure --without-curses --enable-s3 --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && \
    git clone --depth 1 --branch $htslib_version https://github.com/samtools/bcftools && \
    cd bcftools && autoheader && autoconf && \
    ./configure --without-curses --enable-s3 --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && cd ../ && \
    rm -rf bcftools /slivar /htslib /samtools /nim /tmp/*

RUN cd / && \
    wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip -q snpEff_latest_core.zip && \
    rm -rf snpEff_latest_core.zip snpEff/examples/ snpEff/galaxy snpEff/scripts && \
    echo '#!/bin/bash' > /usr/local/bin/snpEff && \
    echo 'java -jar -Xmx4g /snpEff/snpEff.jar "$@"' >> /usr/local/bin/snpEff && \
    chmod +x /usr/local/bin/snpEff && \
    snpEff -version # && snpEff download GRCh38.99

ADD https://github.com/brentp/slivar/releases/download/$slivar_version/slivar /usr/local/bin
ADD https://raw.githubusercontent.com/brentp/slivar/$slivar_version/js/slivar-functions.js /opt/slivar/
RUN chmod +x /usr/local/bin/slivar

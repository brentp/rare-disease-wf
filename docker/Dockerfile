FROM gcr.io/deepvariant-docker/deepvariant:1.1.0
WORKDIR /dv/
RUN cp -r /opt/deepvariant /dv
RUN cp -r /opt/models /dv
ENV PATH="/dv/deepvariant/bin/:$PATH"

FROM brentp/slivar:latest
WORKDIR /slivar
RUN cp /opt/slivar/slivar-functions.js /slivar
RUN cp /usr/bin/slivar /slivar
ENV PATH="/slivar/:$PATH"

COPY CentOS-Base.repo /etc/yum.repos.d/

RUN yum update -y && \
    yum upgrade -y && \
    yum install -y xz git libzip-devel curl wget zlib-devel xz-devel bzip2-devel openssl-devel libcurl-devel

RUN wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo && \
    yum install -y devtoolset-2-gcc

RUN yum install -y devtoolset-2-binutils devtoolset-2-gcc-c++ && \
    source scl_source enable devtoolset-2 && \
    echo "source scl_source enable devtoolset-2" >> ~/.bashrc && \
    echo "source scl_source enable devtoolset-2" >> ~/.bash_profile && \
    wget --quiet https://ftp.gnu.org/gnu/m4/m4-1.4.18.tar.gz  &&\
    source scl_source enable devtoolset-2 && \
    wget --quiet https://ftp.gnu.org/gnu/m4/m4-1.4.18.tar.gz && \
    tar xzf m4-1.4.18.tar.gz && cd m4* && ./configure && make && make install && cd .. && rm -rf m4* && \
    wget --quiet http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz && \
    tar xzf autoconf-2.69.tar.gz && \
    cd autoconf* && ./configure && make && make install && cd .. && rm -rf autoconf* && \
    git clone --depth 1 https://github.com/ebiggers/libdeflate.git && \
    cd libdeflate && make -j 2 CFLAGS='-fPIC -O3' libdeflate.a && \
    cp libdeflate.a /usr/local/lib && cp libdeflate.h /usr/local/include && \
    cd .. && rm -rf libdeflate && \
    git clone --recursive https://github.com/samtools/htslib && \
    cd htslib && git checkout 1.12 && autoheader && autoconf && \
    ./configure --enable-libcurl --with-libdeflate && \
    cd .. && make -j4 CFLAGS="-fPIC -O3" -C htslib install && \
    echo "/usr/local/lib" >> /etc/ld.so.conf && \
    ldconfig

RUN source scl_source enable devtoolset-2 && \
    git clone https://github.com/samtools/bcftools && \
    cd bcftools && git checkout 1.12 && autoheader && autoconf && \
    ./configure --enable-s3 --enable-libcurl --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && \
    git clone https://github.com/samtools/samtools && \
    cd samtools && git checkout 1.12 && autoheader && autoconf && \
    ./configure --without-curses --enable-s3 --enable-libcurl --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && rm -rf bcftools /slivar/htslib samtools /nim /tmp/*

RUN yum install -y sudo yum install java-1.8.0-openjdk-headless unzip && \
    yum clean all && \
    rm -rf /var/cache/yum

RUN cd /slivar && \
    wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip -q snpEff_latest_core.zip && \
    rm -rf snpEff_latest_core.zip snpEff/examples/ snpEff/galaxy snpEff/scripts && \
    echo '#!/bin/bash' > /usr/local/bin/snpEff && \
    echo 'java -jar -Xmx4g /slivar/snpEff/snpEff.jar "$@"' >> /usr/local/bin/snpEff && \
    chmod +x /usr/local/bin/snpEff && \
    snpEff -version
    #&& \ snpEff download GRCh38.99

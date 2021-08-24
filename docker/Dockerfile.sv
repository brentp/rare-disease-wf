FROM ubuntu:16.04

ARG tiwih_version=v0.1.6
ARG duphold_version=v0.2.3
ARG slivar_version=v0.2.5
ARG dysgu_version=v1.2.7

ADD https://github.com/brentp/slivar/releases/download/$slivar_version/slivar /usr/local/bin
ADD https://raw.githubusercontent.com/brentp/slivar/$slivar_version/js/slivar-functions.js /opt/slivar/
ADD https://repo.anaconda.com/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh /opt
ADD https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 /opt
ADD https://github.com/brentp/tiwih/releases/download/$tiwih_version/tiwih /usr/local/bin
ADD https://github.com/brentp/duphold/releases/download/$duphold_version/duphold /usr/local/bin

RUN apt-get update -qq \
    && apt-get -qy install --no-install-recommends build-essential bzip2 python2.7 git ca-certificates parallel \
    && chmod a+x /usr/local/bin/tiwih \
    && chmod a+x /usr/local/bin/duphold \
    && chmod a+x /usr/local/bin/slivar

WORKDIR /opt
RUN  cp /usr/bin/python2.7 /usr/bin/python2 \
	     && tar xjf manta-1.6.0.centos6_x86_64.tar.bz2 \
	     && rm manta-1.6.0.centos6_x86_64.tar.bz2 \
	     && rm -rf manta-1.6.0.centos6_x86_64/manta-1.6.0.centos6_x86_64/share
ADD configManta.py.ini /opt/manta-1.6.0.centos6_x86_64/bin/

ENV PATH=/opt/manta-1.6.0.centos6_x86_64/bin:/opt/manta-1.6.0.centos6_x86_64/libexec/:/opt/:manta-1.6.0.centos6_x86_64/lib/python/:$PATH

## conda + jasminesv
WORKDIR /opt
ENV PATH=/opt/miniconda/bin:$PATH
RUN sh Miniconda3-py37_4.9.2-Linux-x86_64.sh -b -p /opt/miniconda/ && \
	    rm Miniconda3-py37_4.9.2-Linux-x86_64.sh && \
      conda config --add channels bioconda && \
      conda config --add channels conda-forge && \
      conda update --freeze-installed -n base -yc defaults conda && \
      conda install --freeze-installed -yc anaconda nomkl gxx_linux-64 gcc_linux-64 autoconf make && \
      conda install --freeze-installed -c conda-forge nomkl curl && \
      conda init bash && \
      conda install --freeze-installed -yc bioconda nomkl jasminesv">=1.1.2" samtools">=1.10" pysam bcftools paragraph">=2.3" htslib">=1.10" && \
      git clone -b $dysgu_version --recursive https://github.com/kcleal/dysgu.git && \
      true

## dysgu && svpack
SHELL ["/bin/bash", "-c"]
RUN \
    cd dysgu/dysgu/htslib && \
    git checkout 1.11 && \
    conda init bash && \
    . /opt/miniconda/etc/profile.d/conda.sh && \
    eval "$(conda shell.bash hook)" && \
    conda activate root && \
    autoheader && autoconf && ./configure && \
    make -j 4 install && cd ../.. && \
    conda install --freeze-installed -y --file requirements.txt  && \
    for f in dysgu/*.pyx; do cython --cplus $f; done && \
    python setup.py install && \
    dysgu --version && \
    git clone https://github.com/amwenger/svpack/ &&  \
    cp ./svpack/svpack /usr/local/bin && \
    chmod +x /usr/local/bin/svpack && \
    svpack -h

RUN for f in $(find   / -type f -name vcfgraph.py); do \
	    curl -SsLo $f https://raw.githubusercontent.com/brentp/paragraph/bp-dev/src/python/lib/grm/vcfgraph/vcfgraph.py;  \
	    chmod a+x $f;   \
	    d=$(dirname $(dirname $f)); \
	    np="$d/vcf2paragraph/__init__.py" ; \
            curl -SsLo $np https://raw.githubusercontent.com/brentp/paragraph/bp-dev/src/python/lib/grm/vcf2paragraph/__init__.py; \
    done


RUN \
    conda clean -afy && \
    find /opt/miniconda/ -follow -type f -name '*.a' -delete && \
		find / -type f -name "*.a" -delete

FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY conda.yml .
RUN conda env update -n root -f conda.yml && conda clean -a

RUN conda install -c bioconda samtools=1.9
RUN conda install -c bioconda picard=2.20.6
RUN conda install -c bioconda mosdepth=0.2.6
RUN conda install -c bioconda fastp=0.20.0
RUN conda install -c conda-forge fd-find
RUN conda install -c bioconda sambamba
RUN conda install -c bioconda star=2.7.9a
RUN conda install -c bioconda spades
RUN conda install -c bioconda blast

RUN apt-get --allow-releaseinfo-change update && \
   apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*
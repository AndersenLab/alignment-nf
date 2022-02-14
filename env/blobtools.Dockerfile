FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

# COPY conda.yml .
# RUN conda env update -n root -f conda.yml && conda clean -a

RUN conda install -c bioconda star=2.7.9a
RUN conda install -c bioconda spades
RUN conda install -c bioconda blobtools
RUN conda install -c bioconda blast

RUN apt-get --allow-releaseinfo-change update && \
   apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*
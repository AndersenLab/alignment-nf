FROM continuumio/miniconda3
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY conda.yml .
RUN conda env update -n root -f conda.yml && conda clean -a

# RUN conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
# RUN conda install -c bioconda pysam --update-deps
# RUN conda install -c bioconda star=2.7.9a
# RUN conda install -c bioconda spades 
# # spades=3.15.3
# RUN conda install -c bioconda blast=2.12
RUN conda install -c bioconda blobtools


# RUN cd blobtools/

RUN apt-get --allow-releaseinfo-change update && \
   apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*
FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY conda.yml .
RUN conda env update -n root -f conda.yml && conda clean -a

# RUN git clone https://github.com/DRL/blobtools.git

RUN conda install -c anaconda python=3.8.12
RUN conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
# RUN conda install -c bioconda pysam --update-deps
RUN conda install -c bioconda star=2.7.9a
RUN conda install -c bioconda spades 
# spades=3.15.3
RUN conda install -c bioconda blast=2.12
# RUN conda install -c bioconda blobtools

# get blobtools
RUN wget https://github.com/DRL/blobtools/archive/refs/tags/blobtools_v1.1.1.zip
RUN unzip blobtools_v1.1.1.zip
RUN export PATH=$PWD/blobtools-blobtools_v1.1.1/:$PATH

# RUN cd blobtools/

RUN apt-get --allow-releaseinfo-change update && \
   apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*
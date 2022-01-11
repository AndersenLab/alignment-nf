# old align docker
# FROM continuumio/miniconda3:4.7.12
# RUN apt-get update && apt-get install -y procps && apt-get clean -y
# RUN conda create -n env bioconda::bwa bioconda::samtools bioconda::bcftools bioconda::picard=2.22.0 wget && \
#     conda clean -a
# RUN wget https://github.com/biod/sambamba/releases/download/v0.7.1/sambamba-0.7.1-linux-static.gz
# RUN gzip -d sambamba-0.7.1-linux-static.gz && \
#     chmod +x sambamba-0.7.1-linux-static && \
#     mv sambamba-0.7.1-linux-static /opt/conda/envs/env/bin/sambamba
# RUN echo "source activate env" > ~/.bashrc
# ENV PATH /opt/conda/envs/env/bin:$PATH

# new alignment docker
FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

RUN conda install bioconda::bwa=0.7.17
RUN conda install bioconda::sambamba
RUN conda install bioconda::samtools=1.9
# RUN conda install bamtools=2.5.1
RUN conda install bioconda::picard=2.20.6
# RUN conda install bcftools=1.9
# RUN conda install fastqc=0.11.8
# RUN conda install telseq=0.0.2
RUN conda install bioconda::multiqc
RUN conda install r::r=3.6.0
RUN conda install bioconda::mosdepth=0.2.6
RUN conda install bioconda::star=2.7.9a
RUN conda install bioconda::spades
RUN conda install bioconda::blobtools
RUN conda install bioconda::blast

RUN Rscript -e "install.packages('tidyverse', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('plotly', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  

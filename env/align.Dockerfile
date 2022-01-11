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

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

RUN Rscript -e "install.packages('tidyverse', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('plotly', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"

RUN apt-get --allow-releaseinfo-change update && \
	apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*

ARG TARGETPLATFORM=linux/amd64
FROM --platform=${TARGETPLATFORM} ubuntu:22.04

# Download and install dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    perl \
    r-cran-devtools \
    wget \
    unzip \
    python3 \
    r-base \
    build-essential \
    zlib1g-dev \
    locales \
    ghostscript \
    parallel \
    bcftools

# Configure locale
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Download and install PLINK 1.9 and 2
## NOTE: The PLINK2 AWS resource URL is subject to change; it looks like they remove old versions
##       once a year or so. Check the bucket URL(s) if either link doesn't work: 
##       https://s3.amazonaws.com/plink1-assets or https://s3.amazonaws.com/plink2-assets

RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip && \
    unzip plink_linux_x86_64_20230116.zip && \
    mv plink /bin/plink1.9 && \
    rm plink_linux_x86_64_20230116.zip

RUN wget -O plink2_data.zip https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20250609.zip && \
    unzip plink2_data.zip && \
    mv plink2 /bin/plink2 && \
    rm plink2_data.zip

WORKDIR /usr/src

COPY . .

# symbolic link to 'old' plink 
RUN ln -s /bin/plink1.9 /bin/plink

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh
    
# Add conda to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Create env with compatible python version
RUN conda create -n sim python=3.11 -y

# Make sure the compadre environment is activated by default
#SHELL ["conda", "run", "-n", "sim", "/bin/bash", "-c"]

#################################################

ARG REFERENCE_VERSION=v0.1.0

# Install reference VCF files
RUN apt-get update && apt-get install -y curl jq && \
    curl -s https://api.github.com/repos/belowlab/pedigree-simulator/releases/tags/${REFERENCE_VERSION} | \
    jq -r '.assets[] | select(.name | endswith(".vcf.gz")) | .browser_download_url' | \
    xargs -I {} wget {} && \
    mkdir -p /usr/src/data/reference/AMR /usr/src/data/reference/EUR && \
    mv *AMR*.vcf.gz /usr/src/data/reference/AMR/ && \
    mv *EUR*.vcf.gz /usr/src/data/reference/EUR/ && \
    apt-get remove -y jq && apt-get autoremove -y

RUN wget https://github.com/belowlab/pedigree-simulator/releases/download/${REFERENCE_VERSION}/all.frq.zip && \
    unzip all.frq.zip && \
    rm all.frq.zip && \
    mv all.frq /usr/src/data/all.frq

RUN wget https://github.com/belowlab/pedigree-simulator/releases/download/${REFERENCE_VERSION}/map.gz.zip && \
    unzip map.gz.zip && \
    rm map.gz.zip && \
    mv map.gz /usr/src/data/map.gz

RUN wget https://github.com/belowlab/pedigree-simulator/releases/download/${REFERENCE_VERSION}/extract_mega.bim.zip && \
    unzip extract_mega.bim.zip && \
    rm extract_mega.bim.zip && \
    mv extract_mega.bim /usr/src/data/extract_mega.bim

RUN wget https://github.com/belowlab/pedigree-simulator/releases/download/${REFERENCE_VERSION}/hapmap3.zip && \
    unzip hapmap3.zip && \
    rm hapmap3.zip && \
    mv hapmap3 /usr/src/dependencies/PRIMUS_v1.9.0/lib/hapmap3

RUN wget https://github.com/belowlab/pedigree-simulator/releases/download/${REFERENCE_VERSION}/KDE_data.zip && \
    unzip KDE_data.zip && \
    rm KDE_data.zip && \
    mv KDE_data /usr/src/dependencies/PRIMUS_v1.9.0/lib/KDE_data

SHELL ["conda", "run", "-n", "sim", "/bin/bash", "-c"]

##################################################

# clean up empty folders
RUN rm -r __MACOSX

# Install the KernSmooth R package
RUN Rscript -e "install.packages('KernSmooth', repos='http://cran.rstudio.com/')"

# Add perl path to where primus is expecting it 
RUN mkdir -p /usr/src/perl
ENV PERL_PATH=/usr/src/perl
ENV PERL5LIB=""
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:/usr/src/lib/perl_modules/:/usr/src/lib/perl_modules/PRIMUS/:$PERL5LIB

WORKDIR /usr/src/src
# This is slightly redundant, but singularity doesn't work otherwise
ENTRYPOINT ["perl", "/usr/src/src/main.pl"] 
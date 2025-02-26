SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools 
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(shell uname)-$(shell uname -m).sh" && bash Miniforge3-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Miniforge3-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y --override-channels -c conda-forge -c bioconda gseapy matplotlib numpy pandas r-azimuth bioconductor-biomart bioconductor-clusterprofiler bioconductor-complexheatmap r-cowplot r-data.table bioconductor-deseq2 bioconductor-dittoseq r-dplyr bioconductor-enhancedvolcano r-forcats r-ggplot2 r-ggpubr r-ggrepel r-ggsignif r-hrbrthemes bioconductor-org.hs.eg.db r-pals r-patchwork r-pheatmap r-rcolorbrewer r-readr r-readxl r-seurat r-stringi r-stringr r-tidyverse r-viridis && touch .tools

clean:
	rm -rf mamba/ .mamba .tools

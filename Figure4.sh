#!/bin/bash

## Example #1. Overlay, intron shrinkage, gene annotation, PDF output, custom size and colors
# -b: tsv file stating BAMS to use for the plotting
# -c: chromosome coordinates
# -g: GTF annotation file
# -M: minimum read coverage for drawing junctions
# -C: Index of column with color levels (1-based)
# -O: Index of column with overlay levels (1-based)
# -A: value to plot on top of junctions
# --shrink: Shrink the junctions by a factor for nicer display 
# --alpha: Transparency level for density histogram 
# --base-size: Base font size of the plot in pch
# --ann-height: Height of annotation plot in inches
# --height: Height of the individual signal plot in inches
# --widht: Width of the plot in inches
# -P: color palette file
# -R: image resolution

### package github page:  https://github.com/guigolab/ggsashimi

# Figure 4g: sashimi plot of RPL27A

region='11:8682588-8686000' 

gtf="/g/korbel/zafferani/MicroExonator/integrated_patients/Report/out.high_quality_calibrated.gtf"
gene='RPL27A'
#file in which bam files of patients are stored
input=input_bams_${gene}.tsv

python ggsashimi.py -b $input -c $region -F 'png' -M 3 -R 600 -A mean -g $gtf \
-C 3 -O 3 --alpha 1 --base-size=0.01 --height=3 --ann-height 5 --width=12 -P palette.txt -o sashimi_4g


##################################################################################################


# Figure 4h: sashimi plot of FOS

region='14:75278826-75282230'

gene='FOS'
input=input_bams_${gene}.tsv

python ggsashimi.py -b $input -c $region -F 'png' -M 3 -R 600 -A mean -g $gtf \
-C 3 -O 3 --alpha 1 --base-size=0.01 --height=3 --ann-height 5 --width=12 -P palette.txt -o sashimi_4h

#!/usr/bin/bash
OUTDIR="output"
mkdir -p $OUTDIR

# Example 01
echo "Run pyMSAviz CLI example 01..."
pymsaviz -i HIGD2A.fa -o ${OUTDIR}/cli_example01.png --dpi 100

# Example 02
echo "Run pyMSAviz CLI example 02..."
pymsaviz -i MRGPRG.fa -o ${OUTDIR}/cli_example02.png --wrap_length 80 --dpi 100 \
         --color_scheme Taylor --show_consensus --show_count

# Example 03
echo "Run pyMSAviz CLI example 03..."
pymsaviz -i MRGPRG.fa -o ${OUTDIR}/cli_example03.png --start 100 --end 160 --dpi 100 \
         --color_scheme Flower --show_grid --show_consensus --consensus_color tomato

echo -e "\nFinished all example CLI run."

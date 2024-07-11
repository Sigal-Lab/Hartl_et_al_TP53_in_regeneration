#!/bin/bash

ROI=/sc-scratch/sc-scratch-cc13-multiplex_ihc-pipeline/HB/GATK/mm10/call_regions_mm10_without_MaSat.bed


# controlfreec uses pathToSamtools_ + " view -@ " + SambambaThreads_ + " " +mateFileName;
# so $1 = view, $2=-@, $3 = threads, $4=filename

samtools $1 $2 $3 -L $ROI $4



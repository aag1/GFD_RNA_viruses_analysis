#!/bin/bash


sed '1d; 2d; $d' ICTV_GFD_picobirna_RdRp_msa_outgr.FigTreeMidpoint.nexus |
sed "s/\[&label=//g" |
sed "s/\]//g" |
sed "s/'//g" |
cut -d' ' -f5 > ICTV_GFD_picobirna_RdRp_msa_outgr.FigTreeMidpoint.nwk

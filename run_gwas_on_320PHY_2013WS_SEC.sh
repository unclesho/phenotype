#!/bin/bash

./gwas_pipeline.sh am1 320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned
./gwas_pipeline.sh am2 320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned
./gwas_pipeline.sh app_am 320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned
./gwas_pipeline.sh app_ap 320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned
./gwas_pipeline.sh mc_ap 320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned
./gwas_pipeline.sh sc_ap 320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned


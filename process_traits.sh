#!/bin/bash

emmax -v -d 10 -t ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned -p am1/am1.txt -k ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned.hIBS.kinf -o am1/am1_emmax | tee am1/am1_emmax.log
emmax -v -d 10 -t ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned -p am2/am2.txt -k ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned.hIBS.kinf -o am2/am2_emmax | tee am2/am2_emmax.log
emmax -v -d 10 -t ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned -p app_am/app_am.txt -k ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned.hIBS.kinf -o app_am/app_am_emmax | tee app_am/app_am_emmax.log
emmax -v -d 10 -t ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned -p app_ap/app_ap.txt -k ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned.hIBS.kinf -o app_ap/app_ap_emmax | tee app_ap/app_ap_emmax.log
emmax -v -d 10 -t ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned -p mc_ap/mc_ap.txt -k ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned.hIBS.kinf -o mc_ap/mc_ap_emmax | tee mc_ap/mc_ap_emmax.log
emmax -v -d 10 -t ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned -p sc_ap/sc_ap.txt -k ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned.hIBS.kinf -o sc_ap/sc_ap_emmax | tee sc_ap/sc_ap_emmax.log
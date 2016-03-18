#!/bin/bash

sed '/^$/d' am1/am1_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps
sed '/^$/d' am2/am2_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' > am2/am2_emmax.ps.qqman.emmax_200kb.clumped.snps
sed '/^$/d' app_am/app_am_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' > app_am/app_am_emmax.ps.qqman.emmax_200kb.clumped.snps
sed '/^$/d' app_ap/app_ap_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' > app_ap/app_ap_emmax.ps.qqman.emmax_200kb.clumped.snps
sed '/^$/d' mc_ap/mc_ap_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' > mc_ap/mc_ap_emmax.ps.qqman.emmax_200kb.clumped.snps
sed '/^$/d' sc_ap/sc_ap_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' > sc_ap/sc_ap_emmax.ps.qqman.emmax_200kb.clumped.snps

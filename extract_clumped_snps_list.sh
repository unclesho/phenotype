#!/bin/bash

sed '/^$/d' am1/am1_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps_list
sed '/^$/d' am2/am2_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > am2/am2_emmax.ps.qqman.emmax_200kb.clumped.snps_list
sed '/^$/d' app_am/app_am_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > app_am/app_am_emmax.ps.qqman.emmax_200kb.clumped.snps_list
sed '/^$/d' app_ap/app_ap_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > app_ap/app_ap_emmax.ps.qqman.emmax_200kb.clumped.snps_list
sed '/^$/d' mc_ap/mc_ap_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > mc_ap/mc_ap_emmax.ps.qqman.emmax_200kb.clumped.snps_list
sed '/^$/d' sc_ap/sc_ap_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > sc_ap/sc_ap_emmax.ps.qqman.emmax_200kb.clumped.snps_list

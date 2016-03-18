#!/bin/bash

{ printf 'SNP\tCHR\tBP\tBETA\tP\n'; awk '{split($1,a,"_"); print $1 "\t" a[2]*1 "\t" a[3] "\t" $2 "\t" $3}' am1/am1_emmax.ps; } > am1/am1_emmax.ps.qqman
{ printf 'SNP\tCHR\tBP\tBETA\tP\n'; awk '{split($1,a,"_"); print $1 "\t" a[2]*1 "\t" a[3] "\t" $2 "\t" $3}' am2/am2_emmax.ps; } > am2/am2_emmax.ps.qqman
{ printf 'SNP\tCHR\tBP\tBETA\tP\n'; awk '{split($1,a,"_"); print $1 "\t" a[2]*1 "\t" a[3] "\t" $2 "\t" $3}' app_am/app_am_emmax.ps; } > app_am/app_am_emmax.ps.qqman
{ printf 'SNP\tCHR\tBP\tBETA\tP\n'; awk '{split($1,a,"_"); print $1 "\t" a[2]*1 "\t" a[3] "\t" $2 "\t" $3}' app_ap/app_ap_emmax.ps; } > app_ap/app_ap_emmax.ps.qqman
{ printf 'SNP\tCHR\tBP\tBETA\tP\n'; awk '{split($1,a,"_"); print $1 "\t" a[2]*1 "\t" a[3] "\t" $2 "\t" $3}' mc_ap/mc_ap_emmax.ps; } > mc_ap/mc_ap_emmax.ps.qqman
{ printf 'SNP\tCHR\tBP\tBETA\tP\n'; awk '{split($1,a,"_"); print $1 "\t" a[2]*1 "\t" a[3] "\t" $2 "\t" $3}' sc_ap/sc_ap_emmax.ps; } > sc_ap/sc_ap_emmax.ps.qqman

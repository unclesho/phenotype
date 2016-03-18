#!/bin/bash

# Roslen Anacleto
# 20160316

# Prepare the emmax files
plink2 --bfile ../genotype/$1 --recode12 --output-missing-genotype 0 --transpose --out ../genotype/$1
emmax-kin -v -h -s -d 10 ../genotype/$1

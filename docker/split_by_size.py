#!/usr/bin/env python
import sys
import re
import argparse


p = argparse.ArgumentParser()

p.add_argument("--chunk_size", type=int, default=30000000)
p.add_argument("fai", help="fasta index matching fasta used to align and call variants")
p.add_argument("-e", help="chromosome exclude pattern", default="^chrEBV$|^GL|^NC|_random$|Un_|^HLA\-|_alt$|hap\d$")

a = p.parse_args()

excl = re.compile(a.e)
chunk_size = a.chunk_size

for line in open(a.fai):
    toks = line.split("\t")
    chrom_len = int(toks[1])
    if excl.search(toks[0]):
        continue

    for chrom_start in range(0, chrom_len + 1, chunk_size):
        chrom_end = min(chrom_len, chrom_start + chunk_size)
        print(f"{toks[0]}:{chrom_start+1}-{chrom_end}")

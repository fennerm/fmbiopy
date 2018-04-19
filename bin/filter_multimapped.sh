#!/usr/bin/env bash
## Filter multimapped reads from a bam file.
samtools view -h "$1" | grep -P "(NH:i:1|^@)" | samtools view -Sb -

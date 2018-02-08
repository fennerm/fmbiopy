#!/usr/bin/env bash
## Intersect two files
## Usage: Intersect.sh FILE1 FILE2
sort "$1" "$2" | uniq -d

#!/usr/bin/awk -f
## Determine the maximum column number in a tsv file with a variable number of
## columns.
BEGIN {
    max=0
}
{
    if (NF > max) {
        max=NF
    }
}
END {
    print max
}

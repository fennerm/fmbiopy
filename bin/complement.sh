#!/usr/bin/env bash
## Return the items in FILE1 that are not in FILE2
comm -23 <(sort "$1") <(sort "$2")

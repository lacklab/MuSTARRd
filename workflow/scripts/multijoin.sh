#!/bin/sh

# multijoin - join multiple files

echo -en "BC\t"
for arg in "$@"; do
    bname=${arg##*/}
    echo -en ${bname%.count.tsv}"\t"
done
echo -en "\n"

# https://stackoverflow.com/a/17649149

join_rec() {
    if [ $# -eq 1 ]; then
        join - "$1"
    else
        f=$1; shift
        join - "$f" | join_rec "$@"
    fi
}

if [ $# -le 2 ]; then
    join "$@" | tr " " "\t"
else
    f1=$1; f2=$2; shift 2
    join "$f1" "$f2" | join_rec "$@" | tr " " "\t"
fi

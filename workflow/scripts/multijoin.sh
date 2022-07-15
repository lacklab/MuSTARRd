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
        join -t $'\t' -a1 -a2 -e0 -o auto - "$1"
    else
        f=$1; shift
        join -t $'\t' -a1 -a2 -e0 -o auto - "$f" | join_rec "$@"
    fi
}

if [ $# -le 2 ]; then
    join -t $'\t' -a1 -a2 -e0 -o auto "$@"
else
    f1=$1; f2=$2; shift 2
    join -t $'\t' -a1 -a2 -e0 -o auto "$f1" "$f2" | join_rec "$@"
fi

{
    where1 = match($10, fwd)
    where2 = match($10, rev)
    if (where1 != 0 && where2 != 0 && where1 - 1 == umi1 && (length($10) - where2 - RLENGTH + 1) == umi2 && (is + umi1 + umi2) == length($10))
        print $0
}

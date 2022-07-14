{
    if ($9 >= is) {
        where1 = match($10, fwd)
        if (where1 != 0 && where1 - 1 == umi1) {
            print $0
        }
    }
    else if ($9 <= -is) {
        where2 = match($10, rev)
        if (where2 != 0 && (length($10) - where2 - RLENGTH + 1) == umi2) {
            print $0
        }
    }
}

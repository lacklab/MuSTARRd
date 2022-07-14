{
    nm = $5;
    md = $6;
    sub(/NM:i:/, "", nm);
    sub(/MD:Z:/, "", md);
    
    if (nm == 0) {
         printf "%s\t%d\t%d\t%s", $1, $7, $8, "WT\n";
    }
    else {
        split(md, a, "[ACTGN]", seps);

        delete a[length(a)];        

        printf "%s\t%d\t%d\t", $1, $7, $8;
        for (i in a) {
            pos = 0
            for (j = 1; j <= i; j++) {
                pos += a[j]+1
            }
            printf "%s:%d:%s>%s ", $3, pos, seps[i], substr($4, pos, 1);
        }
        printf "\n";
    }
}  

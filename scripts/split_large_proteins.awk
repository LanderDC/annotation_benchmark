/^>/ {
    if (seq) {
        L = length(seq);
        F = int((L + 2499) / 2500);
        R = int(L / F);
        for (i = 0; i < F; i++) {
            start = i * R + 1;    
            end = (i == F - 1) ? L : start + R - 1; 
            fragment = substr(seq, start, end - start + 1);
            print header "_split" (i + 1) "\n" fragment; 
        }
        seq = "";  
    }
    header = $0; 
    next;
}
{
    seq = seq $0;   
}
END {
    if (seq) {
        L = length(seq);
        F = int((L + 2499) / 2500);
        R = int(L / F);
        for (i = 0; i < F; i++) {
            start = i * R + 1;    
            end = (i == F - 1) ? L : start + R - 1; 
            fragment = substr(seq, start, end - start + 1);
            print header "_split" (i + 1) "\n" fragment; 
        }
    }
}
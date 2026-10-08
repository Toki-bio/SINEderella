# Plurality consensus of an aligned FASTA, degapped: a column gives its most frequent base (A, C, G, T)
# when that base is held by at least `plur` sequences, otherwise nothing. Ties go to the base that
# comes first in A, C, G, T. EMBOSS cons took more than 17 minutes on 900 sequences x 3,100 columns;
# this takes seconds. usage: awk -v plur=225 -v name=ref -f colcons.awk aligned.fa
function flush(   i, c) {
    if (seq == "") return
    seq = toupper(seq)
    for (i = 1; i <= length(seq); i++) {
        c = substr(seq, i, 1)
        if (c == "A") a[i]++; else if (c == "C") cc[i]++; else if (c == "G") g[i]++; else if (c == "T") t[i]++
    }
    if (length(seq) > L) L = length(seq)
    seq = ""
}
/^>/ { flush(); next }
{ gsub(/[ \t\r]/, ""); seq = seq $0 }
END {
    flush()
    out = ""
    for (i = 1; i <= L; i++) {
        best = a[i] + 0; b = "A"
        if (cc[i] + 0 > best) { best = cc[i] + 0; b = "C" }
        if (g[i] + 0 > best) { best = g[i] + 0; b = "G" }
        if (t[i] + 0 > best) { best = t[i] + 0; b = "T" }
        if (best >= plur) out = out b
    }
    print ">" name
    print out
}

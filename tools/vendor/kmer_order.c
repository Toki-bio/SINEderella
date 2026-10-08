/*
 * kmer_order: k-mer guide-tree leaf order of a FASTA file (SubFam 1.2.0 ordering step).
 *
 * A C port of the Python/numpy script that SubFam.sh used to embed (KMER_ORDER_PY), itself a port
 * of kmer-tree.js (ViewAlign, Toki-bio/MSA-viewer, MIT). Same arithmetic, same tie-breaking, so
 * the output is byte-identical to the Python version (tools/vendor/test_kmer_order.sh checks it).
 *
 *   kmer_order in.fasta out.fasta K CANONICAL(0/1) ORIENT(0/1)
 *
 * Weighted Jaccard distance on k-mer counts (1 - sum(min) / sum(max)), float32 matrix, UPGMA
 * (average linkage) with kmer-tree.js tie-breaking, at every merge the two clusters are flipped
 * so that their closest ends meet. ORIENT=1: each sequence is reverse-complemented if that shares
 * more k-mers with its predecessor in the order; such sequences get MAFFT's "_R_" id prefix.
 * Letters other than A C G T U (either case) are skipped inside k-mers, as in kmer-tree.js.
 *
 * Build: cc -O2 -ffp-contract=off [-fopenmp] -o kmer_order kmer_order.c -lm
 * (-ffp-contract=off: no fused multiply-add, so the float arithmetic matches numpy's.)
 * Memory: n*n*4 bytes for the distance matrix, as in the Python version.
 * Not reproduced: Python's str.strip() also strips a few exotic Unicode spaces; here only ASCII
 * whitespace is stripped from ids and sequence lines. CR, CRLF and LF all end a line, as in Python.
 */
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static void die(const char *msg) { fprintf(stderr, "kmer_order: %s\n", msg); exit(1); }
static void *xmalloc(size_t n) { void *p = malloc(n ? n : 1); if (!p) die("out of memory"); return p; }
static void *xcalloc(size_t n, size_t s) { void *p = calloc(n ? n : 1, s ? s : 1); if (!p) die("out of memory"); return p; }

/* ---------- FASTA ---------- */
typedef struct { char *name; char *seq; size_t len; } Rec;

static Rec *read_fasta(const char *path, size_t *nout) {
    FILE *fh = fopen(path, "rb");
    if (!fh) die("cannot open input");
    size_t cap = 1 << 20, len = 0;
    char *buf = xmalloc(cap + 1);
    for (;;) {
        if (len == cap) { cap *= 2; buf = realloc(buf, cap + 1); if (!buf) die("out of memory"); }
        size_t r = fread(buf + len, 1, cap - len, fh);
        if (r == 0) break;
        len += r;
    }
    fclose(fh);
    buf[len] = 0;
    /* universal newlines: CRLF and lone CR become LF */
    size_t w = 0;
    for (size_t r = 0; r < len; r++) {
        if (buf[r] == '\r') { buf[w++] = '\n'; if (r + 1 < len && buf[r + 1] == '\n') r++; }
        else buf[w++] = buf[r];
    }
    len = w;
    buf[len] = 0;

    size_t nrec = 0, reccap = 1024;
    Rec *recs = xmalloc(reccap * sizeof(Rec));
    size_t seqcap = 0;
    int have = 0;
    size_t pos = 0;
    while (pos < len) {
        size_t e = pos;
        while (e < len && buf[e] != '\n') e++;
        size_t a = pos, b = e;
        while (a < b && isspace((unsigned char)buf[a])) a++;
        while (b > a && isspace((unsigned char)buf[b - 1])) b--;
        if (e > pos && buf[pos] == '>') {          /* header: the line itself starts with '>' (as startswith) */
            size_t na = pos + 1, nb = e;
            while (na < nb && isspace((unsigned char)buf[na])) na++;
            while (nb > na && isspace((unsigned char)buf[nb - 1])) nb--;
            if (nrec == reccap) { reccap *= 2; recs = realloc(recs, reccap * sizeof(Rec)); if (!recs) die("out of memory"); }
            recs[nrec].name = xmalloc(nb - na + 1);
            memcpy(recs[nrec].name, buf + na, nb - na);
            recs[nrec].name[nb - na] = 0;
            recs[nrec].seq = xmalloc(1);
            recs[nrec].seq[0] = 0;
            recs[nrec].len = 0;
            seqcap = 1;
            nrec++;
            have = 1;
        } else if (have) {
            Rec *r = &recs[nrec - 1];
            size_t add = b - a;
            if (r->len + add + 1 > seqcap) {
                seqcap = (r->len + add + 1) * 2;
                r->seq = realloc(r->seq, seqcap);
                if (!r->seq) die("out of memory");
            }
            memcpy(r->seq + r->len, buf + a, add);
            r->len += add;
            r->seq[r->len] = 0;
        }
        pos = e + 1;
    }
    free(buf);
    *nout = nrec;
    return recs;
}

/* ---------- k-mers ---------- */
static signed char CODE[256];
static unsigned char COMP[256];

static void init_tables(void) {
    memset(CODE, -1, sizeof CODE);
    const char *l = "ACGTUacgtu";
    const int v[10] = {0, 1, 2, 3, 3, 0, 1, 2, 3, 3};
    for (int i = 0; i < 10; i++) CODE[(unsigned char)l[i]] = (signed char)v[i];
    for (int i = 0; i < 256; i++) COMP[i] = (unsigned char)i;
    const char *from = "ACGTUacgtuRYKMBVDHrykmbvdh", *to = "TGCAAtgcaaYRMKVBHDyrmkvbhd";
    for (int i = 0; from[i]; i++) COMP[(unsigned char)from[i]] = (unsigned char)to[i];
}

/* k-mer codes of seq (valid letters only); writes into out (capacity >= len), returns the count */
static size_t kmer_codes(const char *seq, size_t len, int k, int canonical, uint32_t *out) {
    uint32_t mask = (k == 16) ? 0xFFFFFFFFu : ((1u << (2 * k)) - 1);
    uint32_t h = 0, r = 0;
    int have = 0;
    size_t m = 0;
    for (size_t i = 0; i < len; i++) {
        int c = CODE[(unsigned char)seq[i]];
        if (c < 0) continue;
        h = ((h << 2) | (uint32_t)c) & mask;
        r = (r >> 2) | ((uint32_t)(3 - c) << (2 * (k - 1)));
        if (++have >= k) out[m++] = (canonical && r < h) ? r : h;
    }
    return m;
}

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}

/* profile with kmers seen in >= 2 sequences only (the others never add to sum(min)) */
typedef struct { uint32_t *code; uint32_t *cnt; size_t m; } Prof;

static double shared_min(const Prof *a, const Prof *b) {
    size_t i = 0, j = 0;
    double s = 0;
    while (i < a->m && j < b->m) {
        uint32_t x = a->code[i], y = b->code[j];
        if (x < y) i++;
        else if (x > y) j++;
        else { s += (double)(a->cnt[i] < b->cnt[j] ? a->cnt[i] : b->cnt[j]); i++; j++; }
    }
    return s;
}

static size_t N;
static Prof *PROF;
static double *TOTAL;

/* one entry of the distance matrix, exactly as stored (float32) */
static float dist_pair(size_t a, size_t b) {
    if (a == b) return 0.0f;
    double sh = shared_min(&PROF[a], &PROF[b]);
    double un = (TOTAL[a] + TOTAL[b]) - sh;
    return (float)(un > 0 ? 1 - sh / un : 1.0);
}

/* ---------- UPGMA ---------- */
static float *CD;
static char *ACTIVE;
static double *NND;
static long *NNJ;

static void scan(size_t i) {
    float best = INFINITY;
    long bj = -1;
    const float *row = CD + i * N;
    for (size_t j = i + 1; j < N; j++)
        if (ACTIVE[j] && row[j] < best) { best = row[j]; bj = (long)j; }
    if (bj >= 0) { NND[i] = (double)best; NNJ[i] = bj; }
    else { NND[i] = INFINITY; NNJ[i] = -1; }
}

static int *upgma_order(size_t *olen) {
    size_t n = N;
    int *order = xmalloc((n ? n : 1) * sizeof(int));
    if (n < 2) { for (size_t i = 0; i < n; i++) order[i] = (int)i; *olen = n; return order; }
    int **cl = xmalloc(n * sizeof(int *));
    size_t *clen = xmalloc(n * sizeof(size_t));
    for (size_t i = 0; i < n; i++) { cl[i] = xmalloc(sizeof(int)); cl[i][0] = (int)i; clen[i] = 1; }
    ACTIVE = xmalloc(n);
    memset(ACTIVE, 1, n);
    NND = xmalloc(n * sizeof(double));
    NNJ = xmalloc(n * sizeof(long));
    for (size_t i = 0; i < n; i++) scan(i);

    for (size_t it = 0; it + 1 < n; it++) {
        long ci = -1;
        double bd = INFINITY;
        for (size_t i = 0; i < n; i++)
            if (ACTIVE[i] && NNJ[i] >= 0 && NND[i] < bd) { bd = NND[i]; ci = (long)i; }
        if (ci < 0) break;
        long cj = NNJ[ci];
        int *A = cl[ci], *B = cl[cj];
        size_t sI = clen[ci], sJ = clen[cj];
        float opts[4] = {dist_pair((size_t)A[sI - 1], (size_t)B[0]), dist_pair((size_t)A[sI - 1], (size_t)B[sJ - 1]),
                         dist_pair((size_t)A[0], (size_t)B[0]), dist_pair((size_t)A[0], (size_t)B[sJ - 1])};
        int best = 0;
        for (int t = 1; t < 4; t++) if (opts[t] < opts[best]) best = t;
        if (best >= 2) for (size_t x = 0, y = sI - 1; x < y; x++, y--) { int t = A[x]; A[x] = A[y]; A[y] = t; }
        if (best == 1 || best == 3) for (size_t x = 0, y = sJ - 1; x < y; x++, y--) { int t = B[x]; B[x] = B[y]; B[y] = t; }
        A = realloc(A, (sI + sJ) * sizeof(int));
        if (!A) die("out of memory");
        memcpy(A + sI, B, sJ * sizeof(int));
        free(B);
        cl[ci] = A; clen[ci] = sI + sJ; cl[cj] = NULL; clen[cj] = 0;
        ACTIVE[cj] = 0;

        double dI = (double)sI, dJ = (double)sJ, dS = (double)(sI + sJ);
        for (size_t j = 0; j < n; j++) {
            if (!ACTIVE[j] || (long)j == ci) continue;
            float v = (float)(((double)CD[ci * n + j] * dI + (double)CD[cj * n + j] * dJ) / dS);
            CD[ci * n + j] = v;
            CD[j * n + ci] = v;
        }
        scan((size_t)ci);
        for (size_t i = 0; i < n; i++) {
            if (!ACTIVE[i] || (long)i == ci) continue;
            if (NNJ[i] == ci || NNJ[i] == cj) scan(i);
            else if ((long)i < ci) {
                double v = (double)CD[i * n + ci];
                if (v < NND[i] || (v == NND[i] && ci < NNJ[i])) { NND[i] = v; NNJ[i] = ci; }
            }
        }
    }
    size_t f = 0;
    while (f < n && !ACTIVE[f]) f++;
    if (f == n) die("internal error: no active cluster");
    *olen = clen[f];
    memcpy(order, cl[f], clen[f] * sizeof(int));
    for (size_t i = 0; i < n; i++) free(cl[i]);
    free(cl); free(clen);
    return order;
}

int main(int argc, char **argv) {
    if (argc != 6) die("usage: kmer_order in.fasta out.fasta K CANONICAL(0/1) ORIENT(0/1)");
    int k = atoi(argv[3]);
    int canonical = strcmp(argv[4], "1") == 0, orient = strcmp(argv[5], "1") == 0;
    if (k < 3 || k > 12) die("K must be from 3 to 12");
    init_tables();
    size_t n;
    Rec *recs = read_fasta(argv[1], &n);
    N = n;

    /* k-mer codes, counts, and how many sequences hold each k-mer */
    size_t nk = (size_t)1 << (2 * k);
    uint32_t *seen = xcalloc(nk, sizeof(uint32_t));
    uint32_t **codes = xmalloc((n ? n : 1) * sizeof(uint32_t *));
    size_t *ncodes = xmalloc((n ? n : 1) * sizeof(size_t));
    TOTAL = xmalloc((n ? n : 1) * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        codes[i] = xmalloc((recs[i].len + 1) * sizeof(uint32_t));
        ncodes[i] = kmer_codes(recs[i].seq, recs[i].len, k, canonical, codes[i]);
        TOTAL[i] = (double)ncodes[i];
        qsort(codes[i], ncodes[i], sizeof(uint32_t), cmp_u32);
        for (size_t a = 0; a < ncodes[i]; a++)
            if (a == 0 || codes[i][a] != codes[i][a - 1]) seen[codes[i][a]]++;
    }
    PROF = xmalloc((n ? n : 1) * sizeof(Prof));
    for (size_t i = 0; i < n; i++) {
        Prof *p = &PROF[i];
        p->code = xmalloc((ncodes[i] + 1) * sizeof(uint32_t));
        p->cnt = xmalloc((ncodes[i] + 1) * sizeof(uint32_t));
        p->m = 0;
        for (size_t a = 0; a < ncodes[i];) {
            size_t b = a;
            while (b < ncodes[i] && codes[i][b] == codes[i][a]) b++;
            if (seen[codes[i][a]] >= 2) { p->code[p->m] = codes[i][a]; p->cnt[p->m] = (uint32_t)(b - a); p->m++; }
            a = b;
        }
        free(codes[i]);
    }
    free(codes); free(ncodes); free(seen);

    /* distance matrix */
    CD = xmalloc((n * n ? n * n : 1) * sizeof(float));
    #pragma omp parallel for schedule(dynamic, 16)
    for (long i = 0; i < (long)n; i++) {
        CD[(size_t)i * n + (size_t)i] = 0.0f;
        for (size_t j = (size_t)i + 1; j < n; j++) {
            float v = dist_pair((size_t)i, j);
            CD[(size_t)i * n + j] = v;
            CD[j * n + (size_t)i] = v;
        }
    }

    size_t olen;
    int *snapshot_order = upgma_order(&olen);
    int *order = snapshot_order;
    if (olen != n) die("internal error: order length");

    char *flip = xcalloc(n, 1);
    if (orient) {
        uint32_t *mark = xcalloc(nk, sizeof(uint32_t));
        size_t maxlen = 1;
        for (size_t i = 0; i < n; i++) if (recs[i].len > maxlen) maxlen = recs[i].len;
        uint32_t *cf = xmalloc(maxlen * sizeof(uint32_t)), *cr = xmalloc(maxlen * sizeof(uint32_t)), *cp = xmalloc(maxlen * sizeof(uint32_t));
        char *rc = xmalloc(maxlen + 1);
        size_t np = 0;
        int have_prev = 0;
        for (size_t t = 0; t < n; t++) {
            size_t i = (size_t)order[t];
            size_t L = recs[i].len;
            size_t nf = kmer_codes(recs[i].seq, L, k, 0, cf);
            for (size_t x = 0; x < L; x++) rc[x] = (char)COMP[(unsigned char)recs[i].seq[L - 1 - x]];
            rc[L] = 0;
            size_t nr = kmer_codes(rc, L, k, 0, cr);
            if (have_prev) {
                size_t fwd = 0, rev = 0;
                uint32_t stamp = (uint32_t)t;           /* marks of the predecessor carry stamp t */
                for (size_t x = 0; x < nf; x++) fwd += (mark[cf[x]] == stamp);
                for (size_t x = 0; x < nr; x++) rev += (mark[cr[x]] == stamp);
                flip[i] = rev > fwd;
            }
            uint32_t *chosen = flip[i] ? cr : cf;
            size_t nc = flip[i] ? nr : nf;
            memcpy(cp, chosen, nc * sizeof(uint32_t));
            np = nc;
            for (size_t x = 0; x < np; x++) mark[cp[x]] = (uint32_t)(t + 1);
            have_prev = 1;
        }
        free(mark); free(cf); free(cr); free(cp); free(rc);
    }

    FILE *out = fopen(argv[2], "wb");
    if (!out) die("cannot open output");
    for (size_t t = 0; t < n; t++) {
        size_t i = (size_t)order[t];
        if (flip[i]) {
            fprintf(out, ">_R_%s\n", recs[i].name);
            for (size_t x = 0; x < recs[i].len; x++) fputc((int)COMP[(unsigned char)recs[i].seq[recs[i].len - 1 - x]], out);
            fputc('\n', out);
        } else {
            fprintf(out, ">%s\n%s\n", recs[i].name, recs[i].seq);
        }
    }
    if (fclose(out) != 0) die("write error");
    return 0;
}

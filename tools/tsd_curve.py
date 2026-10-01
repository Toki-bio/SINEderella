import sys, math
f=sys.argv[1]; left=int(sys.argv[2]); right=int(sys.argv[3])
seqs=[]
for l in open(f):
    l=l.rstrip()
    if l.startswith(">"): seqs.append("")
    else: seqs[-1]+=l
rows=seqs[2:] if seqs[1].startswith("-")*0 else seqs[1:]   # row1 consensus; plate may have a second seed row
rows=[r for r in seqs[1:] if sum(c.isupper() for c in r)>0]
isb=lambda c: c.upper() in "ACGT"
def flank(r):
    up=""; 
    for x in range(left-2,-1,-1):
        c=r[x]
        if isb(c): up=c.upper()+up
        if len(up)>=30: break
    dn=""
    for x in range(right,len(r)):
        c=r[x]
        if isb(c): dn+=c.upper()
        if len(dn)>=56: break
    return up,dn
def mm(a,b): return sum(0.5 if "N" in (x,y) else (x!=y) for x,y in zip(a,b))
def best(up,dn,mn,rslk):
    TL=0; TSC=-1e9; ts=0
    for L in range(mn,21):
        if len(up)<L or len(dn)<L: continue
        cl=len(up)-L
        for uo in range(0,min(cl,4)+1):
            us=up[cl-uo:cl-uo+L]
            for ds in range(0,min(len(dn)-L,rslk)+1):
                dv=mm(us,dn[ds:ds+L])/L
                if dv>0.20: continue
                sc=(1-dv)*math.sqrt(L)-uo*0.04-ds*0.13
                if sc>TSC or (sc==TSC and L>TL): TSC=sc;TL=L;ts=ds
    return TL,ts
fl=[flank(r) for r in rows]; n=len(fl)
print("copies",n,"body columns",left,right)
for rslk in (3,25):
    print("3' slack",rslk)
    for mn in (4,6,8,10,12,14):
        hr=sum(1 for u,d in fl if best(u,d,mn,rslk)[0]); hs=sum(1 for i in range(n) if best(fl[i][0],fl[(i+1)%n][1],mn,rslk)[0])
        print("  min %2d: real %5.1f%%  shuffled %5.1f%%  excess %5.1f"%(mn,100*hr/n,100*hs/n,100*(hr-hs)/n))
# where the real TSD sits at min 8
import collections
c=collections.Counter(best(u,d,8,25)[1] for u,d in fl if best(u,d,8,25)[0]); print("3' start offset after the body (min 8):",sorted(c.items())[:12])

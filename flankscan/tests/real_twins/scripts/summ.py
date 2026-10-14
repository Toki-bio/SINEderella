import sys, statistics as st
def load(p):
    d={}
    for i,l in enumerate(open(p)):
        if i==0: continue
        f=l.rstrip("\n").split("\t"); k=(f[0],f[1]); fl,b=float(f[3]),float(f[6])
        if k not in d or fl>d[k][0]: d[k]=(fl,b,f[5])
    return d
rnd=load(sys.argv[1]); rb=sorted(v[1] for v in rnd.values()); p95=rb[int(0.95*len(rb))]
print(f"background: random pairs median {st.median(rb):.1f}, p95 {p95}")
sets=[(n,load(p)) for n,p in zip(sys.argv[2::2], sys.argv[3::2])]
seen=set()
for n,d in sets:
    if n.startswith("only:"): d={k:v for k,v in d.items() if k not in seen}; 
    seen|=set(d)
    b=[v[1] for v in d.values()]; a=sum(1 for x in b if x>p95); arr=sum(1 for v in d.values() if v[2]=="array")
    print(f"{n:26s} pairs {len(d):5d} (array {arr:4d})  body median {st.median(b):5.1f}  body > p95: {a:5d} ({100*a/max(1,len(d)):3.0f}%)")

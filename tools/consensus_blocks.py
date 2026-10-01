# consensus_blocks.py BANK.fa OUT.json : all-against-all ungapped similarity blocks (>=20 bp, >=78 %, both strands) between consensuses; used for the rsi block matrix (rsi page).
import numpy as np, re, json, sys
def rd(f):
    n=[];s=[]
    for l in open(f):
        l=l.strip()
        if not l: continue
        if l[0]=='>': n.append(l[1:]);s.append('')
        else: s[-1]+=l.upper()
    return n,s
comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}
rc=lambda x:''.join(comp[c] for c in reversed(x))
def mask(s):
    a=list(s)
    for m in re.finditer(r'(A{7,}|T{7,}|C{7,}|G{7,}|(?:AT){5,}|(?:TA){5,}|(?:CA){5,}|(?:TG){5,}|(?:GA){5,}|(?:TC){5,})',s):
        for i in range(m.start(),m.end()): a[i]='N'
    return ''.join(a)
def segs(v,T):
    # all maximal-scoring segments with score>=T, by repeated max-subarray
    out=[]; stack=[(0,len(v))]
    while stack:
        lo,hi=stack.pop()
        if hi-lo<T: continue
        best=0;bs=be=0;cur=0;cs=lo
        for i in range(lo,hi):
            if cur<=0: cur=0; cs=i
            cur+=v[i]
            if cur>best: best=cur;bs=cs;be=i+1
        if best>=T:
            out.append((bs,be,best)); stack.append((lo,bs)); stack.append((be,hi))
    return out
def blocks(a,b,strand,T=14,minlen=20,minid=0.78):
    A=np.frombuffer(a.encode(),dtype=np.uint8); B=np.frombuffer(b.encode(),dtype=np.uint8)
    res=[]
    for d in range(-(len(A)-1),len(B)):
        i0=max(0,-d); j0=i0+d; L=min(len(A)-i0,len(B)-j0)
        if L<minlen: continue
        x=A[i0:i0+L]; y=B[j0:j0+L]
        v=np.where((x==y)&(x!=78),1,-2)
        if (v>0).sum()<T: continue
        for s,e,sc in segs(v.tolist(),T):
            seg=v[s:e]; m=int((seg>0).sum()); ln=e-s
            if ln>=minlen and m/ln>=minid: res.append((i0+s,i0+e,j0+s,j0+e,m,ln,strand))
    return res
def pick(cands):
    cands.sort(key=lambda c:-c[4]); acc=[]
    def ov(a0,a1,b0,b1): return max(0,min(a1,b1)-max(a0,b0))
    for c in cands:
        dup=False
        for k in acc:
            if k[6]!=c[6]: continue
            if ov(c[0],c[1],k[0],k[1])>0.5*(c[1]-c[0]) and ov(c[2],c[3],k[2],k[3])>0.5*(c[3]-c[2]): dup=True;break
        if not dup: acc.append(c)
    acc.sort(key=lambda c:c[0]); return acc
def merge(acc):
    acc=sorted(acc,key=lambda c:(c[6],c[0])); out=[]
    for c in acc:
        if out and out[-1][6]==c[6] and 0<=c[0]-out[-1][1]<=8 and ((c[6]=='+' and 0<=c[2]-out[-1][3]<=8) or (c[6]=='-' and 0<=out[-1][2]-c[3]<=8)):
            k=out[-1]; out[-1]=(k[0],c[1],min(k[2],c[2]),max(k[3],c[3]),k[4]+c[4],(c[1]-k[0]),k[6])
        else: out.append(c)
    return out
def compute(path):
    """all-against-all blocks of the consensuses in a FASTA: {'names', 'lens', 'blocks': [i, j, strand, a0, a1, b0, b1, identity]}"""
    n,s=rd(path); ms=[mask(x) for x in s]
    data={'names':n,'lens':[len(x) for x in s],'blocks':[]}
    for i in range(len(n)):
        for j in range(i+1,len(n)):
            c=blocks(ms[i],ms[j],'+')
            r=rc(ms[j]); cr=blocks(ms[i],r,'-')
            # convert minus coords to B forward coordinates
            L=len(r); cr=[(a0,a1,L-b1,L-b0,m,ln,'-') for a0,a1,b0,b1,m,ln,_ in cr]
            bl=merge(pick(c+cr))
            for a0,a1,b0,b1,m,ln,st in bl:
                data['blocks'].append([i,j,st,a0+1,a1,b0+1,b1,round(100*m/ln,1)])
    return data
if __name__=='__main__':
    data=compute(sys.argv[1])
    json.dump(data,open(sys.argv[2],'w'))
    print(len(data['blocks']),'blocks')

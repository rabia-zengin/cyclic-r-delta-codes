"""Compare every optimal parameter set of global type found by the search
with the known constructions listed in Table II of the paper.

Usage:  python3 python/compare_literature.py
Output: results/literature_comparison.csv

For each parameter set (q, N, K, d, r_H, delta_H) the script checks the
parameter conditions of
  cyclic:     Sun-Zhu-Wang 2019 (Thms 1-5), Tian-Zhu 2022 (Thms 4-7),
              Qian-Zhang 2020, Chen-Kai 2022, Fang-Fu 2020, Luo-Xing-Yuan 2019
  non-cyclic: Zhang-Liu 2020 / Kong-Wang-Ge 2021 (Cor. III.2),
              Wang-Yang 2023 (Thm 1 with nested MDS codes of length <= q+1),
              Cai-Miao-Schwartz-Tang 2020 (Cor. 4, Construction A),
              Xing-Yuan 2022 / Kong-Wang-Ge 2021 (necessary set conditions),
              Galindo-Hernando-Martin-Cruz 2024 (Cor. 4.5, 4.9, Thm 5.11),
              Luo-Ezerman-Ling 2023 (Thm 4, Cor. 1, 2)
and records the first construction whose conditions hold.

The remaining constructions listed in the paper are not tested one by one,
because their basic conditions never hold for the parameter sets considered;
this is checked at the end of the script:
  locality r = 2 (Fang et al. 2023, Qiu et al. 2022): no set has r_H = 2;
  repeated-root codes (Zhao et al.): gcd(N, q) = 1 for every set;
  lengths dividing q - 1 or q + 1 (Tamo et al., Chen et al., Qiu et al.):
      excluded by the selection N does not divide q +- 1;
  curve-based codes (Gao-Shi-Fang 2025, Huang-Zhao 2025): need q to be a
      square, or give lengths shorter than N for prime q.
"""
import os
import csv,glob
from math import gcd
sets={}
for p in glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data','search_q*.csv')):
    for r in csv.DictReader(open(p)):
        if not r['d']: continue
        q,N,K,d,rH,dH=int(r['q']),int(r['N']),int(r['K']),int(r['d']),int(r['r_H']),int(r['delta_H'])
        if d==int(r['B_H']) and d>=3 and rH>=2 and K>rH and (q-1)%N and (q+1)%N:
            sets.setdefault((q,N,K,d,rH,dH),[]).append(r)
def qzck(q,N,K,d,rH,dH):
    if q%2==0 or dH%2 or d!=dH+2: return None
    if N==2*(q+1) and rH+dH-1==q+1 and 2<=dH<=q-1 and K==2*(q+1)-2*dH: return 'QZ'
    if N==2*(q+1) and q>=7 and q%4==3 and rH+dH-1==(q+1)//2 and 2<=dH<=(q-3)//2 and K==2*(q+1)-4*dH+2: return 'QZ'
    if N==3*(q+1) and q%6==1 and rH+dH-1==q+1 and K==3*(q+2-dH)-2: return 'CK'
    return None
def szw(q,N,K,d,rH,dH):
    r,dl=rH,dH; L=r+dl-1
    if (q+1)%L or N%L or r<2 or dl<2: return None
    u=N//L; g=gcd(u,L)
    if dl%2==0:
        e=(d-dl)//2
        if (d-dl)%2==0 and 1<=e<=dl//2 and r>=2*e+1 and g==1 and K==r*u-2*e: return 'SZW1'
        if r>=3 and (dl//2)%g==0 and d==dl+2 and K==r*u-2: return 'SZW2'
    else:
        if dl>=3:
            e=(d-dl)//2
            if (d-dl)%2==0 and 1<=e<=(dl-1)//2 and r>=2*e+1 and g==1 and K==r*u-2*e: return 'SZW3c'
            if r>=3 and dl%g==0 and d==dl+2 and K==r*u-2: return 'SZW4c'
            for e in (0,1):
                if r>=dl+2 and r%2==1 and g==1 and d==2*dl+2*e-2 and K==r*u-(dl-2+2*e): return 'SZW5'
    return None
def ffq(q,N,K,d,rH,dH):
    """Fang-Fu (q-1 families, and the q+1, d=6 family), Luo-Xing-Yuan; Tian-Zhu treated broadly"""
    r,dl=rH,dH; L=r+dl-1
    if N%L: return None
    u=N//L; g=gcd(u,L)
    if (q-1)%L==0:
        if d==dl+1 and K==r*u-1: return 'FF'
        if r>=3 and dl%g==0 and d==dl+2 and K==r*u-2: return 'FF'
        if r>=dl+1 and g==1 and d==2*dl and K==r*u-dl: return 'FF'
        if dl<d<=2*dl and K==r*u-(d-dl): return 'TZ?'   # Tian-Zhu: conditions not checked
    if (q+1)%L==0 and N%2==1 and dl==3 and r>=5 and g==1 and d==6 and K==r*u-3: return 'FF'
    return None
def divisors(x): return [d for d in range(1,x+1) if x%d==0]
def tz(q,N,K,d,rH,dH):
    """Tian-Zhu 2022, Theorems 4-7 ((r+delta-1) | gcd(N,q+1), rho = rho1*rho2, rho2*(r+delta-1) | q+1)"""
    r,dl=rH,dH; L=r+dl-1
    if (q+1)%L or N%L: return None
    rho=N//L
    for r2 in divisors(rho):
        if (q+1)%(r2*L): continue
        g=gcd(rho//r2,L)
        if dl%2==0:
            if r>=3 and (dl//2)%g==0 and d==dl+2 and K==r*rho-2: return 'TZ4'
            e=(d-dl)//2
            if (d-dl)%2==0 and 1<=e<=dl//2 and r>=2*e+1 and g==1 and K==r*rho-2*e: return 'TZ5'
        else:
            if r>=3 and r%2==1 and ((r-1)//2)%g==0 and d==dl+2 and K==r*rho-2: return 'TZ6'
            e=(d-dl)//2
            if (d-dl)%2==0 and 1<=e<=(dl-1)//2 and r>=2*e+1 and r%2==1 and g==1 and K==r*rho-2*e: return 'TZ7'
    return None
def cai4(q,N,K,d,rH,dH):
    """Cai-Miao-Schwartz-Tang 2020, Corollary 4 (non-cyclic): X=r-v<=delta, q>=L+X"""
    L=rH+dH-1
    if N%L: return None
    w=N//L; X=w*rH-K
    if 0<=X<=dH and X<=rH-1 and d==dH+X and q>=L+X: return 'Cai4'
    return None
def xy_possible(q,N,K,d,rH,dH):
    """Xing-Yuan 2022: necessary conditions for their Vandermonde construction (d>=2delta+1)"""
    L=rH+dH-1; m=N//L
    if d<2*dH+1: return False
    if L>q or d>rH+dH: return False
    from math import comb
    return m<=comb(q,dH)//comb(L,dH)
def chen22_ok(q,N,K,d,rH,dH):
    """upper bound of Chen et al. ICCC 2022 (Lemma 3), valid for d>=2delta+1 with disjoint groups"""
    from math import comb
    L=rH+dH-1
    if d<2*dH+1: return True
    u=d-1-((d-dH)//L+1)*(dH-1)
    return N//L <= (q**u-1)//((q-1)*comb(L,dH))
def caiA_possible(q,N,K,d,rH,dH):
    """necessary conditions for Construction A of Cai et al. (any corollary of it)"""
    L=rH+dH-1; w=N//L; X=w*rH-K
    if not (0<=X<=rH-1) or q<L+X: return False
    if X<=dH: return True
    return (w*(w-1)//2)*(dH-1) >= w*L-(q-X)
def xy_possible2(q,N,K,d,rH,dH):
    L=rH+dH-1; m=N//L
    if d<2*dH+1 or L>q or d>rH+dH: return False
    return (m*(m-1)//2)*(dH-1) >= m*L-q
def kong(q,N,K,d,rH,dH):
    """Kong-Wang-Ge Cor. III.2 / Zhang-Liu Thm: d in [delta+1,2delta], r>=d-delta, L<=q, any length"""
    L=rH+dH-1; m=N//L
    return N%L==0 and dH+1<=d<=2*dH and rH>=d-dH and L<=q and K==rH*m-(d-dH)
def wy(q,N,K,d,rH,dH):
    """Wang-Yang Thm 1 with nested MDS codes of length L<=q+1 (cyclic/constacyclic BCH-type MDS codes)"""
    L=rH+dH-1; m=N//L
    return N%L==0 and dH<=d<=2*dH and rH>d-dH and L<=q+1 and K==rH*m-(d-dH)

def _factorable(P,q,parts):
    # can P be written as a product of `parts` integers in [2,q]?
    if parts==1: return 2<=P<=q
    return any(P%f==0 and _factorable(P//f,q,parts-1) for f in range(2,min(P,q)+1))
def galindo(q,N,K,d,rH,dH):
    """Galindo-Hernando-Martin-Cruz 2024: monomial-Cartesian codes over F_q (Cor. 4.5, 4.9)
       and subfield-subcodes with block p^h+1 (Thm 5.11, odd characteristic, h=1)"""
    L=rH+dH-1
    if N%L: return None
    n1=L; P=N//L; i=rH-1
    if 2<=n1<=q:
        # bivariate, n2=P
        n2=P
        if 2<=n2<=q:
            if i==0 and dH==n1:
                for j in range(n2):
                    if K==j+1 and d==n1*(n2-j): return 'G4.5(1)'
            if 1<=i<=n1-2 and dH==n1-i:
                if K==(i+1)*n2 and d==n1-i: return 'G4.5(1)'
                for s in range(max(0,2*i-n1),i):
                    if K==(i+1)*(n2-1)+s+1 and d==n1-s: return 'G4.5(2)'
                for j in range(1,n2-1):
                    if j>=max(1,(i*(n2+1)-n1)/i) and K==(i+1)*j+1 and d==n1*(n2-j): return 'G4.5(3)'
        # multivariate: P = product of m-1>=2 factors in [2,q]
        if any(_factorable(P,q,t) for t in range(2,8)):
            if 0<=i<=n1-2 and dH==n1-i and K==P*(i+1) and d==n1-i: return 'G4.9(1)'
            for s in range(max(1,2*i-n1+1),i+1):
                if 0<=i<=n1-2 and dH==n1-i and K==P*(i+1)-(i-s+1) and d==n1-s+1: return 'G4.9(2)'
            if i==0 and dH==n1 and K==P-1 and d==2*n1: return 'G4.9(3)'
    # subfield-subcodes (Thm 5.11, p odd, h=1): L=q+1, r=2z+1, delta=q-2z+1
    if L==q+1 and q%2==1 and rH%2==1:
        z=(rH-1)//2; n2=P
        if 0<=z<=q//2-1:
            for a in [z]+[t for t in range(0,z) if 2*t>=max(0,4*z-q-1)]:
                if K==(n2-1)*(2*z+1)+2*a+1 and d==q+1-2*a: return 'G5.11'
    return None

def luo_mp(q,N,K,d,rH,dH):
    """Luo-Ezerman-Ling 2023, Thm 4 / Cor. 1, 2 (matrix-product codes):
       [Nb*L, (M-1)r+1, (Nb-M+1)L], q >= Nb > M > 1; L<=q (Cor 1) or L=q+1, r odd, (r-1)(Nb-M+1)+r<=q+2 (Cor 2)"""
    L=rH+dH-1
    if N%L: return None
    Nb=N//L
    if not (Nb<=q): return None
    for M in range(2,Nb):
        if K==(M-1)*rH+1 and d==(Nb-M+1)*L:
            if L<=q: return 'LEL-Cor1'
            if L==q+1 and rH%2==1 and (rH-1)*(Nb-M+1)+rH<=q+2: return 'LEL-Cor2'
    return None

def classify(k):
    for name,f in (('Sun-Zhu-Wang',szw),('Tian-Zhu',tz),('Qian-Zhang/Chen-Kai',qzck),('Fang-Fu/Luo-Xing-Yuan',ffq)):
        if f(*k): return 'cyclic',name
    for name,f in (('Zhang-Liu/Kong-Wang-Ge',kong),('Wang-Yang',wy),('Cai et al.',cai4),('Galindo et al.',galindo),('Luo-Ezerman-Ling',luo_mp)):
        if f(*k): return 'non-cyclic',name
    if caiA_possible(*k) or xy_possible2(*k): return 'possibly non-cyclic','Cai et al./Xing-Yuan'
    return 'new','-'

if __name__=='__main__':
    out=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','results','literature_comparison.csv')
    G=sorted(k for k,v in sets.items() if v[0]['type']=='G')
    from collections import Counter
    cnt=Counter()
    with open(out,'w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['q','N','K','d','r_H','delta_H','s','n','delta','a','b','status','construction'])
        for k in G:
            r=min(sets[k],key=lambda x:(int(x['s']),int(x['delta']),int(x['a']),int(x['b'])))
            st,name=classify(k); cnt[st]+=1
            w.writerow(list(k)+[r['s'],r['n'],r['delta'],r['a'],r['b'],st,name])
    print(len(G),'optimal parameter sets of type G:',dict(cnt))
    from math import gcd as _gcd
    assert all(k[4] != 2 for k in G), "a set with r_H = 2"
    assert all(_gcd(k[1], k[0]) == 1 for k in G), "a set with gcd(N, q) > 1"
    assert all((k[0] - 1) % k[1] and (k[0] + 1) % k[1] for k in G), "a length dividing q +- 1"
    print("checked: no set has r_H = 2, gcd(N, q) > 1 or N | (q +- 1)")

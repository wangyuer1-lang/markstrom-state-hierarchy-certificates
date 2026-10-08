"""Finite checks for Strict State Hierarchies for Graph Homomorphism Polynomials.
Python 3.9+; standard library only.
"""
from itertools import permutations, product, combinations
from math import factorial, prod
from collections import deque
from pathlib import Path
import hashlib, json

P=1000003
HERE=Path(__file__).resolve().parent

def parity(p):
    return sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2

def glued(r, plus, minus):
    adj=[set() for _ in range(r)]
    def edge(u,v):
        assert u!=v and v not in adj[u]
        adj[u].add(v);adj[v].add(u)
    for kind,count in ((0,plus),(1,minus)):
        for _ in range(count):
            for p in permutations(range(r)):
                if parity(p)!=kind: continue
                c=len(adj);adj.append(set())
                for i,t in enumerate(p,1):
                    for _ in range(i):
                        w=len(adj);adj.append(set());edge(c,w);edge(w,t)
    return adj

def parts(adj):
    side={0:0};queue=[0]
    for u in queue:
        for v in adj[u]:
            if v not in side: side[v]=1-side[u];queue.append(v)
            assert side[v]!=side[u]
    assert len(side)==len(adj)
    return [{v for v,c in side.items() if c==j} for j in (0,1)]

def connected(adj,removed=()):
    removed=set(removed);start=next(v for v in range(len(adj)) if v not in removed)
    seen={start};todo=[start]
    for u in todo:
        for v in adj[u]:
            if v not in removed and v not in seen:seen.add(v);todo.append(v)
    return len(seen)==len(adj)-len(removed)

def augment(adj,t):
    left,right=parts(adj);n=len(adj)
    U=range(n,n+t);V=range(n+t,n+2*t)
    result=[set(x) for x in adj]+[set() for _ in range(2*t)]
    for u in U:
        for v in right|set(V):result[u].add(v);result[v].add(u)
    for v in V:
        for u in left:result[u].add(v);result[v].add(u)
    return result

def matrix(q):
    return [[(min(i,j)+2)*(max(i,j)+3)+3*(i==j)+1 for j in range(q)] for i in range(q)]

def uv_function(r,p,b,A,mod=P):
    q=len(A)
    reduce=(lambda x:x%mod) if mod else (lambda x:x)
    M=[[reduce(sum(A[i][k]*b[k]*A[k][j] for k in range(q))) for j in range(q)] for i in range(q)]
    powers=[[[reduce(M[c][a]**i) for i in range(r+1)] for a in range(q)] for c in range(q)]
    perms=[(perm,parity(perm)) for perm in permutations(range(r))]
    def uv(a):
        result=[1,1]
        for perm,kind in perms:
            L=reduce(sum(p[c]*prod(powers[c][a[perm[i]]][i+1] for i in range(r)) for c in range(q)))
            result[kind]=reduce(result[kind]*L)
        return result
    return uv

def conditional_glued(r,plus,minus,p,b,A,mod=P):
    reduce=(lambda x:x%mod) if mod else (lambda x:x)
    uv=uv_function(r,p,b,A,mod)
    total=0
    for a in product(range(len(A)),repeat=r):
        U,V=uv(a)
        total=reduce(total+prod(p[c] for c in a)*U**plus*V**minus)
    return total

def eliminate(adj,h,A):
    q=len(h)
    factors=[((v,),{(s,):h[s] for s in range(q)}) for v in range(len(adj))]
    for v,ns in enumerate(adj):
        for w in ns:
            if v<w:factors.append(((v,w),{(i,j):A[i][j] for i in range(q) for j in range(q)}))
    remaining=set(range(len(adj)))
    while remaining:
        neighbors={v:set() for v in remaining}
        for scope,_ in factors:
            for v in scope:neighbors[v].update(x for x in scope if x!=v)
        v=min(remaining,key=lambda x:(len(neighbors[x]),x))
        selected=[f for f in factors if v in f[0]]
        factors=[f for f in factors if v not in f[0]]
        scope=tuple(sorted(set().union(*(set(s) for s,_ in selected))-{v}))
        positions=[tuple(-1 if x==v else scope.index(x) for x in s) for s,_ in selected]
        table={}
        for assignment in product(range(q),repeat=len(scope)):
            table[assignment]=sum(prod(f[tuple(c if i<0 else assignment[i] for i in ix)]
                                       for (_,f),ix in zip(selected,positions))
                                  for c in range(q))%P
        factors.append((scope,table));remaining.remove(v)
    return prod(f[()] for _,f in factors)%P

def diagonal_gap(r,eta,d):
    vals=[1,1]
    for perm in permutations(range(r)):
        e=[perm.index(c)+1 for c in range(r)]
        vals[parity(perm)]*=sum(eta[c]*d[c]**e[c] for c in range(r))
    return vals[0]-vals[1]

def integer_digest(n):
    s=str(n)
    return {'digits':len(s.lstrip('-')),'sign':(n>0)-(n<0),'mod_P':n%P,
            'sha256_decimal':hashlib.sha256(s.encode()).hexdigest()}


def verify_archived_inputs():
    graphs = json.loads((HERE/'graphs.json').read_text(encoding='utf-8'))
    for name, data in graphs.items():
        r, t, d = data['r'], data['t'], data['d']
        adj = glued(r, t+d, t-d)
        if data['augmentation']:
            adj = augment(adj, data['augmentation'])
        edges = [[u,v] for u in range(len(adj)) for v in sorted(adj[u]) if u<v]
        assert data['vertices'] == len(adj), name
        assert data['edges'] == edges, name
        assert data['bipartition'] == [sorted(s) for s in parts(adj)], name
    examples = json.loads((HERE/'exact_values.json').read_text(encoding='utf-8'))
    for example in examples:
        r, t, A, h = example['r'], example['t'], example['A'], example['h']
        values = [conditional_glued(r, t+d, t-d, h, h, A, mod=None) for d in range(t+1)]
        assert [str(v) for v in values] == example['values_decimal']
        assert [str(b-a) for a,b in zip(values,values[1:])] == example['adjacent_differences_decimal']

def run():
    report={'prime':P,'q1_endpoint':[],'many_way_split':[],'potts_augmentation':[],
            'diagonal_certificate':[]}
    r=2
    pair=[glued(r,2,0),glued(r,1,1)]
    for adj in pair:
        assert len(adj)==10 and sum(map(len,adj))//2==12
        parts(adj)
        assert all(connected(adj,(v,)) for v in range(len(adj)))
    for q in (1,2):
        A=matrix(q);h=[1]*q
        direct=[eliminate(g,h,A) for g in pair]
        expected=[conditional_glued(r,a,b,h,h,A) for a,b in ((2,0),(1,1))]
        assert direct==expected and (direct[0]==direct[1])==(q==1)
        report['q1_endpoint'].append({'q':q,'A':A,'direct':direct,'vertices':10,'edges':12})
    for q in (2,3,4):
        A=[[5 if i==j else 1 for j in range(q)] for i in range(q)];h=[1]*q
        assert eliminate(pair[0],h,A)==eliminate(pair[1],h,A)

    r=3;t=2
    family=[glued(r,t+d,t-d) for d in range(t+1)]
    for q in (2,3):
        A=matrix(q);h=[1]*q
        direct=[eliminate(g,h,A) for g in family]
        exact=[conditional_glued(r,t+d,t-d,h,h,A,mod=None) for d in range(t+1)]
        assert direct==[v%P for v in exact]
        if q<r:assert len(set(exact))==1
        else:
            assert all(a<b for a,b in zip(exact,exact[1:]))
            U,V=uv_function(r,h,h,A,mod=None)(tuple(range(r)))
            for d in range(t):
                gap=factorial(r)//2*prod(h)*(U*V)**(t-d-1)*(U-V)**2*sum(U**(2*d-i)*V**i for i in range(2*d+1))
                assert exact[d+1]-exact[d]==gap
        report['many_way_split'].append({'r':r,'t':t,'q':q,'vertices':len(family[0]),
            'direct':direct,'exact_integer_values':[integer_digest(v) for v in exact]})
    A=[[5 if i==j else 1 for j in range(4)] for i in range(4)];h=[1]*4
    vals=[eliminate(g,h,A) for g in family]
    assert len(set(vals))==1
    report['many_way_split'].append({'r':r,'t':t,'q':4,'Potts':True,'values':vals})

    # Tutte-equivalence does not survive all bipartite augmentations.
    # r=3,t=1 is outside the sufficient bound r>2t+1.
    for r,q in ((3,3),(4,4)):
        pair=[glued(r,2,0),glued(r,1,1)];t=1
        A=[[5 if i==j else 1 for j in range(q)] for i in range(q)];h=[1]*q
        augmented=[augment(g,t) for g in pair]
        direct=[eliminate(g,h,A) for g in augmented]
        predicted=[0,0]
        for u in range(q):
            for v in range(q):
                p=[A[i][v] for i in range(q)];b=[A[i][u] for i in range(q)]
                w=A[u][v]
                for j,(plus,minus) in enumerate(((2,0),(1,1))):
                    predicted[j]=(predicted[j]+w*conditional_glued(r,plus,minus,p,b,A))%P
        assert direct==predicted
        assert (direct[0]==direct[1]) == (r>2*t+1)
        report['potts_augmentation'].append({'r':r,'t':t,'q':q,'v':4,
            'direct':direct,'vertices':len(augmented[0]),'edges':sum(map(len,augmented[0]))//2})

    # Explicit Kronecker specialization. Not necessary for the general proof.
    for r in (2,3):
        k=factorial(r)//2;B=r*k+1;T=r**k+2
        d=[T**(B**i) for i in range(r)]
        D=diagonal_gap(r,[1]*r,d)
        assert D!=0 and all(x>r-1 for x in d)
        report['diagonal_certificate'].append({'r':r,'B':B,'T':T,'d':d,'D':integer_digest(D)})

    verify_archived_inputs()
    report['status']='passed'
    (HERE/'results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()

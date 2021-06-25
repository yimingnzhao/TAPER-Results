import sys
import operator

def read(f):
    a=dict()
    b=dict()
    for line in open(f,'r'):
            if line.startswith(">"):
                    name = line.strip()
            else:
                s = line.strip().replace('X','-')
                a[name]=list(s)
                b[name]=s.replace('-','')
    return (a,b)

r1,r2 = read(sys.argv[1])

d1,d2 = read(sys.argv[2])

cm = [dict() for i in range(0,len(r1[list(r1.keys())[0]]))]
#print (len(cm))
for k,v in r2.items():
    assert v == d2[k], "%s sequences do not match \n%s\n%s" %(k, v, d2[k])
    rs = r1[k]
    ds = d1[k]
    j = 0
    for i,c in enumerate(rs):
        if c == "-":
            continue
        while ds[j] == "-":
            j = j + 1
        assert c == ds[j], "%s %d %s %d %s" %(k,i,c,j,ds[j])
        cm[i][j] = cm[i].get(j,0) + 1
        j = j + 1

#print (cm[18]) 
kc = set()
for i,m in enumerate(cm):
    #print(max(m.values(), default=0))
    kc.add(max(m.items(), key=operator.itemgetter(1),default=(-1,-1))[0])

for k,v in d1.items():
    #print("%s\n%s" %(k,''.join(['X' if i not in kc and c != "-" else c for (i,c) in enumerate(v)])))
    rs = r1[k]
    j = 0
    s=[]
    for (i,c) in enumerate(v):
        if i not in kc:
            continue
        if rs[j] == c:
            s.append(c)
        elif c == '-':
            s.append('X')
        else:
            raise ("Error %i %j %s %s" %(i,j,c,rs[j]))
        j = j + 1
    print("%s\n%s" %(k,''.join(s)))

    

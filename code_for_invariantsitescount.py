with open('fullAlignment.fasta') as f:
    g=f.read()


K=g.split('>')[1:]
LL=[]
for j in K:
    LL.append(j.split('\n')[0])
LL[4]
del(K[4])

O=['Barcelona_JK3030', 'Barcelona', '6330', 'Justinian_A120','0.PE3_Angola', 'outgroup_Y.pseudo', 'Ellwangen']
for k in O:
    LL=[]
    for j in K:
        LL.append(j.split('\n')[0])
    del(K[LL.index(k)])

print('With Reference')
H=[]
S=[]
for i in K:
    H.append(i.split('\n')[0])
    S.append(''.join(i.split('\n')[1:]))

print('No of Sequences:',len(S))
hah={'A':0,'T':0,'G':0,'C':0,'N':0}
l=[]
for m in S:
    l.append(len(m))

for n in range(0,list(set(l))[0],1):
    val=[]
    for i in S:
        val.append(i[n])
    if len(list(set(val)))==1:
        hah[list(set(val))[0]]=hah[list(set(val))[0]]+1

print(hah)

print('Without Reference')
del(K[0])

H=[]
S=[]
for i in K:
    H.append(i.split('\n')[0])
    S.append(''.join(i.split('\n')[1:]))

print('No of Sequences:',len(S))
hah={'A':0,'T':0,'G':0,'C':0,'N':0}
l=[]
for m in S:
    l.append(len(m))

for n in range(0,list(set(l))[0],1):
    val=[]
    for i in S:
        val.append(i[n])
    if len(list(set(val)))==1:
        hah[list(set(val))[0]]=hah[list(set(val))[0]]+1

print(hah)





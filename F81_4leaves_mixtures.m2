restart

K=frac(QQ[p_1,p_2,p_3,p_4]);
RG=K[l_(1,1)..l_(5,2)]
gens RG
-- 12|34 with identity at the leaves
qbar=value get "F81_4leaves_tensor_id.txt";
-- 12|34 no evolution point (star tree with identity at the leaves)
qbarId=sub(qbar,{l_(5,1)=>1,l_(5,2)=>1});

-- Basis of the model
D={0,sub(diagonalMatrix{1,0,0,0},RG),sub(diagonalMatrix{0,1,1,1},RG)}
S=sort elements (set toList(1..4))^**4/splice/splice;

--Replace all 4 and 3 by 2
Ss=apply(S,i->new MutableList from i);
for j from 0 to 255 do for i from 0 to 3 do if (Ss_j)#i==4 then (Ss_j)#i=2
for j from 0 to 255 do for i from 0 to 3 do if (Ss_j)#i==3 then (Ss_j)#i=2
Ss=unique apply(Ss,i->toSequence i)
length Ss --16
netList Ss

--Points obtained by using all combinations of D1,D2 at the leaves
-- with a generic matrix at the interior edge
P=time toList apply(Ss,i->(D_(i_0)**D_(i_1)**D_(i_2)**D_(i_3))*qbar);
-- with the identity matrix at the interior edge (no evolution point)
PId=time toList apply(Ss,i->(D_(i_0)**D_(i_1)**D_(i_2)**D_(i_3))*qbarId);

--12 points (non-zero coordinates) obtained by using all combinations of D1,D2 at the leaves
L=toList apply(0..15,i->S_(positions(flatten entries P_i,j->j!=0)))
netList select(L,i->i!={})
length select(L,i->i!={}) --12
netList Ss_(positions(L,i->i!={}))

LId=toList apply(0..15,i->S_(positions(flatten entries PId_i,j->j!=0)));
length select(LId,i->i!={}) --12
netList select(LId,i->i!={})
netList Ss_(positions(LId,i->i!={}))

--tests
Ss_(positions(LId,i->i!={}))==Ss_(positions(L,i->i!={})) --true
select(LId,i->i!={})==select(L,i->i!={}) --false

--Linearly independent points obtained with the star tree
points=PId_(positions(LId,i->i!={}));
length points

M=points_0;
for i from 1 to 11  do {M=M|points_i};
rank M --12

-- All 12 points are linearly independent (regardless of whether we use Id or a generic matrix at the edge)
-- because they all have different non-zero coordinates

----------------------------------------------------------------------------
-- Additional linearly independent points obtained with tree topology 12|34
----------------------------------------------------------------------------

D1Id=transpose matrix {apply(flatten entries qbar,i->sub(i,{l_(5,1)=>1,l_(5,2)=>0}))};
D2Id=transpose matrix {apply(flatten entries qbar,i->sub(i,{l_(5,1)=>0,l_(5,2)=>1}))};

--Can we obtain new points with D1Id?
PId11=time toList apply(Ss_(positions(LId,i->i!={})),i->(D_(i_0)**D_(i_1)**D_(i_2)**D_(i_3))*D1Id);
apply(0..11,i->rank(points_i|PId11_i))

--Can we obtain new points with D2Id?
PId12=time toList apply(Ss_(positions(LId,i->i!={})),i->(D_(i_0)**D_(i_1)**D_(i_2)**D_(i_3))*D2Id);
apply(0..11,i->rank(points_i|PId12_i))

netList S_(positions(flatten entries PId11_11,i->i!=0))
netList S_(positions(flatten entries PId12_11,i->i!=0))

M=M|PId11_11;
rank(M)--13

----------------------------------------------------------------------------
-- Additional linearly independent points obtained with tree topology 13|24
----------------------------------------------------------------------------
-- 13|24 with identity at the leaves
qbar2=value get "F81_4leaves_tensor_id_1324.txt";
P2=(D_2**D_2**D_2**D_2)*qbar2;
netList S_(positions(flatten entries P2,j->j!=0))
length S_(positions(flatten entries P2,j->j!=0))--25
N=M|P2;
rank(N)--14

--This is for generic qbar2 (i.e. generic Lambda_5)
-- but we can do it for l_(5,1)=1:
D1Id2=transpose matrix {apply(flatten entries qbar2,i->sub(i,{l_(5,1)=>1,l_(5,2)=>0}))};
P2=(D_2**D_2**D_2**D_2)*D1Id2;
netList S_(positions(flatten entries P2,j->j!=0))
length S_(positions(flatten entries P2,j->j!=0))--9
M=M|P2;
rank M --14


----------------------------------------------------------------------------
-- Additional linearly independent points obtained with tree topology 14|23
----------------------------------------------------------------------------
-- 14|23 with identity at the leaves
qbar3=value get "F81_4leaves_tensor_id_1423.txt";
P3=(D_2**D_2**D_2**D_2)*qbar3;
netList S_(positions(flatten entries P3,j->j!=0))
length S_(positions(flatten entries P3,j->j!=0))--25
N=M|P3;
rank(N)--15

--This is for generic qbar2 (i.e. generic Lambda_5)
-- but we can do it for l_(5,1)=1:
D1Id3=transpose matrix {apply(flatten entries qbar3,i->sub(i,{l_(5,1)=>1,l_(5,2)=>0}))};
P3=(D_2**D_2**D_2**D_2)*D1Id3;
netList S_(positions(flatten entries P3,j->j!=0))
length S_(positions(flatten entries P3,j->j!=0))--9
M=M|P3;
rank M --15

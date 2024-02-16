restart

K=frac(QQ[p_1,p_2,p_3,p_4]);
RG=K[l_(1,1)..l_(5,3)]
-- 12|34 with identity at the leaves
qbar=value get "F84_4leaves_tensor_id.txt";
-- 12|34 no evolution point (star tree with identity at the leaves)
qbarId=sub(qbar,{l_(5,1)=>1,l_(5,2)=>1,l_(5,3)=>1});

-- Basis of the model
D={0,sub(diagonalMatrix{1,0,0,0},RG),sub(diagonalMatrix{0,1,0,0},RG),sub(diagonalMatrix{0,0,1,1},RG)}
S=sort elements (set toList(1..4))^**4/splice/splice;

--Replace all 4 by 3
Ss=apply(S,i->new MutableList from i);
for j from 0 to 255 do for i from 0 to 3 do if (Ss_j)#i==4 then (Ss_j)#i=3
Ss=unique apply(Ss,i->toSequence i)
length Ss --81

--Points obtained by using all combinations of D1,D2 and D3 at the leaves
-- with a generic matrix at the interior edge
P=time toList apply(Ss,i->(D_(i_0)**D_(i_1)**D_(i_2)**D_(i_3))*qbar);
-- with the identity matrix at the interior edge (no evolution point)
PId=time toList apply(Ss,i->(D_(i_0)**D_(i_1)**D_(i_2)**D_(i_3))*qbarId);

--45 points (non-zero coordinates) obtained by using all combinations of D1,D2 and D3 at the leaves
L=toList apply(0..80,i->S_(positions(flatten entries P_i,j->j!=0)));
length select(L,i->i!={}) --45
netList select(L,i->i!={})
netList Ss_(positions(L,i->i!={}))

LId=toList apply(0..80,i->S_(positions(flatten entries PId_i,j->j!=0)));
length select(LId,i->i!={}) --45
netList select(LId,i->i!={})
netList Ss_(positions(LId,i->i!={}))

--tests
Ss_(positions(LId,i->i!={}))==Ss_(positions(L,i->i!={})) --true
select(LId,i->i!={})==select(L,i->i!={}) --false

--tests
PId_(position(Ss,i->i==(3,3,3,3)))
netList (flatten entries PId_(position(Ss,i->i==(3,3,3,3))))_(positions(S,i->i==(3,3,4,4)))
netList (flatten entries P_(position(Ss,i->i==(3,3,3,3))))_(positions(S,i->i==(3,3,4,4)))
PId_(position(Ss,i->i==(2,3,2,3)))

--Linearly independent points obtained with the star tree
points=PId_(positions(LId,i->i!={}));
length points

M=points_0;
for i from 1 to 44  do {M=M|points_i};
rank M

-- All 45 points are linearly independent (regardless of whether we use Id or a generic matrix at the edge)
-- because they all have different non-zero coordinates

----------------------------------------------------------------------------
-- Additional linearly independent points obtained with tree topology 12|34
----------------------------------------------------------------------------

D1Id=transpose matrix {apply(flatten entries qbar,i->sub(i,{l_(5,1)=>1,l_(5,2)=>0,l_(5,3)=>0}))};
D2Id=transpose matrix {apply(flatten entries qbar,i->sub(i,{l_(5,1)=>0,l_(5,2)=>1,l_(5,3)=>0}))};
D3Id=transpose matrix {apply(flatten entries qbar,i->sub(i,{l_(5,1)=>0,l_(5,2)=>0,l_(5,3)=>1}))};

--For coordinates of type (3,3,3,3)
P=(D_3**D_3**D_3**D_3)*qbarId;
netList (select(flatten entries P,i->i!=0))
S_(positions(flatten entries P,j->j!=0))
P1=(D_3**D_3**D_3**D_3)*D1Id
netList (select(flatten entries P1,i->i!=0))
S_(positions(flatten entries P1,j->j!=0))
P2=(D_3**D_3**D_3**D_3)*D2Id;
netList (select(flatten entries P2,i->i!=0))
S_(positions(flatten entries P2,j->j!=0))
P3=(D_3**D_3**D_3**D_3)*D3Id;
netList (select(flatten entries P3,i->i!=0))
S_(positions(flatten entries P3,j->j!=0))

M=M|P1|P2;
rank M --47

-- Adding P3 does not increase the rank because is a multiple of the star tree for (3,3,3,3)

--For coordinates of type (2,2,3,3)
Q1=(D_2**D_2**D_3**D_3)*D1Id;
Q2=(D_3**D_3**D_2**D_2)*D1Id;
M=M|Q1|Q2;
rank M --49

----------------------------------------------------------------------------
-- Additional linearly independent points obtained with tree topology 13|24
----------------------------------------------------------------------------
-- 13|24 with identity at the leaves
qbar2=value get "F84_4leaves_tensor_id_1324.txt";
P3=(D_3**D_3**D_3**D_3)*qbar2;
S_(positions(flatten entries P3,j->j!=0))
netList select(flatten entries P3,j->j!=0)

--This is for generic qbar2 (i.e. generic Lambda_5)
-- but we can do it for l_(5,1)=1:
D1Id2=transpose matrix {apply(flatten entries qbar2,i->sub(i,{l_(5,1)=>1,l_(5,2)=>0,l_(5,3)=>0}))};
P3=(D_3**D_3**D_3**D_3)*D1Id2;
S_(positions(flatten entries P3,j->j!=0))
netList select(flatten entries P3,j->j!=0)


Q3=(D_2**D_3**D_2**D_3)*D1Id2;
S_(positions(flatten entries Q3,j->j!=0)) -- (2,3,2,3),(2,4,2,4)
netList select(flatten entries Q3,j->j!=0)

--tests
Q3gen=(D_2**D_3**D_2**D_3)*qbar2;
S_(positions(flatten entries Q3gen,j->j!=0)) -- (2,3,2,3),(2,4,2,4)
rank(Q3gen|PId_(position(Ss,i->i==(2,3,2,3)))) --2: We could use for Lambda_5: generic, Id, D1

--tests
Q3gen=(D_2**D_2**D_3**D_3)*qbar2;
S_(positions(flatten entries Q3gen,j->j!=0)) -- (2,2,3,3),(2,2,4,4)
rank(Q3gen|PId_(position(Ss,i->i==(2,2,3,3)))) --1


Q4=(D_3**D_2**D_3**D_2)*D1Id2;
S_(positions(flatten entries Q4,j->j!=0)) --(3,2,3,2),(4,2,4,2)
netList select(flatten entries Q4,j->j!=0)

M=M|P3|Q3|Q4;
rank M --52

----------------------------------------------------------------------------
-- Additional linearly independent points obtained with tree topology 14|23
----------------------------------------------------------------------------
-- 14|23 with identity at the leaves
qbar3=value get "F84_4leaves_tensor_id_1423.txt";
P4=(D_3**D_3**D_3**D_3)*qbar3;
S_(positions(flatten entries P4,j->j!=0))
netList select(flatten entries P4,j->j!=0)

--This is for generic qbar2 (i.e. generic Lambda_5)
-- but we can do it for l_(5,1)=1:
D1Id3=transpose matrix {apply(flatten entries qbar3,i->sub(i,{l_(5,1)=>1,l_(5,2)=>0,l_(5,3)=>0}))};
P4=(D_3**D_3**D_3**D_3)*D1Id3;
S_(positions(flatten entries P4,j->j!=0))
netList select(flatten entries P4,j->j!=0)

Q5=(D_2**D_3**D_3**D_2)*D1Id3;
S_(positions(flatten entries Q5,j->j!=0)) -- (2,3,3,2),(2,4,4,2)
netList select(flatten entries Q5,j->j!=0)

--tests
Q5gen=(D_2**D_3**D_3**D_2)*qbar3;
S_(positions(flatten entries Q5gen,j->j!=0)) -- (2,3,3,2),(2,4,4,2)
rank(Q5gen|PId_(position(Ss,i->i==(2,3,3,2)))) --2: We could use for Lambda_5: generic, Id, D1
--tests
Q5gen=(D_2**D_2**D_3**D_3)*qbar3;
S_(positions(flatten entries Q5gen,j->j!=0)) -- (2,2,3,3),(2,2,3,3)
rank(Q5gen|PId_(position(Ss,i->i==(2,2,3,3)))) --1

Q6=(D_3**D_2**D_2**D_3)*D1Id3;
S_(positions(flatten entries Q6,j->j!=0)) --(3,2,2,3),(4,2,2,4)
netList select(flatten entries Q6,j->j!=0)

M=M|P4|Q5|Q6;
rank M --55



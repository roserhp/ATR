restart
--Ring definition for QUARTETS with identity at the leaves and inner egde matrix
-- on p1..p4 and EIGENVALUES l1..l4
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K
gens R

--Diagonalizing change of basis
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})

--No evolution point
M1=id_(R^4)
M2=M1
M3=M1
M4=M1
M5=M1
M6=M1
M7=M1

st=toList(1..4)
S=sort elements (set st)^**5/splice/splice/splice
netList S

--Tensor for configuration 12|34
qq=mutableMatrix(R,4^5,1)
for i to 1023 do (qq_(i,0)=sum sum flatten toList apply(st,k->apply(st,kk->apply(st,kkk->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))*M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))*M3_(position(st,l->l==kk),position(st,l->l==(S_i)_2))*M4_(position(st,l->l==kkk),position(st,l->l==(S_i)_3))*M5_(position(st,l->l==kkk),position(st,l->l==(S_i)_4))*M6_(position(st,l->l==k),position(st,l->l==kk))*M7_(position(st,l->l==kk),position(st,l->l==kkk)))))) 
qq=matrix qq;

--Tensor with change of basis
H5=(transpose H)**(transpose H)**(transpose H)**(transpose H)**(transpose H);
pbar0=time H5*qq;
-- used 0.5 seconds

nonZeroEntries=S_(positions(flatten entries pbar0,i->i!=0))

length nonZeroEntries

netList apply(nonZeroEntries,i->{i,pbar0_(position(S,j->j==i),0)})

pbar0_(position(S,j->j==(1,1,1,3,3)),0)
pbar0_(position(S,j->j==(1,1,3,3,3)),0)
pbar0_(position(S,j->j==(1,3,3,3,3)),0)
pbar0_(position(S,j->j==(3,3,3,3,3)),0)


pbar0_(position(S,j->j==(1,2,2,3,3)),0)
pbar0_(position(S,j->j==(1,2,3,3,3)),0)
pbar0_(position(S,j->j==(2,2,3,3,3)),0)
factor denominator pbar0_(position(S,j->j==(2,2,3,3,3)),0)
pbar0_(position(S,j->j==(2,2,2,3,3)),0)
factor denominator pbar0_(position(S,j->j==(2,2,2,3,3)),0)


pbar0_(position(S,j->j==(2,2,2,2,2)),0)
factor numerator pbar0_(position(S,j->j==(2,2,2,2,2)),0)
factor denominator pbar0_(position(S,j->j==(2,2,2,2,2)),0)


pbar0_(position(S,j->j==(3,3,3,4,4)),0)

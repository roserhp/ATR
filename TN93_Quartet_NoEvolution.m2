restart
--Ring definition for QUARTETS with identity at the leaves and inner egde matrix
-- on p1..p4 and EIGENVALUES l1..l4
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K
gens R

--Diagonalizing change of basis
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})

--sanity check
inverse (transpose H) --A in paper up to p1+p2+p3+p4=1

--Identity at the leaves
M=id_(R^4)
M1=id_(R^4)
M2=M1
M3=M1
M4=M1

st=toList(1..4)
S=sort elements (set st)^**4/splice/splice
netList S

--Tensor for configuration 12|34
qq=mutableMatrix(R,256,1)
for i to 255 do (
	qq_(i,0)=sum flatten toList apply(st,k->apply(st,kk->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))*M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))*M3_(position(st,l->l==kk),position(st,l->l==(S_i)_2))*M4_(position(st,l->l==kk),position(st,l->l==(S_i)_3))*M_(position(st,l->l==k),position(st,l->l==kk))))
) 
qq=matrix qq;
netList (flatten entries qq)

--Tensor with change of basis
H4=(transpose H)**(transpose H)**(transpose H)**(transpose H);
pbar0=time H4*qq;
-- used 0.234375 seconds

"TN93_4leaves_noEvol.txt" << toString pbar0 << endl << close


nonZeroEntries=S_(positions(flatten entries pbar0,i->i!=0))

netList apply(nonZeroEntries,i->{i,pbar0_(position(S,j->j==(2,3,2,3)),0)})


pbar0_(position(S,j->j==(1,3,3,3)),0)
pbar0_(position(S,j->j==(2,2,2,2)),0)
factor numerator pbar0_(position(S,j->j==(2,2,2,2)),0)
factor denominator pbar0_(position(S,j->j==(2,2,2,2)),0)

(p_1+p_2)^2-(p_1+p_2)*(p_3+p_4)+(p_3+p_4)^2==p_1^2+2*p_1*p_2+p_2^2-p_1*p_3-p_2*p_3+p_3^2-p_1*p_4-p_2*p_4+2*p_3*
      p_4+p_4^2
pbar0_(position(S,j->j==(2,2,3,3)),0)
factor denominator pbar0_(position(S,j->j==(2,2,3,3)),0)
pbar0_(position(S,j->j==(2,2,3,3)),0)
factor denominator pbar0_(position(S,j->j==(2,2,4,4)),0)


pbar0_(position(S,j->j==(3,3,3,3)),0)
pbar0_(position(S,j->j==(4,4,4,4)),0)

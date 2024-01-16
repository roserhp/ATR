--TN93: build matrix M with entries in p_1..p_4 and parameters b,c,d
restart
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[b,c,d];
m=matrix{{0,c,b,b},{c,0,b,b},{b,b,0,d},{b,b,d,0}}
Dpi=sub(diagonalMatrix{p_1,p_2,p_3,p_4},R);
mm=mutableMatrix(m*Dpi)
for i to 3 do mm_(i,i)=1-sum flatten entries (matrix mm)^{i}

M=matrix mm
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
D=(inverse H)*M*H

--TN93: build matrix Ml with entries in p_1..p_4 and eigenvalues l1..l4
restart
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_1,l_2,l_3,l_4];
D=diagonalMatrix(gens R)
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
Ml=H*D*(inverse H)
Ml_(0,0)

--------------------------------------------------------------
--------------------------------------------------------------
--------------------------------------------------------------

--TN93: TRIPOD (on the eigenvalues)
restart

--Ring declaration. Warning: for some operations it will be needed to consider a ring on QQ including parameters pi
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(3,4)];
gens R

--Matrices with identity at the leaves
M1=id_(R^4)
M2=M1
M3=M1

--Building tensor 
st={1,2,3,4}
S=sort elements (set st)^**3/splice
netList S

qq=mutableMatrix(R,64,1)
for i to 63 do (
	qq_(i,0)=sum apply({1,2,3,4},k->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))
	                           *M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))
		                   *M3_(position(st,l->l==k),position(st,l->l==(S_i)_2)))) 
qq=matrix qq;
--How to check specific entries?
qq_(position(S,i->i==(2,2,2)),0)
qq_(position(S,i->i==(2,2,3)),0)

--Double-checking base change matrices appearing in paper: H and A
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
A=inverse(transpose H) --not exactly the same matrix because we are using p1+p2+p3+p4=1
B=sub(A,K)*diagonalMatrix{p_1+p_2+p_3+p_4,p_1+p_2+p_3+p_4,1,1} --using p1+p2+p3+p4=1, matches paper

--Matrices with anything at the leaves

D1=diagonalMatrix toList(l_(1,1)..l_(1,4))
D2=diagonalMatrix toList(l_(2,1)..l_(2,4))
D3=diagonalMatrix toList(l_(3,1)..l_(3,4))

M1g=H*D1*(inverse H)
M2g=H*D2*(inverse H)
M3g=H*D3*(inverse H)

--Change of basis in the tensor
qbar=((transpose H)**(transpose H)**(transpose H))*qq
"Tripod_qbar.txt" << toString qbar<< endl << close

--From identity at the leaves to general case
pbar=(D1**D2**D3)*qbar

--Double-checking relation between identity at the leaves and general case
qbar_(position(S,i->i==(2,2,2)),0)*l_(1,2)*l_(2,2)*l_(3,2)==pbar_(position(S,i->i==(2,2,2)),0)

--Table for pi-inner product of standard basis and basis B
Dpi=diagonalMatrix{p_1,p_2,p_3,p_4}
E=id_(K^4)
g=mutableMatrix(K,4,4);
for i to 3 do (for j to 3 do (
	g_(i,j)=(E^{i}*inverse(Dpi)*B_{j})_(0,0)/(transpose(B_{j})*inverse(Dpi)*B_{j})_(0,0)
    ))
g
G=(matrix g)*diagonalMatrix{p_1+p_2+p_3+p_4,p_1+p_2+p_3+p_4,1,1} --using p1+p2+p3+p4=1
entries G==entries H --true
for j to 3 do print factor (transpose B_{j}*inverse(Dpi)*B_{j})_(0,0) -- factors of <u_j,u_j>
--adjust with factors, matches with table in the paper
H*diagonalMatrix{1,(p_1+p_2)*(p_3+p_4),p_3*p_4/(p_3+p_4),p_1*p_2/(p_1+p_2)} --REVISAR!!!


---VANISHING IDEAL FOR THE TRIPOD

--declare ring with p variables (19 nonvanishing tensor coordinates)
nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0))
PBAR=flatten entries pbar^(positions(flatten entries pbar,i->i!=0))
netList PBAR
varp=toList apply(nonZeroEntries,i->(symbol p)_i);
Rp=K[varp]; 
gens Rp
--Vanishing ideal as kernel of the map
f=map(R,Rp,PBAR);
f(p_(1,1,1))
f(p_(4,4,4))
I=time trim kernel f;
-- used 4.88333 seconds
betti I

"3leavesVanishingIdeal_FinalBasis.txt" << toString I << endl << close

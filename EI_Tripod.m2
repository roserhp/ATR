--EI: build matrix M with entries in p_1..p_4 and parameters b
restart
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[b];
m=matrix{{0,b,b,b},{b,0,b,b},{b,b,0,b},{b,b,b,0}}
Dpi=sub(diagonalMatrix{p_1,p_2,p_3,p_4},R);
mm=mutableMatrix(m*Dpi)
for i to 3 do mm_(i,i)=1-sum flatten entries (matrix mm)^{i}

M=matrix mm
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
D=(inverse H)*M*H
D_(2,2)==D_(1,1)
D_(3,3)==D_(1,1)

--EI: build matrix Ml with entries in p_1..p_4 and eigenvalues l1,l2
restart
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_1,l_2];
D=matrix{{l_1,0,0,0},{0,l_2,0,0},{0,0,l_2,0},{0,0,0,l_2}}
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
Ml=H*D*(inverse H)

--------------------------------------------------------------
--------------------------------------------------------------
--------------------------------------------------------------

--EI: TRIPOD (on the eigenvalues)
restart

--Ring declaration. Warning: for some operations it will be needed to consider a ring on QQ including parameters pi
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(3,2)];
gens R

--Retrieve tensor for no evolution point in the new basis
qbar=value get "Tripod_qbar.txt";

--Index set
st={1,2,3,4}
S=sort elements (set st)^**3/splice;

--Matrices with anything at the leaves
D1=matrix{{l_(1,1),0,0,0},{0,l_(1,2),0,0},{0,0,l_(1,2),0},{0,0,0,l_(1,2)}}
D2=matrix{{l_(2,1),0,0,0},{0,l_(2,2),0,0},{0,0,l_(2,2),0},{0,0,0,l_(2,2)}}
D3=matrix{{l_(3,1),0,0,0},{0,l_(3,2),0,0},{0,0,l_(3,2),0},{0,0,0,l_(3,2)}}

--From identity at the leaves to general case
pbar=(D1**D2**D3)*qbar;

--Double-checking relation between identity at the leaves and general case
qbar_(position(S,i->i==(2,2,2)),0)*l_(1,2)*l_(2,2)*l_(3,2)==pbar_(position(S,i->i==(2,2,2)),0)

---VANISHING IDEAL FOR THE TRIPOD

--declare ring with p variables (19 nonvanishing tensor coordinates)
nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0))
length nonZeroEntries
length unique nonZeroEntries --same number of unique nonzero entries as in TN93
PBAR=flatten entries pbar^(positions(flatten entries pbar,i->i!=0));
varp=toList apply(nonZeroEntries,i->(symbol p)_i);
Rp=K[varp]; 
gens Rp
--Vanishing ideal as kernel of the map
f=map(R,Rp,PBAR);
f(p_(1,1,1))
f(p_(4,4,4))
I=time trim kernel f;
-- used 0.34375 seconds
"EI_3leavesVanishingIdeal_FinalBasis.txt" << toString I << endl << close
betti I --15 eq: 14 linear,1 cubic 
dim I, codim I, degree I --(4,15,3) --> complete intersection
netList I_*
support I_14
apply(0..14,i->length terms I_i) --> binomial ideal
netList toList apply(0..14,i->support I_i)

-- 64-45 vanishing variables=19
-- 19-14 linear eq=5
-- 31-1 cubic=4

--EI has 3 parametres (one for each of the 3 transition matrices). Hence dim V=3, dim CV=4.

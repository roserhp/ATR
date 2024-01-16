--HKY: build matrix M with entries in p_1..p_4 and parameters b,c
restart
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[b,c];
m=matrix{{0,c,b,b},{c,0,b,b},{b,b,0,c},{b,b,c,0}}
Dpi=sub(diagonalMatrix{p_1,p_2,p_3,p_4},R);
mm=mutableMatrix(m*Dpi)
for i to 3 do mm_(i,i)=1-sum flatten entries (matrix mm)^{i}

M=matrix mm
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
D=(inverse H)*M*H
D_(3,3)==(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*D_(1,1)+(p_1+p_2)*D_(2,2))

--Is HKY multiplicatively closed?

restart
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[b_1,c_1,b_2,c_2];
m1=matrix{{0,c_1,b_1,b_1},{c_1,0,b_1,b_1},{b_1,b_1,0,c_1},{b_1,b_1,c_1,0}}
Dpi=sub(diagonalMatrix{p_1,p_2,p_3,p_4},R);
mm1=mutableMatrix(m1*Dpi)
for i to 3 do mm1_(i,i)=1-sum flatten entries (matrix mm1)^{i}
M1=matrix mm1

m2=matrix{{0,c_2,b_2,b_2},{c_2,0,b_2,b_2},{b_2,b_2,0,c_2},{b_2,b_2,c_2,0}}
mm2=mutableMatrix(m2*Dpi)
for i to 3 do mm2_(i,i)=1-sum flatten entries (matrix mm2)^{i}
M2=matrix mm2

M1*M2==M2*M1

R1=M1*(inverse Dpi)

M=M1*M2
M_(0,1)==M_(2,3)
M_(0,1)
M_(2,3)

R=M1*M2*(inverse Dpi)
R==transpose R
R_(0,1)==R_(2,3)
R_(0,1)
R_(2,3)

R_(0,1)*p_2==M_(0,1)
R_(2,3)*p_4==M_(2,3)

R_(0,2)==R_(0,3)
R_(0,2)==R_(1,2)
R_(0,2)==R_(1,3)
R_(0,2)


Dpi*M1==(transpose M1)*Dpi --careful: says false but it's true
flatten entries (Dpi*M1)==flatten entries ((transpose M1)*Dpi)

H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
(inverse H)*M1*H
(inverse H)*M2*H
(inverse H)*M*H

D1=(inverse H)*M1*H
D1_(3,3)==(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*D1_(1,1)+(p_1+p_2)*D1_(2,2))


D2=(inverse H)*M2*H
D2_(3,3)==(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*D2_(1,1)+(p_1+p_2)*D2_(2,2))

D=(inverse H)*M*H

D_(3,3)==(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*D_(1,1)+(p_1+p_2)*D_(2,2))

D_(1,1)==D1_(1,1)*D2_(1,1)

D==D1*D2
D==D2*D1

DD1=matrix{{1,0,0,0},{0,b_1-1,0,0},{0,0,-(p_1+p_2)*b_1-(p_3+p_4)*c_1+1,0},{0,0,0,-(p_3+p_4)*b_1-(p_1+p_2)*c_1+1}}
DD2=matrix{{1,0,0,0},{0,b_2-1,0,0},{0,0,-(p_1+p_2)*b_2-(p_3+p_4)*c_2+1,0},{0,0,0,-(p_3+p_4)*b_2-(p_1+p_2)*c_2+1}}

DD12=DD1*DD2

DD12_(3,3)==(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*DD12_(1,1)+(p_1+p_2)*DD12_(2,2))

sub(DD12_(3,3),{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})==sub((1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*DD12_(1,1)+(p_1+p_2)*DD12_(2,2)),{p_1=>1/4,p_2=>1/4,p_3=>1/4,p_4=>1/4})
--true
sub(DD12_(3,3),{p_1=>1/8,p_2=>3/8,p_3=>1/10,p_4=>4/10})==sub((1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*DD12_(1,1)+(p_1+p_2)*DD12_(2,2)),{p_1=>1/8,p_2=>3/8,p_3=>1/10,p_4=>4/10})
--true
sub(DD12_(3,3),{p_1=>1-p_2-p_3-p_4})==sub((1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*DD12_(1,1)+(p_1+p_2)*DD12_(2,2)),{p_1=>1-p_2-p_3-p_4})
--false


--HKY: build matrix Ml with entries in p_1..p_4 and eigenvalues l1..l3
restart
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_1,l_2,l_3];
D=matrix{{l_1,0,0,0},{0,l_2,0,0},{0,0,l_3,0},{0,0,0,(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*l_2+(p_1+p_2)*l_3)}}
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
Ml=H*D*(inverse H)

--------------------------------------------------------------
--------------------------------------------------------------
--------------------------------------------------------------

--HKY: TRIPOD (on the eigenvalues)
restart

--Ring declaration. Warning: for some operations it will be needed to consider a ring on QQ including parameters pi
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(3,3)];

--Retrieve tensor for no evolution point in the new basis
qbar=value get "Tripod_qbar.txt";

--Index set
st={1,2,3,4}
S=sort elements (set st)^**3/splice;

--Matrices with anything at the leaves
D1=matrix{{l_(1,1),0,0,0},{0,l_(1,2),0,0},{0,0,l_(1,3),0},{0,0,0,(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*l_(1,2)+(p_1+p_2)*l_(1,3))}}
D2=matrix{{l_(2,1),0,0,0},{0,l_(2,2),0,0},{0,0,l_(2,3),0},{0,0,0,(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*l_(2,2)+(p_1+p_2)*l_(2,3))}}
D3=matrix{{l_(3,1),0,0,0},{0,l_(3,2),0,0},{0,0,l_(3,3),0},{0,0,0,(1/(p_3+p_4))*((p_3+p_4-p_1-p_2)*l_(3,2)+(p_1+p_2)*l_(3,3))}}

--From identity at the leaves to general case
pbar=(D1**D2**D3)*qbar;

--Double-checking relation between identity at the leaves and general case
qbar_(position(S,i->i==(2,2,2)),0)*l_(1,2)*l_(2,2)*l_(3,2)==pbar_(position(S,i->i==(2,2,2)),0)
qbar_(position(S,i->i==(3,3,3)),0)*l_(1,3)*l_(2,3)*l_(3,3)==pbar_(position(S,i->i==(3,3,3)),0)
qbar_(position(S,i->i==(3,3,3)),0)
pbar_(position(S,i->i==(3,3,3)),0)
qbar_(position(S,i->i==(4,4,4)),0)*D1_(3,3)*D2_(3,3)*D3_(3,3)==pbar_(position(S,i->i==(4,4,4)),0)
qbar_(position(S,i->i==(4,4,4)),0)
pbar_(position(S,i->i==(4,4,4)),0)
netList terms pbar_(position(S,i->i==(4,4,4)),0)
support pbar_(position(S,i->i==(4,4,4)),0)

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
-- used 847.062 seconds
"HKY_3leavesVanishingIdeal_FinalBasis.txt" << toString I << endl << close

-- Fast track to study the ideal
nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0))
varp=toList apply(nonZeroEntries,i->(symbol p)_i);
Rp=K[varp]; 
I=value get "HKY_3leavesVanishingIdeal_FinalBasis.txt";
betti I
I_0
time codim I
select(0..57,i-> degree I_i=={1})
select(0..57,i-> degree I_i=={2})
select(0..57,i-> degree I_i=={3})
select(0..57,i-> degree I_i=={4})
support I_0
length support I_0



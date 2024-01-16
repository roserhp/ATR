restart
--Ring definition for QUARTETS with identity at the leaves and inner egde matrix
-- on p1..p4 and EIGENVALUES l1..l2
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(5,1),l_(5,2)]
gens R

--Diagonalizing change of basis
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})

M=H*diagonalMatrix(R,4,4,{l_(5,1),l_(5,2),l_(5,2),l_(5,2)})*inverse(H)

--Identity at the leaves
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
qbar=time H4*qq;
-- used 2.8125 seconds

"F81_4leaves_tensor_id.txt" << toString qbar << endl << close

--Tensor for configuration 13|24
qq2=mutableMatrix(R,256,1)
for i to 255 do (
	qq2_(i,0)=sum flatten toList apply(st,k->apply(st,kk->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))*M3_(position(st,l->l==k),position(st,l->l==(S_i)_2))*M2_(position(st,l->l==kk),position(st,l->l==(S_i)_1))*M4_(position(st,l->l==kk),position(st,l->l==(S_i)_3))*M_(position(st,l->l==k),position(st,l->l==kk))))
) 
qq2=matrix qq2;
qbar2=time H4*qq2;
-- used 2.84375 seconds
"F81_4leaves_tensor_id_1324.txt" << toString qbar2 << endl << close
--Tensor for configuration 14|23
qq3=mutableMatrix(R,256,1)
for i to 255 do (
	qq3_(i,0)=sum flatten toList apply(st,k->apply(st,kk->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))*M4_(position(st,l->l==k),position(st,l->l==(S_i)_3))*M2_(position(st,l->l==kk),position(st,l->l==(S_i)_1))*M3_(position(st,l->l==kk),position(st,l->l==(S_i)_2))*M_(position(st,l->l==k),position(st,l->l==kk))))
) 
qq3=matrix qq3;
qbar3=time H4*qq3;
-- used 2.8125 seconds
"F81_4leaves_tensor_id_1423.txt" << toString qbar3 << endl << close


---------------------------------
---------------------------------
---------------------------------
--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms qbar_(position(S,j->j==i),0)>1))
netList nonMonomial
length nonMonomial --9
netList apply(nonMonomial,i->support qbar_(position(S,j->j==i),0))

nonZeroEntries=S_(positions(flatten entries qbar,i->i!=0))
length nonZeroEntries 
monomialNonZeroEntries=select(nonZeroEntries,i->(length terms qbar_(position(S,j->j==i),0)==1))
length monomialNonZeroEntries
netList monomialNonZeroEntries
---------------------------------
---------------------------------
---------------------------------

--Flattening 12|34 for identity at the leaves
flattq=mutableMatrix(R,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=p_(st_(i%4))*M_(i%4,j%4));      
    );
flattq=matrix(flattq);
flattq
--Change of basis of the flattening
flattQ=time (transpose(H)**transpose(H))*flattq*(H**H); 
-- used 0.0976014 seconds
--Quasi-block form (reordering according to s in next paragraph of code)
blockQ=flattQ_{0,3,12,7,13,15,5,1,4,10,2,8,6,9,11,14}^{0,3,12,7,13,15,5,1,4,10,2,8,6,9,11,14}
rank blockQ_{1,2,3,4}--1
rank blockQ_{14,15}--0
rank blockQ_{10,11,12,13}--1
rank blockQ_{7,8}--1

--Ring definition for flattening 12|34 general
Rgeneral=K[l_(1,1)..l_(5,4)]
gens Rgeneral
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}
-- Moving from identity at the leaves to general case
flattQ=sub(blockQ,Rgeneral);

flattP=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*flattQ_(i,j)));
aux=time sub(flattP,flatten toList apply(1..5,i->{l_(i,4)=>l_(i,2),l_(i,3)=>l_(i,2)}));
RG=K[l_(1,1)..l_(5,2)]
flattP=sub(aux,RG);
--rk 0
rank flattP_{position(s,i->i==(3,4)),position(s,i->i==(4,3))}--0
--rk 1
rank flattP_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,3))}--1
rank flattP_{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(3,2)),position(s,i->i==(2,3))}--1
rank flattP_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}--1
--rk 2
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}--2
--rk 3
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}--3
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,3))}--3
--rk 4: all
--column generators:
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}--4
--note that it forms a diagonal submatrix with non-zero entries in the diagonal as long as pi is generic and all eigenvalues are non-zero
flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}^{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}


--Tensor in general case
qbar=sub(qbar,Rgeneral);
qbar_(position(S,j->j==(2,2,3,3)),0)
use Rgeneral
pbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*qbar_(position(S,j->j==i),0))};
auxpbar=time sub(pbar,flatten toList apply(1..5,i->{l_(i,4)=>l_(i,2),l_(i,3)=>l_(i,2)}));
pbar=sub(auxpbar,RG);
--test
pbar_(position(S,j->j==(2,2,3,3)),0)
pbar_(position(S,j->j==(2,2,4,4)),0)
pbar_(position(S,j->j==(2,2,3,3)),0)
pbar_(position(S,j->j==(2,2,4,4)),0)
"F81_4leaves_tensor.txt" << toString pbar << endl << close


--Linear invariants arising from tripods
--(1,3,3,3) vs (1,4,4,4): 4
(p_3^2*p_4^2/(p_4^2-p_3^2))*pbar_(position(S,j->j==(1,3,3,3)),0)==(p_1^2*p_2^2/(p_2^2-p_1^2))*pbar_(position(S,j->j==(1,4,4,4)),0)
--(1,1,3,3) vs (1,1,4,4): 6
(p_3*p_4/(p_3+p_4))*pbar_(position(S,j->j==(1,1,3,3)),0)==(p_1*p_2/(p_2+p_1))*pbar_(position(S,j->j==(1,1,4,4)),0)
--(1,2,3,3): 12
p_3*p_4*pbar_(position(S,j->j==(1,2,3,3)),0)==-p_1*p_2*pbar_(position(S,j->j==(1,2,4,4)),0)


--Linear invariants that do not arise from tripods
--(2,3,2,3) vs (2,4,2,4): 4
(p_3*p_4*(p_3+p_4))*pbar_(position(S,j->j==(2,3,2,3)),0)==(p_1*p_2*(p_1+p_2))*pbar_(position(S,j->j==(2,4,2,4)),0)
--(2,3,3,3) vs (2,4,4,4): 4
(p_3^2*p_4^2/(p_3-p_4))*pbar_(position(S,j->j==(2,3,3,3)),0)==(p_1^2*p_2^2/(p_2-p_1))*pbar_(position(S,j->j==(2,4,4,4)),0)


Rx=K[toList apply(S,i->(symbol x)_i)]
g=map(RG,Rx,flatten entries pbar);
g(x_(1,1,1,1))
g(x_(1,1,4,4))

--Quadratic invariants
B1T=flattP_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,2))}
rank B1T

R1T=B1T^(positions(apply(0..15,i->B1T^{i}),j->j!=0))
s_(positions(apply(0..15,i->B1T^{i}),j->j!=0))

          14        j
14     (14,14)   (14,j)
i       (i,14)    (i,j)

rk1T=flatten toList apply({(4,1),(2,4),(4,2)},
    j->apply({(4,1),(2,4),(4,2),(4,4)},i->x_(1,4,1,4)*x_(i|j)-x_((1,4)|j)*x_(i|(1,4))))
length rk1T --12

B1G=flattP_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}
rank B1G

R1G=B1G^(positions(apply(0..15,i->B1G^{i}),j->j!=0))
s_(positions(apply(0..15,i->B1G^{i}),j->j!=0))

R1G=flattP_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}^{position(s,i->i==(4,4)),position(s,i->i==(2,2)),position(s,i->i==(1,2)),position(s,i->i==(2,1)),position(s,i->i==(3,3))}

         12         21
12    (12,12)     (12,21)
i      (i,12)      (i,21)

rk1G=toList apply({(2,1),(3,3),(2,2),(4,4)},i->x_(1,2,1,2)*x_(i|(2,1))-x_(1,2,2,1)*x_(i|(1,2)))
netList apply(rk1G,i->support i)
length rk1G --4
--Double-check
apply(rk1G,i->g(i))
CIrk1G=trim(ideal rk1G+ideal{p_3*p_4*x_(3,3,2,1)+p_1*p_2*x_(4,4,2,1),p_3*p_4*x_(3,3,1,2)+p_1*p_2*x_(4,4,1,2)})
netList CIrk1G_*
betti CIrk1

CIrk1=CIrk1G+ideal rk1T;


--Cubic invariants
g((p_3*p_4/(p_3+p_4))*x_(3,3,1,1)-(p_1*p_2/(p_2+p_1))*x_(4,4,1,1))
g(p_3*p_4*x_(3,3,2,1)+p_1*p_2*x_(4,4,2,1))
g(x_(3,3,2,2))
g(x_(4,4,2,2))

B2G=flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}
rank B2G
R2G=B2G^(positions(apply(0..15,i->B2G^{i}),j->j!=0))
s_(positions(apply(0..15,i->B2G^{i}),j->j!=0))


        11     12     22
   
11     1111   1112   1122

12     1211   1212   1222

i      i11     i12   i22


rk2=toList apply({(2,1),(3,3),(2,2),(4,4)},i->det matrix{{x_(1,1,1,1),x_(1,1,1,2),x_(1,1,2,2)},
	{x_(1,2,1,1),x_(1,2,1,2),x_(1,2,2,2)},{x_(i|(1,1)),x_(i|(1,2)),x_(i|(2,2))}})
netList rk2
netList apply(rk2,i->support i)
-- x_(3,3,2,2) present, won't be able to get rid of the eq with 3's
apply(rk2,i->g(i))
length rk2 --4

CIrk2=trim(ideal rk2+ideal{(p_3*p_4/(p_3+p_4))*x_(3,3,1,1)-(p_1*p_2/(p_2+p_1))*x_(4,4,1,1),p_3*p_4*x_(3,3,2,1)+p_1*p_2*x_(4,4,2,1)});
netList CIrk2_*
betti CIrk2
codim CIrk2 --4

--Quartic invariants
B3T=flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}
rank B3T

R3T=B3T^(positions(apply(0..15,i->B3T^{i}),j->j!=0))
s_(positions(apply(0..15,i->B3T^{i}),j->j!=0))

        11     12     14    44
   
11     1111   1112   1114  1144

12     1211   1212   1214  1244

14     1411   1412   1414  1444

i      i11    i12    i14   i44

rk3T=toList apply({(4,1),(2,4),(4,2),(4,4),(2,2),(2,1),(3,3)},
    i->det matrix{{x_(1,1,1,1),x_(1,1,1,2),x_(1,1,1,4),x_(1,1,4,4)},
	          {x_(1,2,1,1),x_(1,2,1,2),x_(1,2,1,4),x_(1,2,4,4)},
	          {x_(1,4,1,1),x_(1,4,1,2),x_(1,4,1,4),x_(1,4,4,4)},
	          {x_(i|(1,1)),x_(i|(1,2)),x_(i|(1,4)),x_(i|(4,4))}})

netList apply(rk3T,i->support i)
-- x_(3,3,4,4) present, won't be able to get rid of the eq with 3's
-- x_(3,3,1,4) present in the support...how is it possible?
g(x_(3,3,1,4))
netList terms rk3T_6
sub(rk3T_6,{x_(3,3,1,4)=>0,x_(1,1,1,2)=>0,x_(1,1,1,4)=>0,x_(1,2,1,1)=>0,
	x_(1,2,1,4)=>0,x_(1,4,1,1)=>0,x_(1,4,1,2)=>0})
apply(rk3T,i->g(i))
length rk3T --7
CIrk3=trim(ideal rk3T
    +ideal{(p_3*p_4/(p_3+p_4))*x_(3,3,1,1)-(p_1*p_2/(p_2+p_1))*x_(4,4,1,1),
	p_3*p_4*x_(3,3,2,1)+p_1*p_2*x_(4,4,2,1)});

CIedge=trim(CIrk1+CIrk2+CIrk3);
betti CIedge --29: 3 linear, 15 quadrics, 4 cubics, 7 quartics 
numgens CIedge --29
netList CIedge_*

jacEdge=jacobian CIedge;
time rank jacEdge --29 
--here we don't really know if all 29-minors could be 0

--no evolution point
p0=flatten entries sub(sub(sub(qbar,R),apply(gens R,i->i=>1)),Rx);
jacEdgeId=sub(jacEdge,matrix{p0});
time rank jacEdgeId --29
 -- used 0.779021 seconds




-------------------------------------------------------------
--Invariants of quartets arising from invariants of tripods
-------------------------------------------------------------

orderedVars={(1,1,1),(2,2,2),(3,3,3),(4,4,4),(1,2,2),(2,1,2),(2,2,1),(1,3,3),(3,1,3),(3,3,1),(1,4,4),(4,1,4),(4,4,1),(2,3,3),(3,2,3),(3,3,2),(2,4,4),(4,2,4),(4,4,2)}
--netList orderedVars
--length orderedVars --19
varp=toList apply(orderedVars,i->(symbol p)_i);
Rtripod=K[varp]
I=value get "3leavesVanishingIdeal_FinalBasis.txt"

CI=ideal{I_0,I_2,I_6,I_1,I_3,I_7,I_27,I_23,I_13}
codim CI
netList CI_*

s4=apply(orderedVars,i->append(i,1))
s1=apply(orderedVars,i->sequence 1|i);

g4=map(Rx,Rtripod,toList apply(s4,i->x_i));
CI4=g4(CI);
netList CI4_*
apply(flatten entries gens CI4,i->f(i))

g1=map(Rx,Rtripod,toList apply(s1,i->x_i));
CI1=g1(CI);
netList CI1_*
apply(flatten entries gens CI1,i->f(i))

CItripod=trim(CI1+CI4); 
time codim CItripod --14 => not a complete intersection
-- used 174.3 seconds
time codim CI1 --9
 -- used 0.554924 seconds
time codim CI4 --9
  -- used 0.529536 seconds
betti CItripod --18
netList CItripod_*


--------------------------------------------------------------
-- 10 NOV 2023
---------------------------------------------------------------
-- CHANGE SO THAT l_(j,4)=l_(j,3) OR REMOVE
--------------------------------------------------------------

--------------------------------------------------------------------
--------------------------------------------------------------------
-------------------------------------------------------------------
-- Alternative computation of flattenings
--------------------------------------------------------------------
--------------------------------------------------------------------
-------------------------------------------------------------------

restart

K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(5,1),l_(5,2),l_(5,3),l_(5,4)]
qbar=value get "4leaves_tensor_id_FinalBasis.txt";

S=sort elements (set {1,2,3,4})^**4/splice/splice;
--Flattening 12|34
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}
flattQ=mutableMatrix(R,4^2,4^2);
for i to 15 do (for j to 15 do
    flattQ_(i,j)=qbar_(position(S,k->k==((s_i)_0,(s_i)_1,(s_j)_0,(s_j)_1)),0);
    );
flattQ=matrix flattQ;

Rgeneral=K[l_(1,1)..l_(5,4)]
qbar=sub(qbar,Rgeneral);
pbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*qbar_(position(S,j->j==i),0))};
flattQ=sub(flattQ,Rgeneral);
flattP=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*flattQ_(i,j)));

testpbar=value get "4leaves_tensor_FinalBasis.txt";
testflattP=mutableMatrix(Rgeneral,16,16);
for i to 15 do (for j to 15 do
    testflattP_(i,j)=testpbar_(position(S,k->k==((s_i)_0,(s_i)_1,(s_j)_0,(s_j)_1)),0);
    );
testflattP=matrix testflattP;

testflattP==flattP --true


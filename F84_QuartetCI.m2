restart
K=frac(QQ[p_1,p_2,p_3,p_4]);
RG=K[l_(1,1)..l_(5,3)]

--Retrieve tensor on our desired basis
pbar=value get "F84_4leaves_tensor.txt";
qbar=value get "F84_4leaves_tensor_id.txt";

--Tensor coordinate indexes
S=sort elements (set {1,2,3,4})^**4/splice/splice;
--Column and row index in flattening matrix
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}

--Build flattening matrix
flattP=mutableMatrix(RG,16,16);
for i to 15 do (for j to 15 do
    flattP_(i,j)=pbar_(position(S,k->k==((s_i)_0,(s_i)_1,(s_j)_0,(s_j)_1)),0);
    );
flattP=matrix flattP;

-------------------------------------------------------------
--Invariants of quartets: setup
-------------------------------------------------------------
Rx=K[toList apply(S,i->(symbol x)_i)]
g=map(RG,Rx,flatten entries pbar);

-------------------------------------------------------------
--Edge invariants
-------------------------------------------------------------
--Quadratic invariants
B1T=flattP_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,2))}
R1T=B1T^(positions(apply(0..15,i->B1T^{i}),j->j!=0))
rk1T=flatten toList apply({(4,1),(2,4),(4,2)},
    j->apply({(4,1),(2,4),(4,2),(4,4)},i->x_(1,4,1,4)*x_(i|j)-x_((1,4)|j)*x_(i|(1,4))))
length rk1T --12

B1G=flattP_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}
R1G=flattP_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}^{position(s,i->i==(4,4)),position(s,i->i==(2,2)),position(s,i->i==(1,2)),position(s,i->i==(2,1)),position(s,i->i==(3,3))}
rk1G=toList apply({(2,1),(3,3),(2,2),(4,4)},i->x_(1,2,1,2)*x_(i|(2,1))-x_(1,2,2,1)*x_(i|(1,2)))
length rk1G --4

CIrk1=trim(ideal rk1G+ideal rk1T);
betti CIrk1

--Cubic invariants
B2G=flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}
R2G=B2G^(positions(apply(0..15,i->B2G^{i}),j->j!=0))
rk2=toList apply({(2,1),(3,3),(2,2),(4,4)},i->det matrix{{x_(1,1,1,1),x_(1,1,1,2),x_(1,1,2,2)},
	{x_(1,2,1,1),x_(1,2,1,2),x_(1,2,2,2)},{x_(i|(1,1)),x_(i|(1,2)),x_(i|(2,2))}})
length rk2 --4

CIrk2=trim(ideal rk2);
betti CIrk2

--Quartic invariants
B3T=flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}
R3T=B3T^(positions(apply(0..15,i->B3T^{i}),j->j!=0))
rk3T=toList apply({(4,1),(2,4),(4,2),(4,4),(2,2),(2,1),(3,3)},
    i->det matrix{{x_(1,1,1,1),x_(1,1,1,2),x_(1,1,1,4),x_(1,1,4,4)},
	          {x_(1,2,1,1),x_(1,2,1,2),x_(1,2,1,4),x_(1,2,4,4)},
	          {x_(1,4,1,1),x_(1,4,1,2),x_(1,4,1,4),x_(1,4,4,4)},
	          {x_(i|(1,1)),x_(i|(1,2)),x_(i|(1,4)),x_(i|(4,4))}})
length rk3T --7
CIrk3=trim(ideal rk3T);
betti CIrk3

--Linear invariants
rk0={x_(3,4,3,4),x_(3,4,4,3),x_(4,3,3,4),x_(4,3,4,3)}

CIEdge=trim((ideal rk0)+CIrk1+CIrk2+CIrk3);
betti CIEdge --1 of the 16 quadrics will disappear because of the additional linear invariants
             -- includes the 4 linear topology equations from TN93
-------------------------------------------------------------
-- Non-linear (model) invariants of quartets arising from invariants of tripods
-------------------------------------------------------------

orderedVars={(1,1,1),(2,2,2),(3,3,3),(4,4,4),(1,2,2),(2,1,2),(2,2,1),(1,3,3),(3,1,3),(3,3,1),(1,4,4),(4,1,4),(4,4,1),(2,3,3),(3,2,3),(3,3,2),(2,4,4),(4,2,4),(4,4,2)}
--netList orderedVars
--length orderedVars --19
varp=toList apply(orderedVars,i->(symbol p)_i);
Rtripod=K[varp]
I=value get "F84_3leavesVanishingIdeal_FinalBasis.txt"

CITripod=ideal{I_7,I_8,I_9,I_10,I_14};
dim CITripod, codim CITripod, degree CITripod --(14,5,72), 72=2^3*3^2 

s4=apply(orderedVars,i->append(i,1))
s1=apply(orderedVars,i->sequence 1|i);

g4=map(Rx,Rtripod,toList apply(s4,i->x_i));
CI4=g4(CITripod);
netList CI4_*
apply(flatten entries gens CI4,i->g(i))

g1=map(Rx,Rtripod,toList apply(s1,i->x_i));
CI1=g1(CITripod);
netList CI1_*
apply(flatten entries gens CI1,i->g(i))

CItripod=trim(CI1+CI4); 
betti CItripod --10: 6 quadrics + 4 cubics
netList CItripod_*

-------------------------------------------------------------
-- Linear (model) invariants of quartets arising from invariants of tripods
-------------------------------------------------------------

CILinearTripod=ideal{I_0,I_1,I_2,I_3,I_4,I_5,I_6}
s3=apply(orderedVars,i->(i_0,i_1,1,i_2))
s2=apply(orderedVars,i->(i_0,1,i_1,i_2))
g3=map(Rx,Rtripod,toList apply(s3,i->x_i));
CI3L=g3(CILinearTripod);
g2=map(Rx,Rtripod,toList apply(s2,i->x_i));
CI2L=g2(CILinearTripod);
CI1L=g1(CILinearTripod);
CI4L=g4(CILinearTripod);
CITripodLinear=trim(CI1L+CI2L+CI3L+CI4L);
betti CITripodLinear

-------------------------------------------------------------
-- Linear invariants NOT arising from invariants of tripods
-------------------------------------------------------------
--check
g((p_3*p_4*(p_3+p_4))*x_(2,3,2,3)-(p_1*p_2*(p_2+p_1))*x_(2,4,2,4))
g((p_3*p_4*(p_3+p_4))*x_(2,3,3,2)-(p_1*p_2*(p_2+p_1))*x_(2,4,4,2))
g((p_3*p_4*(p_3+p_4))*x_(3,2,3,2)-(p_1*p_2*(p_2+p_1))*x_(4,2,4,2))
g((p_3*p_4*(p_3+p_4))*x_(3,2,2,3)-(p_1*p_2*(p_2+p_1))*x_(4,2,2,4))
g((p_3^2*p_4^2/(p_3-p_4))*x_(2,3,3,3)-(p_1^2*p_2^2/(p_2-p_1))*x_(2,4,4,4))
g((p_3^2*p_4^2/(p_3-p_4))*x_(3,2,3,3)-(p_1^2*p_2^2/(p_2-p_1))*x_(4,2,4,4))
g((p_3^2*p_4^2/(p_3-p_4))*x_(3,3,2,3)-(p_1^2*p_2^2/(p_2-p_1))*x_(4,4,2,4))
g((p_3^2*p_4^2/(p_3-p_4))*x_(3,3,3,2)-(p_1^2*p_2^2/(p_2-p_1))*x_(4,4,4,2))
g(x_(3,3,4,4)-x_(4,4,3,3))


CILinear=ideal{(p_3*p_4*(p_3+p_4))*x_(2,3,2,3)-(p_1*p_2*(p_2+p_1))*x_(2,4,2,4),
(p_3*p_4*(p_3+p_4))*x_(2,3,3,2)-(p_1*p_2*(p_2+p_1))*x_(2,4,4,2),
(p_3*p_4*(p_3+p_4))*x_(3,2,3,2)-(p_1*p_2*(p_2+p_1))*x_(4,2,4,2),
(p_3*p_4*(p_3+p_4))*x_(3,2,2,3)-(p_1*p_2*(p_2+p_1))*x_(4,2,2,4),
(p_3^2*p_4^2/(p_3-p_4))*x_(2,3,3,3)-(p_1^2*p_2^2/(p_2-p_1))*x_(2,4,4,4),
(p_3^2*p_4^2/(p_3-p_4))*x_(3,2,3,3)-(p_1^2*p_2^2/(p_2-p_1))*x_(4,2,4,4),
(p_3^2*p_4^2/(p_3-p_4))*x_(3,3,2,3)-(p_1^2*p_2^2/(p_2-p_1))*x_(4,4,2,4),
(p_3^2*p_4^2/(p_3-p_4))*x_(3,3,3,2)-(p_1^2*p_2^2/(p_2-p_1))*x_(4,4,4,2),
x_(3,3,4,4)-x_(4,4,3,3)}

betti CILinear --9

CI=trim(CIEdge+CILinear+CItripod+CITripodLinear);
betti CI


jac=jacobian CI;
time rank jac --71
--here we don't really know if all 71-minors could be 0

--Check in no evolution point
p0=flatten entries sub(sub(sub(qbar,RG),apply(gens RG,i->i=>1)),Rx);
netList p0
jacId=sub(jac,matrix{p0});
time rank jacId --71


--2 missing equations
g(x_(2,2,3,3)*x_(4,4,2,2)-x_(3,3,2,2)*x_(2,2,4,4))

aux=trim(CI+ideal{x_(2,2,3,3)*x_(4,4,2,2)-x_(3,3,2,2)*x_(2,2,4,4)});
betti aux

jacaux=jacobian aux;
time rank jacaux --72
jacIdaux=sub(jacaux,matrix{p0});
time rank jacIdaux --72


---FALTA UNA EQUACIÓ NO LINEAL!!!


g(x_(2,2,3,3)*x_(4,4,2,2)*x_(3,3,4,4)-x_(2,2,4,4)*x_(3,3,2,2)*x_(4,4,3,3))
--redundant pq x_(3,3,4,4)=x_(4,4,3,3)

g(x_(2,2,3,3)*x_(4,4,2,2)*x_(3,3,3,3)-x_(2,2,4,4)*x_(3,3,2,2)*x_(4,4,4,4))

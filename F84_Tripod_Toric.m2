restart

--vars in the order of Lemma 4.1 (parametrization of bar{q})
orderedVars={(1,1,1),(2,2,2),(4,4,4),(1,2,2),(2,1,2),(2,2,1),(1,4,4),(4,1,4),(4,4,1),(2,4,4),(4,2,4),(4,4,2)}
--netList orderedVars
--length orderedVars --19
varp=toList apply(orderedVars,i->(symbol p)_i);
R=QQ[varp]
--matrix for TN93
C2=matrix{{1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0},{0,1,0,0,0,1,1,0,0,0,0,0,0,1,0,0,1,0,0},{0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0},{0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1},{1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0},{0,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0},{0,0,1,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0},{0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1},{1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0},{0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1},{0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0},{0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0}};
C=C2_{0,1,3,4,5,6,10,11,12,16,17,18}
rank C

A=C^{0}||C^{1}||C^{2}+C^{3}||C^{4}||C^{5}||C^{6}+C^{7}||C^{8}||C^{9}||C^{10}+C^{11}
rank A --7

A2=A^{0,1,2,3,4,6,7}
rank A2 --7

A3=A^{0,1,2,3,5,7,8}
rank A3 --7

loadPackage "Quasidegrees"
IA=time toricIdeal(A,R);
  -- used 0.015625 seconds
betti (trim IA)
netList IA_*
dim IA, codim IA --(7,5)
IA2=time toricIdeal(A2,R);
IA==IA2 --true

B=syz A --Columns v=(v_+)-(v_-) of B form a maximal set of linearly independent vectors of ker_Z A over Q
loadPackage "Binomials"
JB=time latticeBasisIdeal(R,B); --Binomial ideal generated by p^{v+}-p{v-}
betti JB --2 quartics, 2 cubics, 1 sextic
netList JB_*
dim JB, codim JB --(8,4)

B2=syz A2
JB2=time latticeBasisIdeal(R,B2); --Binomial ideal generated by p^{v+}-p{v-}
JB==JB2 --false
betti JB2 --1 quartic, 3 cubics, 1 quartic
netList JB2_*
dim JB2, codim JB2 --(9,3)

B3=syz A3
JB3=time latticeBasisIdeal(R,B3); --Binomial ideal generated by p^{v+}-p{v-}
betti JB3 --2 quartics, 2 cubics, 1 sextic
netList JB2_*
dim JB2, codim JB2 --(9,3)

-- Pq surten dimensions differents als lattice basis ideals???
-- No són interseccions completes!!!


IsubA2=time toricIdeal(subA2,R);
-- error: expected objects in the same ring
-- The number of columns of the matrix doesn't match the number of variables of the ring!

--double-check about non CM:
needsPackage "Depth"
isRegularSequence flatten entries gens JB --false
--we should be able to find a mixed submatrix with strictly less rows than columns:
B^{0,4}_{2,6,8}

isRegularSequence flatten entries gens JB2 --false
B2^{2,7}_{1,5,7,8}

BCI=transpose matrix{{0,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,0,-1,0,0},
		    {0,1,0,0,0,-1,0,0,0,0,0,1,0,0,0,0,0,-1,0},
    		    {0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,0,0,0,-1},
                    {0,0,0,0,0,0,0,-1,0,0,1,0,0,1,0,0,-1,0,0},
		    {0,0,0,0,0,0,0,0,-1,0,0,1,0,0,1,0,0,-1,0},
		    {0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,1,0,0,-1},
		    {-1,0,-2,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0},
		    {-1,0,0,-2,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0},
    		    {0,-1,-2,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0}}

JBCI=time latticeBasisIdeal(R,BCI); 
betti JBCI
netList JBCI_*
isRegularSequence flatten entries gens JBCI --true
IBCI=time saturate(JBCI,product gens R);
IBCI==IA --true

netList minimalPrimes IBCI
minPrimes=minimalPrimes JBCI;
length minPrimes
netList minPrimes
minPrimes_131==IBCI --true

netList apply(minPrimes,i-> support i)

J=time saturate(JBCI,p_(1,1,1));
minPrimesJ=time minimalPrimes J;
length minPrimesJ --39
netList minPrimesJ


trim(minPrimesJ_0+ideal{p_(1,1,1)-1})

J2=time saturate(JBCI,p_(1,1,1)*p_(2,2,2));
minPrimesJ2=time minimalPrimes J2;
length minPrimesJ2
netList minPrimesJ2

J4=time saturate(JBCI,p_(1,1,1)*p_(2,2,2)*p_(3,3,3)*p_(4,4,4));
minPrimesJ4=time minimalPrimes J4;
length minPrimesJ4
netList minPrimesJ4

a1=trim(JBCI+ideal{p_(1,1,1)-1});
a2=trim(IA+ideal{p_(1,1,1)-1});
a1==a2
length minimalPrimes a2 --1
length minimalPrimes a1 --39

sum apply(gens R,i->saturate(JBCI,i))
betti trim oo36

aux=JBCI;

for i in gens R do aux=saturate(aux,i);
betti trim aux
aux==IBCI --true

--To apply Th 2.1 in Dickenstein et al, we need to prove that any mixed submatrix will have more (or equal)
--rows than columns

maxminB=time minors(9,B);
maxminB ---Greatest common divisor of maximal minors of B is 1
-- According to Eisenbud-Sturmfels "Binomial ideals": 
---1. J_B is always CI in the Laurent polynomial ring 
-- 2. since g=1, then I_A and J_B coincide in the Laurent polynomial ring. 

IB=time saturate(JB,product gens R);
-- used 4.92187 seconds

IB==IA --true
isSubset(JB,IA)--true

--Trying to understand the difference between basis ideal and lattice basis ideal...
AQ=promote(A,QQ)
BQ=syz AQ
JBQ=time latticeBasisIdeal(R,BQ); --problematic because of it being in QQ


--Laurent polynomial ring
S=degreesRing 19
describe S
f=map(S,R,gens S);
f(p_(1,1,1))
LJB=f(JB);
netList LJB_*
LIA=f(IA);
netList LIA_*
LJB==LIA
--error: cannot compute Groebner basis of ideal over a Laurent polynomial ring, ie. with Inverses=>true

--generators
sub(p_(1,1,1)*p_(2,2,2)*p_(1,1,2)-p_(1,2,2)*p_(2,1,2)*p_(1,1,1),R)
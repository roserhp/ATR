restart
needsPackage "MultigradedImplicitization";
--needsPackage "Polyhedra";
--Ring definitions: 
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,2)];
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
Rx=K[var];

--Retrieve tensor on our desired basis
use R;
pbar=value get "EI_4leaves_tensor.txt";

F = map(R,Rx,transpose pbar);

--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms pbar_(position(S,j->j==i),0)>1));

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0));
length nonZeroEntries

monomialNonZeroEntries=select(nonZeroEntries,i->(length terms pbar_(position(S,j->j==i),0)==1));


-- substitute in a random root distribution
pbar = sub(pbar, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3});
nonZeroPBar = delete(null, apply(flatten entries pbar, m -> if m != 0 then m));

--Parametrization: almost monomial map

S = QQ[l_(1,1)..l_(5,2)];
T = QQ[apply(nonZeroEntries, ind -> x_ind)];
f = map(S,T, apply(nonZeroPBar, i -> sub(i,S)));

--grading
D = maxGrading f;
D

--dimension of grading
rank(D) --5

--compute linear invariants
G1=time componentsOfKernel(1,f);
   -- used 62.3125 seconds
peek G1
L=flatten toList apply(keys G1,i->G1#i)
length L
"EI_4leaves_linear.txt" << toString L << endl << close


--compute the quadratics
G = time componentsOfKernel(2,f);

H = (flatten delete({},values(G))) / (g -> sub(g, source f))

--double check in the ideal
(unique (H / f)) == {0}


--stick in a file
fileName = "TN93_quartet_quadrics" << "";

for g in H do (
    fileName << g << endl
);

fileName << close;

------------------------------------------
------------------------------------------
restart
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,2)];
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
Rx=K[var];

--Retrieve tensor on our desired basis
use R;
pbar=value get "EI_4leaves_tensor.txt";

F = map(R,Rx,transpose pbar);

--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms pbar_(position(S,j->j==i),0)>1));

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0));
length nonZeroEntries

monomialNonZeroEntries=select(nonZeroEntries,i->(length terms pbar_(position(S,j->j==i),0)==1));


-- substitute in a random root distribution
pbar = sub(pbar, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3});
nonZeroPBar = delete(null, apply(flatten entries pbar, m -> if m != 0 then m));

--Parametrization: almost monomial map

S = QQ[l_(1,1)..l_(5,2)];
T = QQ[apply(nonZeroEntries, ind -> x_ind)];
f = map(S,T, apply(nonZeroPBar, i -> sub(i,S)));

--grading
D = maxGrading f;
D

--dimension of grading
rank(D) --5

--compute linear invariants
G1=time componentsOfKernel(1,f);
   -- used 62.3125 seconds
peek G1
L=flatten toList apply(keys G1,i->G1#i)
length L
"EI_4leaves_linear.txt" << toString L << endl << close


--compute the quadratics
G = time componentsOfKernel(2,f);

H = (flatten delete({},values(G))) / (g -> sub(g, source f))

--double check in the ideal
(unique (H / f)) == {0}


--stick in a file
fileName = "TN93_quartet_quadrics" << "";

for g in H do (
    fileName << g << endl
);

fileName << close;

------------------------------------------
------------------------------------------
restart
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,2)];
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
Rx=K[var];

--Retrieve tensor on our desired basis
use R;
pbar=value get "EI_4leaves_tensor.txt";

F = map(R,Rx,transpose pbar);

--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms pbar_(position(S,j->j==i),0)>1));

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0));
length nonZeroEntries

monomialNonZeroEntries=select(nonZeroEntries,i->(length terms pbar_(position(S,j->j==i),0)==1));


-- substitute in a random root distribution
pbar = sub(pbar, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3});
nonZeroPBar = delete(null, apply(flatten entries pbar, m -> if m != 0 then m));

--Parametrization: almost monomial map

S = QQ[l_(1,1)..l_(5,2)];
T = QQ[apply(nonZeroEntries, ind -> x_ind)];
f = map(S,T, apply(nonZeroPBar, i -> sub(i,S)));

--grading
D = maxGrading f;
D

--dimension of grading
rank(D) --5

--compute linear invariants
G1=time componentsOfKernel(1,f);
   -- used 62.3125 seconds
peek G1
L=flatten toList apply(keys G1,i->G1#i)
length L
"EI_4leaves_linear.txt" << toString L << endl << close


--compute the quadratics
G = time componentsOfKernel(2,f);

H = (flatten delete({},values(G))) / (g -> sub(g, source f))

--double check in the ideal
(unique (H / f)) == {0}


--stick in a file
fileName = "TN93_quartet_quadrics" << "";

for g in H do (
    fileName << g << endl
);

fileName << close;

------------------------------------------
------------------------------------------
restart
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,2)];
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
Rx=K[var];

--Retrieve tensor on our desired basis
use R;
pbar=value get "EI_4leaves_tensor.txt";

F = map(R,Rx,transpose pbar);

--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms pbar_(position(S,j->j==i),0)>1));

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0));
length nonZeroEntries

monomialNonZeroEntries=select(nonZeroEntries,i->(length terms pbar_(position(S,j->j==i),0)==1));


-- substitute in a random root distribution
pbar = sub(pbar, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3});
nonZeroPBar = delete(null, apply(flatten entries pbar, m -> if m != 0 then m));

--Parametrization: almost monomial map

S = QQ[l_(1,1)..l_(5,2)];
T = QQ[apply(nonZeroEntries, ind -> x_ind)];
f = map(S,T, apply(nonZeroPBar, i -> sub(i,S)));

--grading
D = maxGrading f;
D

--dimension of grading
rank(D) --5

--compute linear invariants
G1=time componentsOfKernel(1,f);
   -- used 62.3125 seconds
peek G1
L=flatten toList apply(keys G1,i->G1#i)
length L
"EI_4leaves_linear.txt" << toString L << endl << close


--compute the quadratics
G = time componentsOfKernel(2,f);

H = (flatten delete({},values(G))) / (g -> sub(g, source f))

--double check in the ideal
(unique (H / f)) == {0}


--stick in a file
fileName = "TN93_quartet_quadrics" << "";

for g in H do (
    fileName << g << endl
);

fileName << close;

------------------------------------------
------------------------------------------
restart
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
Rx=QQ[var];
L=value get "EI_4leaves_linear.txt";
IL=ideal L;
betti IL
256-176-67
netList apply(L,i->support i)
length select(L,i->length (terms i)==2)
length select(L,i->length (terms i)==3)
netList select(L,i->length (terms i)==3)

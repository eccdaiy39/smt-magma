z:=2^44+2^22-2^9+2;
x:=z div 7;
r:=Ceiling((z^6+37*z^3+343)/343);
t:=(z^4+16*z+7)div 7;
p:=Ceiling((z^8+5*z^7+7*z^6+37*z^5+188*z^4+259*z^3+343*z^2+1763*z+2401)/21);
f:=29*17628883*2609815993167423713*17092188708684442834305509;
s:=5444521764485075480836860980926216661018;
w:=4368577123740566675106566649785592356161010496818595952740538619772484275351318\
19717017676708913034748308;
Fp:=GF(p);
F3<u>:=ExtensionField<Fp,u|u^3+7>;
F6<v>:=ExtensionField<F3,v|v^2-u>;
F18<w>:=ExtensionField<F6,w|w^3-v>;
E:=EllipticCurve([Fp|0,3]);
/*twisted curve*/
Et:=EllipticCurve([F3|0,3/u]);
ord_E:=#E;
ord_Et:=#Et;
function VectorG1()
B:=RMatrixSpace(Integers(), 2,2)![r,0,-s,1];
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
if GCD(C[1]^2-C[1]*C[2]+C[2]^2, ord_E) eq r then
printf"The vector C for G1 testing on KSS18-P348 is given as\n";
     return C;
else
    V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
     until GCD(C[1]^2-C[1]*C[2]+C[2]^2, ord_E) eq r;
printf"The vector C for G1 testing on KSS18-P348 is given as\n";
 return C;
end if;
end function;

function VectorG2()
B:=RMatrixSpace(Integers(), 6,6)!0;
B[1][1]:=r;
for i:=2 to 6 do
    B[i][1]:=-p^(i-1);B[i][i]:=1;
end for;
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
P<x,m,n>:=PolynomialRing(Integers(),3);
I := ideal<P|x^2-m*x+n>;
R:=quo<P|I>;
b:=R!(C[1]+C[2]*x+C[3]*x^2+C[4]*x^3+C[5]*x^4+C[6]*x^5);
b:=Evaluate(b,2,t);
b:=Evaluate(b,3,p);
b:=Coefficients(b);
if GCD(b[2]^2+b[1]*b[2]*t+b[1]^2*p, ord_Et) eq r then
     printf"The vector C for G2 testing on KSS18-P348 is given as\n";
     return C;
else 
     V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
        b:=R!(C[1]+C[2]*x+C[3]*x^2+C[4]*x^3+C[5]*x^4+C[6]*x^5);
        b:=Evaluate(b,2,t);
        b:=Evaluate(b,3,p);
        b:=Coefficients(b);
     until GCD(b[2]^2+b[1]*b[2]*t+b[1]^2*p, ord_Et) eq r;
printf"The vector C for G2 testing on KSS18-P348 is given as\n";
     return C;
end if;
end function;
function VectorGT()
B:=RMatrixSpace(Integers(), 6,6)!0;
B[1][1]:=r;
for i:=2 to 6 do
    B[i][1]:=-p^(i-1);B[i][i]:=1;
end for;
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;

if GCD(C[1]+C[2]*p+C[3]*p^2+C[4]*p^3+C[5]*p^4+C[6]*p^5, p^6-p^3+1) eq r then
     printf"The vector C for G2 testing on KSS18-P348 is given as\n";
     return C;
else 
     V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
     until GCD(C[1]+C[2]*p+C[3]*p^2+C[4]*p^3+C[5]*p^4+C[6]*p^5+C[7]*p^6+C[8]*p^7, p^8+1) eq r;
     printf"The vector C for GT testing on KSS16-P330 is given as\n";
    return C;
end if;
end function;
VectorG1();
VectorG2();
VectorGT();


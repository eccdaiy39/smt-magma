z:=-2^77+2^50+2^33;
p:=(z-1)^2*(z^4-z^2+1) div 3+z;
t:=z+1;
r:=z^4-z^2+1;
s:=-z^2;
Fp:=GF(p);
F2<u>:=ExtensionField<Fp,u|u^2+1>;
F6<v>:=ExtensionField<F2,v|v^3-u-2>;
F12<w>:=ExtensionField<F6,w|w^2-v>;
E:=EllipticCurve([Fp|0,4]);
/*twisted curve*/
Et:=EllipticCurve([F2|0,4*(2+u)]);
ord_E:=#E;
ord_Et:=#Et;
function VectorG1()
B:=RMatrixSpace(Integers(), 2,2)![r,0,-s,1];
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
if GCD(C[1]^2-C[1]*C[2]+C[2]^2, ord_E) eq r then
printf"The vector C for G1 testing on BLS12-P461 is given as\n";
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
B:=RMatrixSpace(Integers(), 4,4)!0;
B[1][1]:=r;
for i:=2 to 4 do
    B[i][1]:=-p^(i-1);B[i][i]:=1;
end for;
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
P<x,m,n>:=PolynomialRing(Integers(),3);
I := ideal<P|x^2-m*x+n>;
R:=quo<P|I>;
b:=R!(C[1]+C[2]*x+C[3]*x^2+C[4]*x^3);
b:=Evaluate(b,2,t);
b:=Evaluate(b,3,p);
b:=Coefficients(b);
if GCD(b[2]^2+b[1]*b[2]*t+b[1]^2*p, ord_Et) eq r then
     printf"The vector C for G2 testing on BLS12-P461 is given as\n";
     return C;
else 
     V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
        b:=R!(C[1]+C[2]*x+C[3]*x^2+C[4]*x^3);
        b:=Evaluate(b,2,t);
        b:=Evaluate(b,3,p);
        b:=Coefficients(b);
     until GCD(b[2]^2+b[1]*b[2]*t+b[1]^2*p, ord_Et) eq r;
printf"The vector C for G2 testing on BLS12-P461 is given as\n";
     return C;
end if;
end function;

function VectorGT()
B:=RMatrixSpace(Integers(), 4,4)!0;
B[1][1]:=r;
for i:=2 to 4 do
    B[i][1]:=-p^(i-1);B[i][i]:=1;
end for;
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;

if GCD(C[1]+C[2]*p+C[3]*p^2+C[4]*p^3, p^4-p^2+1) eq r then
     printf"The vector C for GT testing on BLS12-P461 is given as\n";
     return C;
else 
     V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
     until GCD(C[1]+C[2]*p+C[3]*p^2+C[4]*p^3, p^4-p^2+1) eq r;
     printf"The vector C for GT testing on BLS12-P461 is given as\n";
    return C;
end if;
end function;
VectorG1();
VectorG2();
VectorGT();

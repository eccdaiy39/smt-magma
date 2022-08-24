z:=-2^34+2^27-2^23+2^20-2^11+1;
r:=(z^8+48*z^4+625) div 61250;
t:=(2*z^5+41*z+35)div 35;
p:=(z^10+2*z^9+5*z^8+48*z^6+152*z^5+240*z^4+625*z^2+2398*z+3125) div 980;
Fp:=GF(p);
F4<u>:=ExtensionField<Fp,u|u^4-2>;
F8<v>:=ExtensionField<F4,v|v^2-u>;
F16<w>:=ExtensionField<F8,w|w^2-v>;
E:=EllipticCurve([Fp|1,0]);
Et:=EllipticCurve([F4|1/u,0]);
ord_E:=#E;
ord_Et:=#Et;

w1:=13975266739012216332523612368319120102936332749475824347676261174209213011\
79434231476966992563652766;
s1:=11676130122426325465158468012623314663795550790646633235959793344757293222\
8306;
function VectorG1()
B:=RMatrixSpace(Integers(), 2,2)![r,0,-s1,1];
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
if GCD(C[1]^2+C[2]^2, ord_E) eq r then
printf"The vector C for G1 testing on KSS16-P330 is given as\n";
     return C;
else
    V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
     until GCD(C[1]^2+C[2]^2, ord_E) eq r;
printf"The vector C for G1 testing is given as\n";
 return C;
end if;
end function;

function VectorG2()
B:=RMatrixSpace(Integers(), 8,8)!0;
B[1][1]:=r;
for i:=2 to 8 do
    B[i][1]:=-p^(i-1);B[i][i]:=1;
end for;
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
P<x,m,n>:=PolynomialRing(Integers(),3);
I := ideal<P|x^2-m*x+n>;
R:=quo<P|I>;
b:=R!(C[1]+C[2]*x+C[3]*x^2+C[4]*x^3+C[5]*x^4+C[6]*x^5+C[7]*x^6+C[8]*x^7);
b:=Evaluate(b,2,t);
b:=Evaluate(b,3,p);
b:=Coefficients(b);
if GCD(b[2]^2+b[1]*b[2]*t+b[1]^2*p, ord_Et) eq r then
     printf"The vector C for G2 testing on KSS16-P330 is given as\n";
     return C;
else 
     V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
        b:=R!(C[1]+C[2]*x+C[3]*x^2+C[4]*x^3+C[5]*x^4+C[6]*x^5+C[7]*x^6+C[8]*x^7);
        b:=Evaluate(b,2,t);
        b:=Evaluate(b,3,p);
        b:=Coefficients(b);
     until GCD(b[2]^2+b[1]*b[2]*t+b[1]^2*p, ord_Et) eq r;
printf"The vector C for G2 testing on KSS16-P330 is given as\n";
     return C;
end if;
end function;
function VectorGT()
B:=RMatrixSpace(Integers(), 8,8)!0;
B[1][1]:=r;
for i:=2 to 8 do
    B[i][1]:=-p^(i-1);B[i][i]:=1;
end for;
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;

if GCD(C[1]+C[2]*p+C[3]*p^2+C[4]*p^3+C[5]*p^4+C[6]*p^5+C[7]*p^6+C[8]*p^7, p^8+1) eq r then
     printf"The vector C for G2 testing on KSS16-P330 is given as\n";
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



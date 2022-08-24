P<u>:=PolynomialRing(RationalField());
r:=P!CyclotomicPolynomial(78);
p:=1/3*(u+1)^2*(u^26-u^13+1)-u^27;
z:=-2224;
p:=1749234309176102157657582860550885176950582224007184238236721873530271444092780\
387026731606667;
r:=2143085360734996117913472445644488914854141302995428209978212956147876055499508\
01;
t:=-72424742659885778123097924206425573989051009199;
Fp:=GF(p);
F13:=ext<Fp|u^13+2>;
E:=EllipticCurve([Fp|0,-17]);
/*w is the scalar multiplication[z^13-1] of GLV*/
w:=1749233248214262447516620809521841187877100017707061238854829525973209048875911\
623336506581723;
s:=-z^13;
ord_E:=#E;
function VectorG1()
B:=RMatrixSpace(Integers(), 2,2)![r,0,-s,1];
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
if GCD(C[1]^2-C[1]*C[2]+C[2]^2, ord_E) eq r then
printf"The vector C for G1 testing on BW13-P310 is given as\n";
     return C;
else
    V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
     until GCD(C[1]^2-C[1]*C[2]+C[2]^2, ord_E) eq r;
printf"The vector C for G1 testing on BW13-P310 is given as\n";
 return C;
end if;
end function;

function VectorGT1()
B:=RMatrixSpace(Integers(), 12,12)!0;
B[1][1]:=r;
for i:=2 to 12 do
    B[i][1]:=-p^(i-1);B[i][i]:=1;
end for;
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
a0:=0;a1:=0;
for i:=1 to 12 do
    a0:=a0+C[i]*p^(i-1);
    a1:=a1+p^(i-1);
end for;
a1:=a1+p^12;
if GCD(a0, a1) eq r then
     printf"The vector C for GT testing on BW13-P310 is given as\n";
     return C;
else 
     V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
        a0:=0;
        for i:=1 to 12 do
            a0:=a0+C[i]*p^(i-1);
        a1:=a1+p^(i-1);
end for;
a1:=a1+p^12;
     until GCD(a0, a1) eq r;
     printf"The vector C for GT testing on BW13-P310 is given as\n";
    return C;
end if;
end function;

function VectorGT2()
a0:=0;a1:=0;
C:=[z^2,-z,1,0,0,0,0,0,0,0,0,0];
#C;
for i:=1 to 12 do
    a0:=a0+C[i]*p^(i-1);
    a1:=a1+p^(i-1);
end for;
a1:=a1+p^12;
if GCD(a0, a1) eq r then
     printf"The vector C for GT testing on BW13-P310 also can be selected as\n";
     return C;
else 
     return 0;

end if;
end function;
VectorG1();
VectorGT1();
VectorGT2();


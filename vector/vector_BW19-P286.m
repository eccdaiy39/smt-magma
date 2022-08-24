z:=-145;
p:=(z+1)^2*(z^38-z^19+1)div 3-z^39;
t:=-z^20+z+1;
r:=Evaluate(CyclotomicPolynomial(114),z);
Fp:=GF(p);
F19<u>:=ExtensionField<Fp,u|u^19+u+31>;
E:=EllipticCurve([Fp|0,31]);
/*compute the order of E(F_{p^13})*/
t2:=t^2-2*p;
t3:=t*t2-t*p;
t4:=t*t3-t2*p;
t5:=t*t4-t3*p;
t6:=t*t5-t4*p;
t7:=t*t6-t5*p;
t8:=t*t7-t6*p;
t9:=t*t8-t7*p;
t10:=t*t9-t8*p;
t11:=t*t10-t9*p;
t12:=t*t11-t10*p;
t13:=t*t12-t11*p;
t14:=t*t13-t12*p;
t15:=t*t14-t13*p;
t16:=t*t15-t14*p;
t17:=t*t16-t15*p;
t18:=t*t17-t16*p;
t19:=t*t18-t17*p;
ord_E:=p^19+1-t19;
/*W1 is the scalar multiplication[z^19-1] of GLV*/
w1:=9561856609697676887108507068813558963674902373779955480974873474191350327753040\
8892292;
s:=-z^19;
function VectorG1()
B:=RMatrixSpace(Integers(), 2,2)![r,0,-s,1];
L:= LatticeWithBasis(B);
C:=ShortestVector(L);
min:=Norm(C);max:=2*min;
if GCD(C[1]^2-C[1]*C[2]+C[2]^2, ord_E) eq r then
printf"The vector C for G1 testing on BW19-P286 is given as\n";
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
function VectorGT()
B:=RMatrixSpace(Integers(), 18,18)!0;
B[1][1]:=r;
for i:=2 to 18 do
    B[i][1]:=-p^(i-1);B[i][i]:=1;
end for;
L:= LatticeWithBasis(B);
C:=ShortestVector(L);

min:=Norm(C);max:=2*min;
a0:=0;a1:=0;
for i:=1 to 18 do
    a0:=a0+C[i]*p^(i-1);
    a1:=a1+p^(i-1);
end for;
a1:=a1+p^18;
if GCD(a0, a1) eq r then
     printf"The vector C for GT testing on BW19-P286 is given as\n";
     return C;
else 
     V:=ShortVectorsProcess(L, min, max);
     repeat
        C:=NextVector(V);
        a0:=0;
        for i:=1 to 12 do
            a0:=a0+C[i]*p^(i-1);
        end for;
     until GCD(a0, a1) eq r;
     printf"The vector C for GT testing on BW19-286 is given as\n";
    return C;
end if;
end function;
 VectorG1();
VectorGT();



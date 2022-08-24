/*The parameters of BLS12-P461*/
z:=-2^77+2^50+2^33;
p:=(z-1)^2*(z^4-z^2+1) div 3+z;
t:=z+1;
r:=z^4-z^2+1;
Fp:=GF(p);
F2<u>:=ExtensionField<Fp,u|u^2+1>;
F6<v>:=ExtensionField<F2,v|v^3-u-2>;
F12<w>:=ExtensionField<F6,w|w^2-v>;
E:=EllipticCurve([Fp|0,4]);
E12:=EllipticCurve([F12|0,4]);
Et:=EllipticCurve([F2|0,4*(2+u)]);
cof1:=(p+1-t) div r;
cof2:=#Et div r;
coft:=(p^4-p^2+1) div r;
w2:=w^2;w3:=w^3;w2_inv:=1/w^2;w3_inv:=1/w^3;
/*GLV endmorphism: (x,y)->(w1*x, y)*/
//Random(Fp)^((p-1) div 3);
w1:=7880400945708837625552799996195794517320839498224072696258089854634610041561273\
1200443912142466655860687269376557054;
function endmo(Q,i)
        R:=Q;
        for t:=1 to i do
        R:=E12![w2_inv*R[1],w3_inv*R[2],1];
        R:=E12![R[1]^p,R[2]^p,1];
        R:=Et![w2*R[1],w3*R[2],1];
        end for;
        return R;
end function;
//*********************************G1 memebrship testing***************************//
//********************[a0, a1]=[1, -z^2]***************//

repeat
Q:=cof1*Random(E);
until Q[3] eq 1;
time_begin:=Cputime();
for i:=1 to 100 do
R0:=-z^2*Q;
R1:=E![w1*Q[1], Q[2],1];
end for;
T:=Cputime(time_begin);
if R0[1] eq R1[1] and R0[2] eq R1[2] then 
printf"PASS! The point is a member of G1\n";
else 
printf"REJECT! The point is NOT a member of G1\n";
end if;  
printf"Timing of the G1 memerbship testing  is %o*10^-2\n",T;

//*********************************G2 memebrship testing***************************//
//********************C=[z,-1,0,0]***************//
repeat
Q:=cof2*Random(Et);
until Q[3] ne 0;
time_begin:=Cputime();
for i:=1 to 100 do
R0:=z*Q;
R1:=endmo(Q,1);
end for;
T:=Cputime(time_begin);
if R0[1] eq  R1[1] and R0[2] eq R1[2] then 
printf"PASS! The point is a member of G2\n";
else 
printf"REJECT! The point is NOT a member of G2\n";
end if;  
printf"Timing of the G2 memerbship testing  is %o*100^-2\n",T;

//****************GT memebrship testing*************//
//****************the vector C=[z,-1,0,0]*************//
repeat
a:=Random(F12);
a:=Frobenius(a,Fp,6)/a;
a:=Frobenius(a,Fp,2)*a;
a:=a^coft;
until a ne 1 and a ne 0;
time_begin:=Cputime();
for i:=1 to 100 do
r0:=a^(-z);
r1:=Frobenius(a,Fp,1);
end for;
T:=Cputime(time_begin);
if r0*r1 eq 1 then 
printf"PASS! The point is a member of GT\n";
else 
printf"REJECT! The point is NOT a member of GT\n";
end if;  
printf"Timing of the GT memerbship testing  is %o*100^-2\n",T;

/*The parameters of BN-P446*/
z:=2^110+2^36+1;
p:=36*z^4+36*z^3+24*z^2+6*z+1;
r:=36*z^4+36*z^3+18*z^2+6*z+1;
t:=6*z^2+1;
Fp:=GF(p);
F2<u>:=ExtensionField<Fp,u|u^2+1>;
F6<v>:=ExtensionField<F2,v|v^3-(u+5)>;
F12<w>:=ExtensionField<F6,w|w^2-v>;
E:=EllipticCurve([Fp|0,257]);
E12:=EllipticCurve([F12|0,257]);

Et:=EllipticCurve([F2|0,257*(u+5)]);
//cof2:=#Et div r;
cof2:=1022116956040697189835203046526938749956395084607296049022800981998028463615288\
44466809084151324502007903115563866166482562610972590189;
coft:=(p^4-p^2+1) div r;

coft:=1067829166080176086040317559589693357779743430571623688800978265918071126018893\
4942914323079009901448294961145721051567198125954207403343433020557422606446451\
6759397147817606766614113492768836156429115376422103117949135800735104875833285\
3663065432412037524572188589296311502796924631324383132758793098124934917899933\
1664419217531309449507946643489895560910055810894881547268271533401441202001658\
31324689;
function endmo(Q,i)
        R:=Q;
        for t:=1 to i do
        R:=E12![1/w^2*R[1],1/w^3*R[2],1];
        R:=E12![R[1]^p,R[2]^p,1];
        R:=Et![w^2*R[1],w^3*R[2],1];
        end for;
        return R;
end function;
//****************G2 memebrship testing*************//
//****************the vector C=[z+1,z,z,-2z]*************//
repeat
Q:=cof2*Random(Et);
until Q[3] ne 0;
time_begin:=Cputime();
for i:=1 to 100 do
R0:=z*Q;
R1:=R0+Q;
R2:=2*R0;
R2:=endmo(R2, 3);
R0:=R1+endmo(R0, 1)+endmo(R0, 2);
end for;
T:=Cputime(time_begin);
if R0[1] eq  R2[1] and R0[2] eq R2[2] then 
printf"PASS! The point is a member of G2\n";
else 
printf"REJECT! The point is NOT a member of G2\n";
end if;  
printf"Timing of the G2 memerbship testing  is %o*100^-2\n",T;

//****************GT memebrship testing*************//
//****************the vector C=[z+1,z,z,-2z]*************//
repeat
a:=Random(F12);
a:=Frobenius(a,Fp,6)/a;
a:=Frobenius(a,Fp,2)*a;
a:=a^coft;
until a ne 1 and a ne 0;
time_begin:=Cputime();
for i:=1 to 100 do
r0:=a^z;
r1:=r0*a;
r2:=r0^2;
r2:=Frobenius(r2,Fp,3);
r0:=r1*Frobenius(r0,Fp,1)*Frobenius(r0,Fp,2);
end for;
T:=Cputime(time_begin);
if r0 eq r2 then 
printf"PASS! The point is a member of GT\n";
else 
printf"REJECT! The point is NOT a member of GT\n";
end if;  
printf"Timing of the GT memerbship testing  is %o*100^-2\n",T;


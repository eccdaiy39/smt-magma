/*The parameters of KSS18-P348*/

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
E18:=EllipticCurve([F18|0,3]);
Et:=EllipticCurve([F3|0,3/u]);
//cof1:=#E div r;
//cof2:= #Et div r;
//coft:= (p^6-p^3+1) div r;
cof1:=5054924237165166565754854231;
cof2:=9647053017183846247407587346851852477723302252014897837190815463447679093245699\
8094365283840474439301860519260707034829445465602656195359034273574581783622598\
6671906536019818492436310427263288984126762401815333515288862503156572441242623;
coft:=8042937371980804667637112974984651020193128778876067186765866624315446330498240\
9902071249596397357877231739095133904803540431506110760234255703769177343234728\
4888601712915242927581102876884081638768470174128386357731434091062436064174974\
7757444556259785571358847946334396657512025108241993640995581116620161079618978\
9164219279085847098125941525115397573630612144566923240985927823359603000429189\
2583551400000451419901531323375872444073045523327785683192213363598436523069570\
26064263566261406864896020325281049918403999417040039997860157845037822412331;
w2:=w^2;w3:=w^3;w2_inv:=1/w^2;w3_inv:=1/w^3;
/*GLV endmorphism: (x,y)->(w1*x, y)*/
w1:=4368577123740566675106566649785592356161010496818595952740538619772484275351318\
19717017676708913034748308;
function endmo(Q,i)
        R:=Q;
        for t:=1 to i do
        R:=E18![w2*R[1],w3*R[2],1];
        R:=E18![R[1]^p,R[2]^p,1];
        R:=Et![w2_inv*R[1],w3_inv*R[2],1];
        end for;
        return R;
end function;
//*********************************G1 memebrship testing***************************//
//********************[a0, a1]=[(z/7)^3, -18a0-1]***************//
a0:=(z  div 7)^3;
repeat
Q:=cof1*Random(E);
until Q[3] eq 1;
time_begin:=Cputime();
for i:=1 to 100 do
R0:=a0*Q;
R1:=(18*R0+Q);
R1:=E![w1*R1[1], R1[2],1];
end for;
T:=Cputime(time_begin);
if R0[1] eq R1[1] and R0[2] eq R1[2] then 
printf"PASS! The point is a member of G1\n";
else 
printf"REJECT! The point is NOT a member of G1\n";
end if;  
printf"Timing of the G1 memerbship testing  is %o*10^-2\n",T;
//*********************************G2 memebrship testing***************************//
//********************C=[2z/7,1,0,z/7,0,0]***************//
repeat
Q:=cof2*Random(Et);
until Q[3] ne 0;
u:=z div 7;
time_begin:=Cputime();
for i:=1 to 100 do
R0:=u*Q;
R1:=2*R0;
R0:=endmo(Q,1)+endmo(R0,3);
end for;
T:=Cputime(time_begin);
if R0[1] eq  R1[1] and R0[2] eq -R1[2] then 
printf"PASS! The point is a member of G2\n";
else 
printf"REJECT! The point is NOT a member of G2\n";
end if;  
printf"Timing of the G2 memerbship testing  is %o*100^-2\n",T;
//*********************************G2 memebrship testing***************************//
//********************C=[2z/7,1,0,z/7,0,0]***************//
u:=z div 7;
repeat
a:=Random(F18);
a:=Frobenius(a,Fp,9)/a;
a:=Frobenius(a,Fp,3)*a;
a:=a^coft;
until a ne 1 and a ne 0;
time_begin:=Cputime();
for i:=1 to 100 do
r0:=a^u;
r1:=r0^2;
r0:=r1*Frobenius(a,Fp,1)*Frobenius(r0,Fp,3);
r2:=Frobenius(a,Fp,6)*a;
r3:=Frobenius(a,Fp,3);
end for;
T:=Cputime(time_begin);
if r0 eq 1 and r2 eq r3 then 
printf"PASS! The point is a member of GT\n";
else 
printf"REJECT! The point is NOT a member of GT\n";
end if;  
printf"Timing of the GT memerbship testing  is %o*100^-2\n",T;




z:=-145;
p:=(z+1)^2*(z^38-z^19+1)div 3-z^39;
t:=-z^20+z+1;
r:=Evaluate(CyclotomicPolynomial(114),z);
Fp:=GF(p);
F19<u>:=ExtensionField<Fp,u|u^19+u+31>;
E:=EllipticCurve([Fp|0,31]);
E19:=EllipticCurve([F19|0,31]);

/*compute the order of E(F_{p^19})*/
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
cof1:=(p+1-t) div r;
cof2:=(p^19+1-t19) div (r*(p+1-t));
coft:=(p^19-1) div ((p-1)*r);
/*w is the scalar multiplication[z^19-1] of GLV*/
w:=9561856609697676887108507068813558963674902373779955480974873474191350327753040\
8892292;
//*********************************G1 memebrship testing***************************//
//********************[a0, a1]=[-(z^10-z)*(z^6-z^3+1)*(z+1),a0*z-1]*********************//
a0:=-(z^10-z)*(z^6-z^3+1)*(z+1);
repeat
Q:=cof1*Random(E);
until Q[3] eq 1;
time_begin:=Cputime();
for i:=1 to 100 do
R1:=a0*Q;R2:=Q-z*R1;
R2:=[w*R2[1],R2[2],1];
end for;
T:=Cputime(time_begin);
if R1[1] eq R2[1] and R1[2] eq R2[2] then 
printf"PASS! The point is a member of G1\n";
else 
printf"REJECT! The point is NOT a member of G1\n";
end if;  
printf"Timing of the G1 memerbship testing  is %o*10^-2\n",T;
//****************G2 memebrship testing*************//

TR19:=function(Q)
   T:=Q;
   for i:=1 to 18 do
       T:=T+E19![Q[1]^(p^i),Q[2]^(p^i),1];
   end for;
return T;
end function;
repeat
Q:=cof2*Random(E19);
until Q[3] ne 0;
Q:=19*Q-TR19(Q);
time_begin:=Cputime();
for i:=1 to 100 do
R0:=E19![Frobenius(Q[1],Fp,1),Frobenius(Q[2],Fp,1),1];
R1:=R0+E19![Frobenius(Q[1],Fp,2),Frobenius(Q[2],Fp,2),1]+E19![Frobenius(Q[1],Fp,3),Frobenius(Q[2],Fp,3),1];
R2:=E19![Frobenius(R1[1],Fp,3),Frobenius(R1[2],Fp,3),1];
R3:=R1+R2;
R4:=E19![Frobenius(R3[1],Fp,6),Frobenius(R3[2],Fp,6),1];
R5:=R3+R4;
R6:=E19![Frobenius(R4[1],Fp,6),Frobenius(R4[2],Fp,6),1];
R7:=Q+R5+R6;
R1:=E19![w*Q[1], Q[2]];
R1:=-z*R1;
end for;
T:=Cputime(time_begin);
if R7[3] eq 0 and R0[1] eq R1[1] and R0[2] eq R1[2] then 
printf"PASS! The point is a member of G2\n";
else 
printf"REJECT! The point is NOT a member of G2\n";
end if;  
printf"Timing of the G2 memerbship testing  is %o*10^-2\n",T;
//****************GT memebrship testing******************//
//****************the vector C=[z^2,-z,1,0...0,]*************//
//repeat
//a:=Random(F19);
//a:=Frobenius(a,Fp,1)/a;
//a:=a^coft;
//until a eq 1;
/*Since the time of expontiation by coft is too long, we obtain an element in advance*/
a:=3570450943411368132427555922028421756517217176140528955317750572144861117533404\
    7333085*u^18 + 547067271442280242142363379064686382717806416137488001323792\
    41211441084545762836174045*u^17 + 26676765215742733596199169134745204158691\
    439644661415826125439934650175231150116666582*u^16 +
    248500267296583379545538678044948126717964079481612252784597162427877650739\
    43203382408*u^15 + 32118208961645511823612468109686903773554609090868889218\
    788104288150738570715478802874*u^14 + 1656319721322359454675636231006510345\
    5747776075955706625253463533591740737142064803943*u^13 +
    722670809254025005292765574242989604453389880317793793682235152871828228649\
    48567282449*u^12 + 76902058074416706251567667630273140627711894336162979271\
    21677996479389049463518947880*u^11 + 35168322758782509768795905174470584952\
    86250169244825326765153129251457613133416517630*u^10 +
    332554728049663330410660821934354614669735337292750055428803083824938378333\
    12725549293*u^9 + 164751342846749297304087720381907895745169509798554146631\
    29206205430549167344714458576*u^8 + 587736609941550047657065138492292966491\
    40254262962312469028190200772835427186575779214*u^7 +
    172099454574376014542680533673903118122082191295265877441241417622851848641\
    27616197119*u^6 + 932183072722109574106467407214254436477556314099508427312\
    83135032779382464211310342576*u^5 + 840362891239573555068475922692224446557\
    54894288801077185629220474795217359885918566611*u^4 +
    964684219565540315118395035433004398512353030761442127629537154563805056217\
    2699737230*u^3 + 7536815887792192090259543732342328758508422224701663436874\
    3999502064390979136418096604*u^2 + 9189968332602972360597743327727538625535\
    6589966024197264679278802882429672558669668211*u +
    230350756350786701532016867268570711653515617643627928566461623557061242874\
    20469059924;
T:=Cputime(time_begin);
time_begin:=Cputime();
for i:=1 to 100 do
r0:=Frobenius(a,Fp,1);r1:=r0*Frobenius(a,Fp,2)*Frobenius(a,Fp,3);
r2:=Frobenius(r1,Fp,3);r3:=r1*r2;
r4:=Frobenius(r3,Fp,6);r5:=r3*r4;
r6:=Frobenius(r4,Fp,6);
r7:=a*r5*r6;
r1:=a^(-z);
r2:=r1^(-z);
r1:=r2*Frobenius(r1,Fp,1)*Frobenius(a,Fp,2);
end for;
T:=Cputime(time_begin);
if r7 eq 1 and r1 eq 1 then 
printf"PASS! The point is a member of GT\n";
else 
printf"REJECT! The point is NOT a member of GT\n";
end if;  
printf"Timing of the GT memerbship testing  is %o*10^-2\n",T;


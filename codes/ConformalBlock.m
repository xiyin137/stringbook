(* ::Package:: *)

BeginPackage["ConformalBlock`"];

(*
Print["ConformalBlock Version 2.01"];
Print["6 January 2016"];
Print[ "\[Copyright] Ying-Hsuan Lin, 2016"];
*)

IP::usage="IP[{\!\(\*SubscriptBox[\(l\), \(1\)]\), \!\(\*SubscriptBox[\(l\), \(2\)]\), ..}, {\!\(\*SubscriptBox[\(q\), \(1\)]\), \!\(\*SubscriptBox[\(q\), \(2\)]\), ...}, {\!\(\*SubscriptBox[\(r\), \(1\)]\), \!\(\*SubscriptBox[\(r\), \(2\)]\), ...}, \!\(\*SubscriptBox[\(d\), \(l\)]\), d, \!\(\*SubscriptBox[\(d\), \(r\)]\), c] computes the inner product \
\[LeftAngleBracket]\!\(\*SubscriptBox[\(d\), \(l\)]\)| \!\(\*SubscriptBox[\(L\), SubscriptBox[\(l\), \(1\)]]\) \!\(\*SubscriptBox[\(L\), SubscriptBox[\(l\), \(2\)]]\) ... (\!\(\*SubscriptBox[\(L\), SubscriptBox[\(q\), \(1\)]]\)\[CenterDot] \!\(\*SubscriptBox[\(L\), SubscriptBox[\(q\), \(2\)]]\)\[CenterDot] ... \[CenterDot] \!\(\*SubscriptBox[\(O\), \(d\)]\)(x)) \!\(\*SubscriptBox[\(L\), SubscriptBox[\(r\), \(1\)]]\) \!\(\*SubscriptBox[\(L\), SubscriptBox[\(r\), \(2\)]]\) ... |\!\(\*SubscriptBox[\(d\), \(r\)]\)\[RightAngleBracket] with the normalization \[LeftAngleBracket]\!\(\*SubscriptBox[\(d\), \(l\)]\)| \!\(\*SubscriptBox[\(O\), \(d\)]\)(x) |\!\(\*SubscriptBox[\(d\), \(r\)]\)\[RightAngleBracket] = 1.";

beta::usage="beta[NN, \!\(\*SubscriptBox[\(d\), \(1\)]\), \!\(\*SubscriptBox[\(d\), \(2\)]\), d, c] computes the coefficients of the level NN descendants of \!\(\*SubscriptBox[\(O\), \(d\)]\) in the \!\(\*SubscriptBox[\(O\), SubscriptBox[\(d\), \(1\)]]\)\!\(\*SubscriptBox[\(O\), SubscriptBox[\(d\), \(2\)]]\) OPE, \
with the normalization <\!\(\*SubscriptBox[\(O\), SubscriptBox[\(d\), \(1\)]]\) \!\(\*SubscriptBox[\(O\), SubscriptBox[\(d\), \(2\)]]\) \!\(\*SubscriptBox[\(O\), \(d\)]\)> = 1, and returns the coefficients in a list.";

S4CBzCoef::usage="S4CBzCoef[NN, \!\(\*SubscriptBox[\(d\), \(1\)]\), \!\(\*SubscriptBox[\(d\), \(2\)]\), \!\(\*SubscriptBox[\(d\), \(3\)]\), \!\(\*SubscriptBox[\(d\), \(4\)]\), d, c] computes the NN-th order coefficient in the z expansion \
of the sphere four-point conformal block with \!\(\*SubscriptBox[\(z\), \(1\)]\)=z, \!\(\*SubscriptBox[\(z\), \(2\)]\)=0, \!\(\*SubscriptBox[\(z\), \(3\)]\)=1, \!\(\*SubscriptBox[\(z\), \(4\)]\)=\[Infinity] (the four legs are ordered counterclockwise).";

S4CBz::usage="S4CBz[NN, \!\(\*SubscriptBox[\(d\), \(1\)]\), \!\(\*SubscriptBox[\(d\), \(2\)]\), \!\(\*SubscriptBox[\(d\), \(3\)]\), \!\(\*SubscriptBox[\(d\), \(4\)]\), d, c, z] gives the sphere four-point conformal block with \!\(\*SubscriptBox[\(z\), \(1\)]\)=z, \!\(\*SubscriptBox[\(z\), \(2\)]\)=0, \!\(\*SubscriptBox[\(z\), \(3\)]\)=1, \!\(\*SubscriptBox[\(z\), \(4\)]\)=\[Infinity] \
(the four legs are ordered counterclockwise) as a z expansion to the NN-th order \!\(\*SuperscriptBox[\(z\), \(d - \*SubscriptBox[\(d\), \(1\)] - \*SubscriptBox[\(d\), \(2\)]\)]\)(1 + ... + # z^NN).";

S4CBq::usage="S4CBq[NN, \!\(\*SubscriptBox[\(d\), \(1\)]\), \!\(\*SubscriptBox[\(d\), \(2\)]\), \!\(\*SubscriptBox[\(d\), \(3\)]\), \!\(\*SubscriptBox[\(d\), \(4\)]\), d, c, z, q] gives the sphere four-point conformal block with \!\(\*SubscriptBox[\(z\), \(1\)]\)=z, \!\(\*SubscriptBox[\(z\), \(2\)]\)=0, \!\(\*SubscriptBox[\(z\), \(3\)]\)=1, \!\(\*SubscriptBox[\(z\), \(4\)]\)=\[Infinity] \
(the four legs are labeled in counterclockwise order) in the form of a prefactor multiplying a q expansion to the NN-th order \
(this form appears in the context of Zamolodchikov's recurrence formula).";

S4CBqH::usage="";

S4CBqRecur::usage="S4CBqRecur[NN,d1,d2,d3,d4,c] gives the sphere four-point conformal block with \!\(\*SubscriptBox[\(z\), \(1\)]\)=z, \!\(\*SubscriptBox[\(z\), \(2\)]\)=0, \!\(\*SubscriptBox[\(z\), \(3\)]\)=1, \!\(\*SubscriptBox[\(z\), \(4\)]\)=\[Infinity] \
(the four legs are labeled in counterclockwise order) in the form of a prefactor multiplying a q expansion to the NN-th order \
using Zamolodchikov's recurrence formula. A function F[d,x,q] is returned.
S4CBqRecur[NN,d1,d2,d3,d4,c,ndig] works with ndig-digit floating numbers.";

S4CBzRecur::usage="";

S4CBqRecurH::usage="";

S4CBqRecurHb::usage="";

T1CBCoef::usage="T1CBCoef[NN, \!\(\*SubscriptBox[\(d\), \(1\)]\), d, c] computes the NN-th order coefficient in the q expansion of the torus one-point \
conformal block where q is the torus modulus.";

T1CB::usage="T1CB[NN, \!\(\*SubscriptBox[\(d\), \(1\)]\), d, c, q] gives the torus one-point conformal block as a q expansion to the NN-th order \
where q is the torus modulus.";

T2tCBCoef::usage="T2tCBCoef[\!\(\*SubscriptBox[\(NN\), \(a\)]\), \!\(\*SubscriptBox[\(NN\), \(b\)]\), \!\(\*SubscriptBox[\(d\), \(1\)]\), \!\(\*SubscriptBox[\(d\), \(2\)]\), \!\(\*SubscriptBox[\(d\), \(a\)]\), \!\(\*SubscriptBox[\(d\), \(b\)]\), c] computes the (\!\(\*SubscriptBox[\(NN\), \(a\)]\), \!\(\*SubscriptBox[\(NN\), \(b\)]\))-th order coefficient of the torus two-point \
conformal block in the t-channel (symmetric in a and b) with \!\(\*SubscriptBox[\(z\), \(1\)]\)=z, \!\(\*SubscriptBox[\(z\), \(2\)]\)=1 \
(glue 0 and \[Infinity] together in the sphere four-point conformal block).";

T2tCB::usage="T2tCB[\!\(\*SubscriptBox[\(NN\), \(a\)]\), \!\(\*SubscriptBox[\(NN\), \(b\)]\), \!\(\*SubscriptBox[\(d\), \(1\)]\), \!\(\*SubscriptBox[\(d\), \(2\)]\), \!\(\*SubscriptBox[\(d\), \(a\)]\), \!\(\*SubscriptBox[\(d\), \(b\)]\), c, \!\(\*SubscriptBox[\(q\), \(a\)]\), \!\(\*SubscriptBox[\(q\), \(b\)]\)] gives the torus two-point conformal block in the t-channel \
(symmetric in a and b) with \!\(\*SubscriptBox[\(z\), \(1\)]\)=z, \!\(\*SubscriptBox[\(z\), \(2\)]\)=1 (glue 0 and \[Infinity] together in the sphere four-point conformal block), as a double expansion \
in \!\(\*SubscriptBox[\(q\), \(a\)]\) and \!\(\*SubscriptBox[\(q\), \(b\)]\) up to the (\!\(\*SubscriptBox[\(NN\), \(a\)]\), \!\(\*SubscriptBox[\(NN\), \(b\)]\))-th order.";

T2sCBCoef::usage="T2sCBCoef[\!\(\*SubscriptBox[\(NN\), \(a\)]\), \!\(\*SubscriptBox[\(NN\), \(b\)]\), \!\(\*SubscriptBox[\(d\), \(1\)]\), \!\(\*SubscriptBox[\(d\), \(2\)]\), \!\(\*SubscriptBox[\(d\), \(a\)]\), \!\(\*SubscriptBox[\(d\), \(b\)]\), c] computes the (\!\(\*SubscriptBox[\(NN\), \(a\)]\), \!\(\*SubscriptBox[\(NN\), \(b\)]\))-th order coefficient of the torus two-point \
conformal block in the s-channel where \!\(\*SubscriptBox[\(x\), \(1\)]\)=\!\(\*SuperscriptBox[\(\[Xi]\), \(1/2\)]\), \!\(\*SubscriptBox[\(x\), \(2\)]\)=\!\(\*SuperscriptBox[\(\[Xi]\), \(\(-1\)/2\)]\), \!\(\*SubscriptBox[\(q\), \(a\)]\)=\!\(\*SubscriptBox[\(x\), \(1\)]\)-\!\(\*SubscriptBox[\(x\), \(2\)]\), and \!\(\*SubscriptBox[\(q\), \(b\)]\) is the torus modulus.";

T2sCB::usage="T2sCB[\!\(\*SubscriptBox[\(NN\), \(a\)]\), \!\(\*SubscriptBox[\(NN\), \(b\)]\), \!\(\*SubscriptBox[\(d\), \(1\)]\), \!\(\*SubscriptBox[\(d\), \(2\)]\), \!\(\*SubscriptBox[\(d\), \(a\)]\), \!\(\*SubscriptBox[\(d\), \(b\)]\), c, \!\(\*SubscriptBox[\(q\), \(a\)]\), \!\(\*SubscriptBox[\(q\), \(b\)]\)] gives the torus two-point conformal block in the s-channel \
(symmetric in a and b) where \!\(\*SubscriptBox[\(x\), \(1\)]\)=\!\(\*SuperscriptBox[\(\[Xi]\), \(1/2\)]\), \!\(\*SubscriptBox[\(x\), \(2\)]\)=\!\(\*SuperscriptBox[\(\[Xi]\), \(\(-1\)/2\)]\), \!\(\*SubscriptBox[\(q\), \(a\)]\)=\!\(\*SubscriptBox[\(x\), \(1\)]\)-\!\(\*SubscriptBox[\(x\), \(2\)]\), and \!\(\*SubscriptBox[\(q\), \(b\)]\) is the torus modulus, as a double expansion \
in \!\(\*SubscriptBox[\(q\), \(a\)]\) and \!\(\*SubscriptBox[\(q\), \(b\)]\) up to the (\!\(\*SubscriptBox[\(NN\), \(a\)]\), \!\(\*SubscriptBox[\(NN\), \(b\)]\))-th order.";


Begin["`Private`"];

(**** COMPUTE INNER PRODUCT <dl|L_{l} (L_{q} dot O_d(x)) L_{r}|dr> USING VIRASORO ALGEBRA ****)
(* position of the operator *)
Protect[x];

(* normalize *)
IP[{},{},{},d__]:=1;

(* if the operator is a descendant, convert L_n acting on the operator to L's acting on the left and on the right *)
IP[{l___},{n_,p___},{r___},d__]:=IP[{l},{n,p},{r},d]=Module[{nl,nr,t1,t2},
nl=Total[{l}];
nr=-Total[{r}];
t1=Sum[(-1)^(m)x^m Binomial[1+n,m]IP[{l,-m+n},{p},{r},d],{m,0,nl+n}];
t2=Sum[(-1)^(m-n)x^(1+n-m) Binomial[1+n,m]IP[{l},{p},{m-1,r},d],{m,0,1+nr}];
t1+t2//Simplify
]

(* shorthand when the operator is a primary *)
IP[{l___},{r___},d__]:=IP[{l},{},{r},d]

(**** For L_n to the left of the operator ****)
(* commute L_n for n>0 pass a primary operator from left to right *)
IP[{l___,n_},{},{r___},dl_,d_,dr_,c_]:=(IP[{l,n},{},{r},dl,d,dr,c]=Module[{},
If[d==0&&(dl+Total[{l}]-dr+Total[{r}]==0),Return[IP[{l},{},{n,r},dl,d,dr,c]]];
IP[{l},{},{n,r},dl,d,dr,c]+x^n(n d+dl+Total[{l}]-dr+Total[{r}])IP[{l},{},{r},dl,d,dr,c]//Simplify
])/;n>0

(* commute L_n for n\[LessEqual]0 all the way to the left *)
IP[{l___,n_,p___},{q___},{r___},dl_,d_,dr_,c_]:=(IP[{l,n,p},{q},{r},dl,d,dr,c]=Module[{t,i,lis},
If[n==0,
Return[(Total[{l}]+dl)IP[{l,p},{q},{r},dl,d,dr,c]]
];
t=0;
If[-n<=Total[{l}],
For[i=1,i<=Length[{l}],i++,
If[{l}[[i]]==-n,
t+=(-2n(Sum[{l}[[j]],{j,1,i-1}]+dl)-c (n^3-n)/12)IP[Join[Delete[{l},i],{p}],{q},{r},dl,d,dr,c],
lis={l};
lis[[i]]=lis[[i]]+n;
t+=({l}[[i]]-n)IP[Join[lis,{p}],{q},{r},dl,d,dr,c]
]
]
];
t//Simplify
])/;n<=0

(**** For L_n to the right of the operator ****)
(* if the operator is a primary, we can take the adjoint of the inner product so that L_{-n} is to the left of the operator,
and reduce to the previous two cases *)
IP[{l___},{},{p___,n_,r___},dl_,d_,dr_,c_]:=IP[{l},{},{p,n,r},dl,d,dr,c]=IP[-Reverse[{p,n,r}],{},-Reverse[{l}],dr,d,dl,c]

(**** A BASIS OF DESCENDANTS AT LEVEL NN ****)
(* for a primary of generic weight (no null descendant) *)
RB[NN_Integer]:=RB[NN]=-IntegerPartitions[NN]
LB[NN_Integer]:=LB[NN]=-Reverse[RB[NN],2]
(* a non-null basis of descendants of a zero-weight primary *)
RBV[NN_Integer]:=RBV[NN]=-IntegerPartitions[NN,All,Range[2,NN]]
LBV[NN_Integer]:=LBV[NN]=-Reverse[RBV[NN],2]

(**** GENERALIZED GRAM MATRIX: Compute <descendants| descendants of O(x)|descendants> ****)
(* compute the inner product of descendants at fixed levels of the bra, the ket, or the operator *)
(* NN can either be a non-negative integer or an empty list *)
(* the difference between NN=0 and NN={} is that NN=0 returns a list with an extra level of nesting *)
Gramx[NNl_,NNo_,NNr_,dl_,d_,dr_,c_]:=Gramx[NNl,NNo,NNr,dl,d,dr,c]=Module[{lb,ob,rb,dim},
{lb,ob,rb}={LB,RB,RB};
If[dl==0,lb=LBV];
If[d==0,ob=RBV];
If[dr==0,rb=RBV];

(* determine the dimension of the list (level of nesting) *)
dim=Boole[IntegerQ[NNl]]+Boole[IntegerQ[NNo]]+Boole[IntegerQ[NNr]];

If[dim==3,
Return[Array[IP[lb[NNl][[#1]],ob[NNo][[#2]],rb[NNr][[#3]],dl,d,dr,c]&,{Length[lb[NNl]],Length[ob[NNo]],Length[rb[NNr]]}]]
];

If[dim==2,
If[!IntegerQ[NNl],
Return[Array[IP[NNl,ob[NNo][[#1]],rb[NNr][[#2]],dl,d,dr,c]&,{Length[ob[NNo]],Length[rb[NNr]]}]]
];
If[!IntegerQ[NNo],
Return[Array[IP[lb[NNl][[#1]],NNo,rb[NNr][[#2]],dl,d,dr,c]&,{Length[lb[NNl]],Length[rb[NNr]]}]]
];
If[!IntegerQ[NNr],
Return[Array[IP[lb[NNl][[#1]],ob[NNo][[#2]],NNr,dl,d,dr,c]&,{Length[lb[NNl]],Length[ob[NNo]]}]]
]
];

If[dim==1,
If[IntegerQ[NNl],
Return[Array[IP[lb[NNl][[#]],NNo,NNr,dl,d,dr,c]&,{Length[lb[NNl]]}]]
];
If[IntegerQ[NNo],
Return[Array[IP[NNl,ob[NNo][[#]],NNr,dl,d,dr,c]&,{Length[ob[NNo]]}]]
];
If[IntegerQ[NNr],
Return[Array[IP[NNl,NNo,rb[NNr][[#]],dl,d,dr,c]&,{Length[rb[NNr]]}]]
]
];

Return[1];
]

(* generalized Gram matrix with operator inserted at x=1 *)
Gram[NNl_,NNo_,NNr_,dl_,d_,dr_,c_]:=Gram[NNl,NNo,NNr,dl,d,dr,c]=Gramx[NNl,NNo,NNr,dl,d,dr,c]/.{x->1}

(* ordinary Gram matrix when the operator is a primary *)
Gram[NNl_,NNr_,dl_,d_,dr_,c_]:=Gram[NNl,{},NNr,dl,d,dr,c]

(* inverse Gram matrix *)
InverseGram[NN_,d_,c_]:=Inverse[Gram[NN,{},NN,d,0,d,c]]

(**** CONFORMAL BLOCKS ****)

beta[NN_,d1_,d2_,d_,c_]:=beta[NN,d1,d2,d,c]=Module[{IG,GR},
If[d==0&&NN==1,Return[0]];
IG=InverseGram[NN,d,c];
GR=Gram[NN,{},{},d,d1,d2,c];
IG.GR
]

S4CBzCoef[NN_Integer,d1_,d2_,d3_,d4_,d_,c_]:=S4CBzCoef[NN,d1,d2,d3,d4,d,c]=Module[{GL,IG,GR},
If[d==0&&NN==1,Return[0]];
GL=Gram[{},NN,d4,d3,d,c];
IG=InverseGram[NN,d,c];
GR=Gram[NN,{},d,d1,d2,c];
GL.IG.GR
]

S4CBz[NN_Integer,d1_,d2_,d3_,d4_,d_,c_,z_]:=S4CBz[NN,d1,d2,d3,d4,d,c]=z^(d-d1-d2)Sum[S4CBzCoef[n,d1,d2,d3,d4,d,c]z^n,{n,0,NN}]

(* invert EllipticNomeQ[x] as a series in q *)
XinQ[NN_Integer]:=Module[{qinx,as,a,xinq,eqn},
qinx=Series[EllipticNomeQ[x],{x,0,NN}]//Normal;
as=Array[a[#]&,{NN}];
xinq=Sum[a[i]q^i,{i,1,NN}];
eqn=Array[SeriesCoefficient[qinx/.{x->xinq},{q,0,#}]&,{NN}]-Join[{1},Array[0&,{NN-1}]];
xinq/.(Solve[eqn==0,as][[1]])
]

S4CBq[NN_Integer,d1_,d2_,d3_,d4_,d_,c_,z0_,q0_]:=S4CBq[NN,d1,d2,d3,d4,d,c,z0,q0]=Module[{a,cb,H},
xinq=XinQ[NN+1];
cb=Sum[S4CBzCoef[n,d1,d2,d3,d4,d,c]xinq^n,{n,0,NN}];
H=Series[cb*(xinq/16/q)^(d-(c-1)/24)/(1-xinq)^((c-1)/24-d1-d3)/EllipticTheta[3,0,q]^(3(c-1)/6-4(d1+d2+d3+d4)),{q,0,NN}]//Normal;
(16q)^(d-(c-1)/24)z^((c-1)/24-d1-d2)(1-z)^((c-1)/24-d1-d3)EllipticTheta[3,0,q]^(3(c-1)/6-4(d1+d2+d3+d4))H/.{z->z0,q->q0}
]

S4CBqH[NN_Integer,d1_,d2_,d3_,d4_,d_,c_,z0_,q0_]:=S4CBqH[NN,d1,d2,d3,d4,d,c,z0,q0]=Module[{a,cb,H},
xinq=XinQ[NN+1];
cb=Sum[S4CBzCoef[n,d1,d2,d3,d4,d,c]xinq^n,{n,0,NN}];
H=Series[cb*(xinq/16/q)^(d-(c-1)/24)/(1-xinq)^((c-1)/24-d1-d3)/EllipticTheta[3,0,q]^(3(c-1)/6-4(d1+d2+d3+d4)),{q,0,NN}]//Simplify//Normal;
H/.{z->z0,q->q0}
]

T1CBCoef[NN_Integer,d1_,d_,c_]:=T1CBCoef[NN,d1,d,c]=Module[{IG,G},
If[d==0&&NN==1,Return[0]];
IG=InverseGram[NN,d,c];
G=Gram[NN,NN,d,d1,d,c];
Tr[IG.G]
]

T1CB[NN_Integer,d1_,d_,c_,q_]:=T1CB[NN,d1,d,c,q]=q^(d-c/24)Sum[T1CBCoef[n,d1,d,c]q^n,{n,0,NN}]

T2tCBCoef[NNa_Integer,NNb_Integer,d1_,d2_,da_,db_,c_]:=T2tCBCoef[NNa,NNb,d1,d2,da,db,c]=Module[{IGa,Gab,IGb,Gba},
If[da==0&&NNa==1,Return[0]];
If[db==0&&NNb==1,Return[0]];
IGa=InverseGram[NNa,da,c];
Gab=Gram[NNa,NNb,da,d1,db,c];
IGb=InverseGram[NNb,db,c];
Gba=Gram[NNb,NNa,db,d2,da,c];
Tr[IGa.Gab.IGb.Gba]
]

T2tCB[NNa_Integer,NNb_Integer,d1_,d2_,da_,db_,c_,qa_,qb_]:=T2tCB[NNa,NNb,d1,d2,da,db,c,qa,qb]=
qa^(da-c/24)qb^(db-c/24)Sum[T2tCBCoef[na,nb,d1,d2,da,db,c]qa^na qb^nb,{na,0,NNa},{nb,0,NNb}]

(* torus 2-pt CB in the s-channel with x2=x *)
T2sCBx[NNa_Integer,NNb_Integer,d1_,d2_,da_,db_,c_]:=T2sCBx[NNa,NNb,d1,d2,da,db,c]=Module[{GO,IGa,IGb,GR,vec,i},
If[da==0&&NNa==1,Return[0]];
If[db==0&&NNb==1,Return[0]];
GO=Transpose[Gramx[NNb,NNa,NNb,db,da,db,c],{2,1,3}];
IGa=Inverse[Gram[NNa,NNa,{},da,da,0,c]];
IGb=InverseGram[NNb,db,c];
GR=Gram[NNa,{},{},da,d1,d2,c];
vec={};
For[i=1,i<=Length[GO],i++,
AppendTo[vec,Tr[IGb.GO[[i]]]]
];
vec.IGa.GR
]

T2sCBCoef[NNa_Integer,NNb_Integer,d1_,d2_,da_,db_,c_]:=T2sCBCoef[NNa,NNb,d1,d2,da,db,c]=Module[{KO,IGa,IGb,KR,vec,i,qa},
If[da==0&&NNa==1,Return[0]];
If[db==0&&NNb==1,Return[0]];
SeriesCoefficient[Sum[((T2sCBx[na,NNb,d1,d2,da,db,c]x^(-da))/.{x->(-qa+Sqrt[4+qa^2])/2})qa^na,{na,0,NNa}],{qa,0,NNa}]
]

T2sCB[NNa_Integer,NNb_Integer,d1_,d2_,da_,db_,c_,qa_,qb_]:=T2sCB[NNa,NNb,d1,d2,da,db,c,qa,qb]=
qb^(db-c/24)Sum[qb^nb Series[Sum[((T2sCBx[na,nb,d1,d2,da,db,c]x^(-da))/.{x->(-qa+Sqrt[4+qa^2])/2})qa^na,{na,0,NNa}],{qa,0,NNa}]//Normal,{nb,0,NNb}]

(**** COMPUTE SPHERE FOUR-POINT CONFORMAL BLOCK VIA ZAMOLODCHIKOV'S RECURRENCE RELATION ****)

(* computes the inverse matrix A(x) of M(x) as an order NN series in x *)
MatrixSeriesInverse[M_,x_,NN_]:=Module[{Mseries,Mcoef,Acoef,i},

Mseries=Series[M,{x,0,NN}];
Mcoef={};
Acoef={};

For[i=0,i<=NN,i++,
AppendTo[Mcoef,SeriesCoefficient[Mseries,{x,0,i}]];
If[i==0,
AppendTo[Acoef,Inverse[Mcoef[[1]]]],
AppendTo[Acoef,-(Sum[Acoef[[j]].Mcoef[[i-j+2]],{j,1,i}]).Acoef[[1]]]
]
];

Sum[x^i Acoef[[i+1]],{i,0,NN}]
];

(* returns the number of integer lattice points (m,n) satisfying m*n\[LessEqual]NN, and a mapping from (m,n) to its ordering *)
Ordmn[NN_]:=Module[{len,m,n,ord,ordmn},
ord=0*IdentityMatrix[NN];
len=1;
For[m=1,m<=NN,m++,
For[n=1,n<=NN/m,n++,
ord[[m,n]]=len;
len++;
]
];
len-=1;
ordmn[m_,n_]:=ord[[m,n]];
{len,ordmn}
];

(**** Zamolodchikov's recurrence ****)
(* analytic *)
S4CBqRecur[NN_Integer,d1_,d2_,d3_,d4_,c_]:=S4CBqRecur[NN,d1,d2,d3,d4,c,{}]

(* numerical *)
S4CBqRecur[NN_Integer,d1_,d2_,d3_,d4_,c_,ndig_]:=S4CBqRecur[NN,d1,d2,d3,d4,c,ndig]=
Module[{Np,Q,b,\[Lambda],\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,d,R,Rvals,dvals,Hmn,H,q,F},
(* work with floating numbers with ndig digits if ndig is specified, otherwise simplify the analytic expression *)
If[NumberQ[ndig],
Np[x_]:=SetPrecision[x,ndig],
Np[x_]:=Simplify[x];
];

Q=Sqrt[(c-1)/6]//Np;
b=(Q-Sqrt[-4+Q^2])/2//Np;
\[Lambda][d_]:=Sqrt[Q^2/4-d]//Np;
\[Lambda]1=\[Lambda][d1];
\[Lambda]2=\[Lambda][d2];
\[Lambda]3=\[Lambda][d3];
\[Lambda]4=\[Lambda][d4];

\[Lambda][m_,n_]:=1/2(m/b+n b)//Np;
d[m_,n_]:=Q^2/4-(m/b+n b)^2/4//Np;
R[m_,n_]:=2 Product[(\[Lambda]1+\[Lambda]2-\[Lambda][r,s])(\[Lambda]1-\[Lambda]2-\[Lambda][r,s])(\[Lambda]3+\[Lambda]4-\[Lambda][r,s])(\[Lambda]3-\[Lambda]4-\[Lambda][r,s]),{r,-m+1,m-1,2},{s,-n+1,n-1,2}]\
/Product[If[(k==0&&l==0)||(k==m&&l==n),1,\[Lambda][k,l]],{k,-m+1,m},{l,-n+1,n}]//Np;

(* prestore values of R[m,n] and d[m,n] for reuse *)
Module[{m,n},
Rvals=0*IdentityMatrix[NN];
dvals=0*IdentityMatrix[NN];
For[m=1,m<=NN,m++,
For[n=1,n<=NN/m,n++,
Rvals[[m,n]]=R[m,n];
dvals[[m,n]]=If[m==1&&n==1,0,d[m,n]];
]
]
];

(* solve for the degenerate basis Hmn *)
Module[{len,ordmn,ones,M,r,s,m,n,denom,hmn},
{len,ordmn}=Ordmn[NN];
ones=Array[1&,{len}];

M=0*IdentityMatrix[len];
For[r=1,r<=NN,r++,
For[s=1,s<=NN/r,s++,
For[m=1,m<=(NN-r s),m++,
For[n=1,n<=(NN-r s)/m,n++,
denom=(dvals[[r,s]]+r s-dvals[[m,n]]);
M[[ordmn[r,s],ordmn[m,n]]]=q^(m n)Rvals[[m,n]]/If[NumberQ[denom]&&denom==0,1,denom];
]
]
]
];

hmn=MatrixSeriesInverse[IdentityMatrix[len]-M,q,NN].ones;
Hmn[m_,n_]:=hmn[[ordmn[m,n]]]
];

(**** H for any d ****)

H[d_]:=Module[{h,m,n,denom,dp},
dp=d//Np;

h=1;
For[m=1,m<=NN,m++,
For[n=1,n<=NN/m,n++,
denom=dp-dvals[[m,n]];
h+=Series[q^(m n) Rvals[[m,n]],{q,0,NN}]Hmn[m,n]/If[denom===0,1,denom];
]
];

h//Simplify//Normal
];

F[d_,x_,q0_]:=(16q)^(d-Q^2/4) x^(-d1-d2)x^(Q^2/4)(1-x)^(Q^2/4-d1-d3)EllipticTheta[3,0,q]^(3Q^2-4(d1+d2+d3+d4))H[d]/.{q->q0};

F
]

S4CBzRecur[NN_Integer,d1_,d2_,d3_,d4_,d_,c_]:=S4CBzRecur[NN,d1,d2,d3,d4,d,c,{}]

(* numerical *)
S4CBzRecur[NN_Integer,d1_,d2_,d3_,d4_,d_,c_,ndig_]:=S4CBzRecur[NN,d1,d2,d3,d4,d,c,ndig]=Module[{F,G,H},
F=S4CBqRecur[NN,d1,d2,d3,d4,c,ndig];
G=x^(d-d1-d2)(Series[x^(d1+d2-d)F[d,x,q]/.{q->EllipticNomeQ[x]},{x,0,NN}]//Normal);
H[x0_]:=G/.{x->x0};

H
]

(* H *)
S4CBqRecurH[NN_Integer,d1_,d2_,d3_,d4_,c_,ndig_]:=S4CBqRecurH[NN,d1,d2,d3,d4,c,ndig]=
Module[{Np,Q,b,\[Lambda],\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,d,R,Rvals,dvals,Hmn,H,q,F},
(* work with floating numbers with ndig digits if ndig is specified, otherwise simplify the analytic expression *)
If[NumberQ[ndig],
Np[x_]:=SetPrecision[x,ndig],
Np[x_]:=Simplify[x];
];

Q=Sqrt[(c-1)/6]//Np;
b=(Q-Sqrt[-4+Q^2])/2//Np;
\[Lambda][d_]:=Sqrt[Q^2/4-d]//Np;
\[Lambda]1=\[Lambda][d1];
\[Lambda]2=\[Lambda][d2];
\[Lambda]3=\[Lambda][d3];
\[Lambda]4=\[Lambda][d4];

\[Lambda][m_,n_]:=1/2(m/b+n b)//Np;
d[m_,n_]:=Q^2/4-(m/b+n b)^2/4//Np;
R[m_,n_]:=2 Product[(\[Lambda]1+\[Lambda]2-\[Lambda][r,s])(\[Lambda]1-\[Lambda]2-\[Lambda][r,s])(\[Lambda]3+\[Lambda]4-\[Lambda][r,s])(\[Lambda]3-\[Lambda]4-\[Lambda][r,s]),{r,-m+1,m-1,2},{s,-n+1,n-1,2}]\
/Product[If[(k==0&&l==0)||(k==m&&l==n),1,\[Lambda][k,l]],{k,-m+1,m},{l,-n+1,n}]//Np;

(* prestore values of R[m,n] and d[m,n] for reuse *)
Module[{m,n},
Rvals=0*IdentityMatrix[NN];
dvals=0*IdentityMatrix[NN];
For[m=1,m<=NN,m++,
For[n=1,n<=NN/m,n++,
Rvals[[m,n]]=R[m,n];
dvals[[m,n]]=If[m==1&&n==1,0,d[m,n]];
]
]
];

(* solve for the degenerate basis Hmn *)
Module[{len,ordmn,ones,M,r,s,m,n,denom,hmn},
{len,ordmn}=Ordmn[NN];
ones=Array[1&,{len}];

M=0*IdentityMatrix[len];
For[r=1,r<=NN,r++,
For[s=1,s<=NN/r,s++,
For[m=1,m<=(NN-r s),m++,
For[n=1,n<=(NN-r s)/m,n++,
denom=(dvals[[r,s]]+r s-dvals[[m,n]]);
M[[ordmn[r,s],ordmn[m,n]]]=q^(m n)Rvals[[m,n]]/If[NumberQ[denom]&&denom==0,1,denom];
]
]
]
];

hmn=MatrixSeriesInverse[IdentityMatrix[len]-M,q,NN].ones;
Hmn[m_,n_]:=hmn[[ordmn[m,n]]]
];

(**** H for any d ****)

H[d_]:=Module[{h,m,n,denom,dp},
dp=d//Np;

h=1;
For[m=1,m<=NN,m++,
For[n=1,n<=NN/m,n++,
denom=dp-dvals[[m,n]];
h+=Series[q^(m n) Rvals[[m,n]],{q,0,NN}]Hmn[m,n]/If[denom===0,1,denom];
]
];

h//Simplify//Normal
];

F[d_,x_,q0_] := H[d]/.{q->q0};

F
]

S4CBqRecurHb[NN_Integer,d1_,d2_,d3_,d4_,b_,ndig_]:=S4CBqRecurHb[NN,d1,d2,d3,d4,c,ndig]=
Module[{Np,Q,\[Lambda],\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,d,R,Rvals,dvals,Hmn,H,q,F},
(* work with floating numbers with ndig digits if ndig is specified, otherwise simplify the analytic expression *)
If[NumberQ[ndig],
Np[x_]:=SetPrecision[x,ndig],
Np[x_]:=Simplify[x];
];

Q = b+1/b;
\[Lambda][d_]:=Simplify[Sqrt[Q^2/4-d]//Np, Assumptions->{b>0, b<1}];
\[Lambda]1=\[Lambda][d1];
\[Lambda]2=\[Lambda][d2];
\[Lambda]3=\[Lambda][d3];
\[Lambda]4=\[Lambda][d4];

\[Lambda][m_,n_]:=1/2(m/b+n b)//Np;
d[m_,n_]:=Q^2/4-(m/b+n b)^2/4//Np;
R[m_,n_]:=2 Product[(\[Lambda]1+\[Lambda]2-\[Lambda][r,s])(\[Lambda]1-\[Lambda]2-\[Lambda][r,s])(\[Lambda]3+\[Lambda]4-\[Lambda][r,s])(\[Lambda]3-\[Lambda]4-\[Lambda][r,s]),{r,-m+1,m-1,2},{s,-n+1,n-1,2}]\
/Product[If[(k==0&&l==0)||(k==m&&l==n),1,\[Lambda][k,l]],{k,-m+1,m},{l,-n+1,n}]//Np;

(* prestore values of R[m,n] and d[m,n] for reuse *)
Module[{m,n},
Rvals=0*IdentityMatrix[NN];
dvals=0*IdentityMatrix[NN];
For[m=1,m<=NN,m++,
For[n=1,n<=NN/m,n++,
Rvals[[m,n]]=R[m,n];
dvals[[m,n]]=If[m==1&&n==1,0,d[m,n]];
]
]
];

(* solve for the degenerate basis Hmn *)
Module[{len,ordmn,ones,M,r,s,m,n,denom,hmn},
{len,ordmn}=Ordmn[NN];
ones=Array[1&,{len}];

M=0*IdentityMatrix[len];
For[r=1,r<=NN,r++,
For[s=1,s<=NN/r,s++,
For[m=1,m<=(NN-r s),m++,
For[n=1,n<=(NN-r s)/m,n++,
denom=(dvals[[r,s]]+r s-dvals[[m,n]]);
M[[ordmn[r,s],ordmn[m,n]]]=q^(m n)Rvals[[m,n]]/If[NumberQ[denom]&&denom==0,1,denom];
]
]
]
];

hmn=MatrixSeriesInverse[IdentityMatrix[len]-M,q,NN].ones;
Hmn[m_,n_]:=hmn[[ordmn[m,n]]]
];

(**** H for any d ****)

H[d_]:=Module[{h,m,n,denom,dp},
dp=d//Np;

h=1;
For[m=1,m<=NN,m++,
For[n=1,n<=NN/m,n++,
denom=dp-dvals[[m,n]];
h+=Series[q^(m n) Rvals[[m,n]],{q,0,NN}]Hmn[m,n]/If[denom===0,1,denom];
]
];

h//Simplify//Normal
];

F[d_,x_,q0_] := H[d]/.{q->q0};

F
]

End[];

EndPackage[];

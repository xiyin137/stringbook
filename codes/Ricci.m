(* :Name: Ricci`Ricci` *)

(* :Title: Ricci *)

(* :Author: Juan M. Aguirregabiria *)

(* :Summary:
  Basic definitions for the calculation, in a coordinate frame, of the 
  contravariant metric tensor, Christoffel symbols, Riemann, Ricci, Einstein 
  and Weyl tensors and the scalar curvature, for a given metric in arbitrary 
  dimension. 
*)

(* :Context: Ricci`Ricci` *)
   
(* :Package Version: 1.44 *)

(* :Mathematica Version: 3.0, 4.0 *)

(* :Copyright: Copyright 1989-2002 Juan M. Aguirregabiria *)

(* :History:
  Version 1.44, July     2002. DetMetric.
  Version 1.43, February 2001. Christoffel is printed as capital gamma.
  Version 1.42, May      2000. Missprint in "Usage:"
  Version 1.41, October  1999. Bug in PrintRiemann
  Version 1.40, May      1998. SetMetric, local variables
  Version 1.34, May      1998. Optional formatting in Print*
                               PrintScalarCurvature[]
  Version 1.33, May      1997. PrintScalarCurvature
  Version 1.32, April    1997. Context names
  Version 1.31, December 1993. Bug in PrintEinstein
  Version 1.3,  December 1991. Weyl tensor and Print...
  Version 1.2,  June     1990. Derivative name changed!
  Version 1.1,  April    1990.
  Version 1.0,  December 1989.
*)

(* :Keywords:
  differential geometry, riemannian geometry, general relativity, 
  Riemann tensor, Ricci tensor, Einstein equations, coordinate frame, 
  affine connection, Christoffel symbols
*)

(* :Source:
  C.W. Misner, K.S. Thorne and J.A. Wheeler, "Gravitation", Freeman, 
  San Francisco (1973).
*)

(* :Warning:
  According to the standard mathematical notation for the Riemann 
  tensor, the contravariant index will be always the first one, 
  instead of being the last one as usual in Mathematica.
*)

(* :Discussion:

    Usage: 
  
  1. Needs["Ricci`Ricci`"].
  2. Introduce the appropriate auxiliary definitions, functions...
  3. Use SetMetric to introduce coordinates and metric or:
     3.1 Set FirstIndex, LastIndex, if their default
           values are not appropriate. 
     3.2 Introduce the name of each Coordinate: Coordinate[n] = ...
     3.3 Introduce the essential non-null components of the MetricTensor:
           MetricTensor[i,j] = ...
  4. Set RicciSimplify, if its default value is not appropriate. 
  5. Define the appropriate transformation rules, ...
     You then dispose directly of InverseMetric[i,j], CovariantChristoffel[i,j,k], 
       Christoffel[i,j,k], CovariantRiemann[i,j,k,l], Riemann[i,j,k,l], 
       Ricci[i,j], ScalarCurvature, Einstein[i,j,k,l], CovariantWeyl[i,j,k,l] and 
       Weyl[i,j,k,l]. 
       The same names without indices and preceded by Print display the 
       essential non-null components. 
       PrintRicciPacakage[] displays all these quantities. 
     The intermediate quantities will be computed as needed. 
     But if you prefer to compute systematically all these quantities: 
  6. Use ComputeRicciPackage[] or, in turn, ComputeInverseMetric[], 
     ComputeCovariantChristoffel[], ComputeChristoffel[], 
     ComputeRiemann[], ComputeRicci[], 
     ComputeScalarCurvature[], ComputeEinstein[] and ComputeWeyl[]. 
     These steps are optional: the first time you use a magnitude (to examine
     or print it, for instance) it will be automatically computed if necessary  
  7. Use SaveRicci["file"] to save the results for later use, if necessary. 
  8. Use ClearRicciPackage[] (or <<Ricci.m) to reset the package before a new 
     problem is analyzed. 

  See the example in RicciX.m, which can also be used as a template file for 
  other problems.            

  This package and the definitions and conventions used are described in the 
  paper entitled "Computing the Ricci and Einstein tensors", by the same 
  author.

  To compute the same quantities in an orthonormal frame use the companion 
  Tetrad package. Since both package have some common names (specially 
  Coordinate), one should use <<Ricci.m to reset the package if Ricci and 
  Tetrad are being used simultaneously.
*)
    
BeginPackage["Ricci`Ricci`"]

SetMetric::usage = "\n
  SetMetric[ds2,vars,first,simplify] define the coordinates and\n
    metric tensor from a line element, ds2, like\n    
    -(Dt[x]-a^2 Dt[y])^2+Dt[z]^2+2 b Cos[t] Dt[t] Dt[z]\n
    (Notice the use of the \"differential\" Dt.)\n
    \"vars\" is an optional list of variables like {t,x,y,z}.\n
            (Useful to choose the orrder of coordinates.)\n
    \"first\" is the label of the first coordinate\n
            (it defaults to 0 in 4 dimensios and to 1 otherwise).\n
    \"simplify\" is an optional function to be applied to each\n
            metric component after computing it."

FirstIndex::usage = "\n
  FirstIndex is the label of the first coordinate.\n
  Usually 0 (the default value) or 1."

LastIndex::usage = "\n
  LastIndex is the label of the last coordinate.\n
  LastIndex = FirstIndex+Dimension-1. By default 3."

Coordinate::usage = "\n
  Coordinate[i] is the name for the i-th coordinate."

MetricTensor::usage = "\n
  MetricTensor[i,j] are the components of the covariant metric."

DetMetric::usage = "\n
  DetMetric is teh determinant of the covariant metric."

InverseMetric::usage = "\n
  InverseMetric[i,j] are the components of the contravariant metric."

CovariantChristoffel::usage = "\n
  CovariantChristoffel[i,j,k] are the connection symbols\n
  with all their indices lowered."

Christoffel::usage = "\n
  Christoffel[i,j,k] are the connection symbols.\n
  The first index is the contravariant one."

CovariantRiemann::usage = "\n
  CovariantRiemann[i,j,k,l] are the covariant components\n
  of the Riemann tensor."

Riemann::usage = "\n
  Riemann[i,j,k,l] are the components of the Riemann tensor.\n
  The first index is the contravariant one."

Ricci::usage = "\n
  Ricci[i,j] are the components of the Ricci tensor."

ScalarCurvature::usage = "\n
  ScalarCurvature is the contraction of the Ricci tensor with the metric."

Einstein::usage = "\n
  Einstein[i,j] are the components of the Einstein tensor."

CovariantWeyl::usage = "\n
  CovariantWeyl[i,j,k,l] are the covariant components of the Weyl tensor."

Weyl::usage = "\n
  Weyl[i,j,k,l] are the components of the Weyl tensor.\n
  The first index is the contravariant one."

RicciSimplify::usage = "\n
  RicciSimplify is the simplification function\n
  to be applied at each step in the Ricci package."
 
ComputeInverseMetric::usage = "\n
  ComputeInverseMetric[] generates the contravariant metric tensor."

ComputeCovariantChristoffel::usage = "\n
  ComputeCovariantChristoffel[] generates the\n
  covariant connection symbols."

ComputeChristoffel::usage = "\n
  ComputeChristoffel[] generates the connection symbols."

ComputeRiemann::usage = "\n
  ComputeRiemann[] generates the (covariant) Riemann tensor."

ComputeRicci::usage = "\n
  ComputeRicci[] generates the Ricci tensor."

ComputeScalarCurvature::usage = "\n
  ComputeScalarCurvature[] generates the scalar curvature."

ComputeEinstein::usage = "\n
  ComputeEinstein[] generates the Einstein tensor."

ComputeWeyl::usage = "\n
  ComputeWeyl[] generates the (covariant) Weyl tensor."

ComputeRicciPackage::usage = "\n
  ComputeRicciPackage[] generates all the quantities in the Ricci package."

SaveRicci::usage = "\n
  SaveRicci[file] saves in file all the quantities in the Ricci package."

PrintCoordinate::usage = "\n
PrintCoordinate[] displays the actual names of Coordinate[i].\n
PrintCoordinate[Short] displays the short form."

PrintMetricTensor::usage = "\n
PrintMetricTensor[] displays the non-null components of the MetricTensor.\n
PrintMetricTensor[Short] displays the short form."

PrintInverseMetric::usage = "\n
PrintInverseMetric[] displays the non-null components of the InverseMetric.\n
PrintInverseMetric[Short] displays the short form."

PrintCovariantChristoffel::usage = "\n
PrintCovariantChristoffel[] displays the non-null components\n
  of CovariantChristoffel.\n
PrintCovariantChristoffel[Short] displays the short form."

PrintChristoffel::usage = "\n
PrintChristoffel[] displays the non-null components of Christoffel.\n
PrintChristoffel[Short] displays the short form."

PrintCovariantRiemann::usage = "\n
PrintCovariantRiemann[] displays the non-null components of CovariantRiemann.\n
PrintCovariantRiemann[Short] displays the short form."

PrintRiemann::usage = "\n
PrintRiemann[] displays the non-null components of Riemann.\n
PrintRiemann[Short] displays the short form."

PrintRicci::usage = "\n
PrintRicci[] displays the non-null components of Ricci.\n
PrintRicci[Short] displays the short form."

PrintScalarCurvature::usage = "\n 
PrintScalarCurvature[] displays the scalar curvature.\n
PrintScalarCurvature[Short] displays the short form."

PrintEinstein::usage = "\n
PrintEinstein[] displays the non-null components of Einstein.\n
PrintEinstein[Short] displays the short form."

PrintCovariantWeyl::usage = "\n
PrintCovariantWeyl[] displays the non-null components of CovariantWeyl.\n
PrintCovariantWeyl[Short] displays the short form."

PrintWeyl::usage = "\n
PrintWeyl[] displays the non-null components of Weyl.\n
PrintWeyl[Short] displays the short form."

PrintRicciPackage::usage = "\n
  PrintRicciPackage[] displays all the quantities computed in the package.\n
PrintRicciPackage[Short] displays the short forms."

ClearRicciPackage::usage = "\n
  ClearRicciPackage[] resets the Ricci package\n
  before a new metric is analyzed."


Begin["`Private`"]

ClearRicciPackage[] := (

FirstIndex = 0;

LastIndex = 3;

RicciSimplify = Together;

Clear[Coordinate];

Clear[MetricTensor];

(*
Module[
  {i,j},
  Do [MetricTensor[i,j] = 0,
      {i,FirstIndex,LastIndex},
      {j,i         ,LastIndex}];];
*)

MetricTensor[i_,j_] := 0;

SetAttributes[MetricTensor,Orderless];

(* :Note:
  Enter metric from line element 
*)

SetMetric::variables = "Variables `1` and differentials `2` do not match.";

SetMetric[ds2_,var_List,frst_Integer,simp_:Identity] := 
  Module[{i,j,l=Length[var]},
    If[Sort[var]=!=Sort[GetVars[ds2]],
      Message[SetMetric::variables,var,Dt /@ GetVars[ds2]]];
    FirstIndex = frst;
    LastIndex  = frst+l-1;
    Do[Coordinate[i+frst-1]            = var[[i]],{i,l}];
    Do[MetricTensor[i+frst-1,j+frst-1] =
         If[i==j, Coefficient[ds2,Dt[var[[i]]],2],
                  Coefficient[Coefficient[ds2,Dt[var[[i]]]],Dt[var[[j]]]]/2
           ] // simp;
      ,{i,l},{j,i,l}];
];

SetMetric[ds2_,var_List,simp_:Identity] :=
	 SetMetric[ds2,var,If[Length[var] == 4,0,1],simp];

SetMetric[ds2_,frst_Integer,simp_:Identity] :=
	 SetMetric[ds2,GetVars[ds2],frst,simp];

SetMetric[ds2_,simp_:Identity] :=
	 SetMetric[ds2,GetVars[ds2],simp];

GetVars[ds2_] := (Select[Variables[ds2],(Head[#] == Dt)&]/. Dt -> Identity);

(* :Note:
  DetMetric is computed only once.
*)

Clear[DetMetric];

DetMetric := DetMetric = Module[
  {i,j,m},
  m = Det[Table[MetricTensor[i,j],
                          {i,FirstIndex,LastIndex},
                          {j,FirstIndex,LastIndex}]]
                          // RicciSimplify];

Clear[InverseMetric];

SetAttributes[InverseMetric,Orderless];

(* :Note:
  The InverseMetric is computed only once.
*)

InverseMetric::singular = "The metric is singular.";

InverseMetric[k_,l_] := Module[
  {i,j,m},
  m = Check[Inverse[Table[MetricTensor[i,j],
                          {i,FirstIndex,LastIndex},
                          {j,FirstIndex,LastIndex}]],
            "RicciMeaninless"];
  If[m == "RicciMeaninless",
     Message[InverseMetric::singular]; 
     m = IdentityMatrix[LastIndex-FirstIndex+1]
    ];
  Do [InverseMetric[i,j] = m[[i-FirstIndex+1,j-FirstIndex+1]] // RicciSimplify,
      {i,FirstIndex,LastIndex},
      {j,i         ,LastIndex}];
  InverseMetric[k,l]];

Clear[CovariantChristoffel];
Clear[Christoffel];

CovariantChristoffel[i_,j_,k_] := CovariantChristoffel[i,k,j]   /; j > k;
Christoffel[i_,j_,k_]          := Christoffel[i,k,j]            /; j > k;

(* :Note:
  Several quantities are computed only once
*)

CovariantChristoffel[p_,q_,r_] := Module[
  {i,j,k,d},
  Do [d[i,j,k] = d[j,i,k] = 
        D[MetricTensor[i,j],Coordinate[k]] // RicciSimplify,
      {i,FirstIndex,LastIndex},
      {j,i         ,LastIndex},
      {k,FirstIndex,LastIndex}];
  Do [CovariantChristoffel[i,j,k] = (d[i,j,k]+d[i,k,j]-d[j,k,i])/2 // 
          RicciSimplify,
      {i,FirstIndex,LastIndex},
      {j,FirstIndex,LastIndex},
      {k,j         ,LastIndex}];  
  CovariantChristoffel[p,q,r]];

Christoffel[p_,q_,r_] := Module[
  {i,j,k,l},
  Do [Christoffel[i,j,k] =
        Sum[InverseMetric[i,l] CovariantChristoffel[l,j,k],
            {l,FirstIndex,LastIndex}] // RicciSimplify,
      {i,FirstIndex,LastIndex},
      {j,FirstIndex,LastIndex},
      {k,j         ,LastIndex}];
  Christoffel[p,q,r]];

Clear[CovariantRiemann];

CovariantRiemann[i_,i_,k_,l_] :=  0;        
CovariantRiemann[i_,j_,k_,k_] :=  0;        
CovariantRiemann[i_,j_,k_,l_] := -CovariantRiemann[j,i,k,l]     /; i > j;
CovariantRiemann[i_,j_,k_,l_] := -CovariantRiemann[i,j,l,k]     /; k > l;
CovariantRiemann[i_,j_,k_,l_] :=  CovariantRiemann[k,l,i,j]     /;
                          i < j && k < l && (i > k || (i == k && j > l));
CovariantRiemann[i_,j_,k_,l_] :=  CovariantRiemann[i,l,k,j] -
                                  CovariantRiemann[i,k,l,j] // RicciSimplify /;
                                                         (i < k < l < j);

CovariantRiemann[i_,j_,k_,l_] := CovariantRiemann[i,j,k,l] = Module[
  {m},
  D[CovariantChristoffel[i,j,l],Coordinate[k]] -
  D[CovariantChristoffel[i,j,k],Coordinate[l]] +
  Sum[Christoffel[m,i,l] CovariantChristoffel[m,j,k] - 
      Christoffel[m,i,k] CovariantChristoffel[m,j,l] ,
      {m,FirstIndex,LastIndex}] // RicciSimplify];

Clear[Riemann];

Riemann[i_,j_,k_,l_] := Module[
  {m},
  Sum[InverseMetric[i,m] CovariantRiemann[m,j,k,l],
      {m,FirstIndex,LastIndex}] // RicciSimplify];

Clear[Ricci];

SetAttributes[Ricci,Orderless];

Ricci[i_,j_] := Ricci[i,j] = Module[
  {k},
  Sum[Riemann[k,i,k,j],
      {k,FirstIndex,LastIndex}] // RicciSimplify];

Clear[ScalarCurvature];

ScalarCurvature := ScalarCurvature = Module[
  {i,j},
  Sum[InverseMetric[i,j] Ricci[i,j], 
      {i,FirstIndex,LastIndex},
      {j,FirstIndex,LastIndex}] // RicciSimplify];

Clear[Einstein];

SetAttributes[Einstein,Orderless];

Einstein[i_,j_] := Einstein[i,j] =
  Ricci[i,j]-ScalarCurvature/2 MetricTensor[i,j] // RicciSimplify;

Clear[CovariantWeyl];

CovariantWeyl[i_,i_,k_,l_] :=  0;        
CovariantWeyl[i_,j_,k_,k_] :=  0;        
CovariantWeyl[i_,j_,k_,l_] := -CovariantWeyl[j,i,k,l]     /; i > j;
CovariantWeyl[i_,j_,k_,l_] := -CovariantWeyl[i,j,l,k]     /; k > l;
CovariantWeyl[i_,j_,k_,l_] :=  CovariantWeyl[k,l,i,j]     /;
                          i < j && k < l && (i > k || (i == k && j > l));
CovariantWeyl[i_,j_,k_,l_] :=  CovariantWeyl[i,l,k,j] -
                               CovariantWeyl[i,k,l,j]     // RicciSimplify  /;
                                                         (i < k < l < j);

CovariantWeyl[i_,j_,k_,l_] := CovariantWeyl[i,j,k,l] = Module[
  {n1},
  If[(n1 = LastIndex-FirstIndex) < 3, 0,
     CovariantRiemann[i,j,k,l]-              (* + ??? *)
     1/(n1-1) (MetricTensor[i,k] Ricci[j,l]-MetricTensor[i,l] Ricci[j,k]+
               MetricTensor[j,l] Ricci[i,k]-MetricTensor[j,k] Ricci[i,l])+
     ScalarCurvature/n1/(n1-1) (MetricTensor[i,k] MetricTensor[j,l]-
                                MetricTensor[i,l] MetricTensor[j,k]) //
     RicciSimplify]];

Clear[Weyl];

Weyl[i_,j_,k_,l_] := Module[
  {m},
  Sum[InverseMetric[i,m] CovariantWeyl[m,j,k,l],
      {m,FirstIndex,LastIndex}] // RicciSimplify];

(* :Note:
  Compute...
*)

ComputeInverseMetric[] :=
  (InverseMetric[FirstIndex,FirstIndex];);

ComputeCovariantChristoffel[] :=
  (CovariantChristoffel[FirstIndex,FirstIndex,FirstIndex];);

ComputeChristoffel[] :=
  (Christoffel[FirstIndex,FirstIndex,FirstIndex];);

ComputeRiemann[] := Module[
  {i,j,k,l},
  Do [If [j <= l, CovariantRiemann[i,j,k,l]],
      {i,FirstIndex,LastIndex-1},
      {j,i+1       ,LastIndex  },
      {k,i         ,LastIndex-1},
      {l,k+1       ,LastIndex  }];];

ComputeRicci[] := Module[
  {i,j},
  Do [Ricci[i,j],
      {i,FirstIndex,LastIndex},
      {j,i         ,LastIndex}];];

ComputeScalarCurvature[] := (ScalarCurvature;);

ComputeEinstein[] := Module[
  {i,j},
  Do [Einstein[i,j],
      {i,FirstIndex,LastIndex},
      {j,i         ,LastIndex}];];

ComputeWeyl[] := Module[
  {i,j,k,l},
  Do [If [j <= l, CovariantWeyl[i,j,k,l]],
      {i,FirstIndex,LastIndex-1},
      {j,i+1       ,LastIndex  },
      {k,i         ,LastIndex-1},
      {l,k+1       ,LastIndex  }];];

ComputeRicciPackage[] := (
  Print["Computing the inverse metric"];
  ComputeInverseMetric[];
  Print["Computing the connection"];
  ComputeCovariantChristoffel[];
  ComputeChristoffel[];
  Print["Computing the Riemann tensor"];
  ComputeRiemann[];
  Print["Computing the Ricci tensor"];
  ComputeRicci[];
  Print["Computing the scalar curvature"];
  ComputeScalarCurvature[];
  Print["Computing the Einstein tensor"];
  ComputeEinstein[];
  Print["Computing the Weyl tensor"];
  ComputeWeyl[];
);

(* :Note:
  Save results
*)

SaveRicci[file_] :=
  Save[file,
       FirstIndex,
       LastIndex,
       Coordinate,
       MetricTensor,
       InverseMetric,
       CovariantChristoffel,
       Christoffel,
       CovariantRiemann,
       Riemann,
       Ricci,
       ScalarCurvature,
       Einstein,
       CovariantWeyl,
       Weyl];

(* :Note:
  Display results
*)

PrintCoordinate[format_:Identity] := Module[
  {i}, 
  Print[];
  SequenceForm @@
    Drop[
      Flatten[Table[{"x",Superscript[i]," = ",Coordinate[i] // format,",  "},
                    {i,FirstIndex,LastIndex}]],
      -1] // Print;];

PrintMetricTensor[format_:Identity] := Module[
  {i,j,z}, 
  z = False;
  Do[If[MetricTensor[i,j] =!= 0,
        z = True;
        Print[];
        Print["g",Subscript[i],Subscript[j]," = ",MetricTensor[i,j] // format]],
     {i,FirstIndex,LastIndex},{j,i,LastIndex}];
  If[z == False, Print[]; Print["g",Subscript["ij"]," = 0"]];];

PrintInverseMetric[format_:Identity] := Module[
  {i,j,z}, 
  z = False;
  Do[If[InverseMetric[i,j] =!= 0,
        z = True;
        Print[]; 
        Print["g",Superscript[i],Superscript[j]," = ",InverseMetric[i,j] // format]],
     {i,FirstIndex,LastIndex},{j,i,LastIndex}];
  If[z == False, Print[]; Print["g",Superscript["ij"]," = 0"]];];

PrintCovariantChristoffel[format_:Identity] := Module[
  {i,j,k,z}, 
  z = False;
  Do[If[CovariantChristoffel[i,j,k] =!= 0,
        z = True;
        Print[]; 
        Print["\[CapitalGamma]",Subscript[i],Subscript[j],Subscript[k]," = ",
              CovariantChristoffel[i,j,k] // format]],
     {i,FirstIndex,LastIndex},{j,FirstIndex,LastIndex},{k,j,LastIndex}];
  If[z == False, Print[]; Print["\[CapitalGamma]",Subscript["ijk"]," = 0"]];];

PrintChristoffel[format_:Identity] := Module[
  {i,j,k,z}, 
  z = False;
  Do[If[Christoffel[i,j,k] =!= 0,
        z = True;
        Print[]; 
        Print["\[CapitalGamma]",Superscript[i],Subscript[j],Subscript[k]," = ",
              Christoffel[i,j,k] // format]],
     {i,FirstIndex,LastIndex},{j,FirstIndex,LastIndex},{k,j,LastIndex}];
  If[z == False, Print[]; Print["\[CapitalGamma]",Superscript["i"],Subscript["jk"]," = 0"]];];

PrintCovariantRiemann[format_:Identity] := Module[
  {i,j,k,l,z}, 
  z = False;
  Do[If[CovariantRiemann[i,j,k,l] =!= 0,
        z = True;
        Print[]; 
        Print["R",Subscript[i],Subscript[j],Subscript[k],Subscript[l]," = ",
              CovariantRiemann[i,j,k,l] // format]],
     {i,FirstIndex,LastIndex-1},{j,    i+1   ,LastIndex},
     {k,i         ,LastIndex-1},{l,Max[k+1,j],LastIndex  }];
  If[z == False, Print[]; Print["R",Subscript["ijkl"]," = 0"]];];

PrintRiemann[format_:Identity] := Module[
  {i,j,k,l,z}, 
  z = False;
  Do[If[Riemann[i,j,k,l] =!= 0,
        z = True;
        Print[]; 
        Print["R",Superscript[i],Subscript[j],Subscript[k],Subscript[l]," = ",
              Riemann[i,j,k,l] // format]],
     {i,FirstIndex,LastIndex-1},{j,    i+1   ,LastIndex},
     {k,i         ,LastIndex-1},{l,Max[k+1,j],LastIndex  }];
  If[z == False, Print[]; Print["R",Superscript["i"],Subscript["jkl"]," = 0"]];];

PrintRicci[format_:Identity] := Module[
  {i,j,z}, 
  z = False;
  Do[If[Ricci[i,j] =!= 0,
        z = True;
        Print[]; 
        Print["R",Subscript[i],Subscript[j]," = ",Ricci[i,j] // format]],
     {i,FirstIndex,LastIndex},{j,i,LastIndex}];
  If[z == False, Print[]; Print["R",Subscript["ij"]," = 0"]];];

PrintScalarCurvature[format_:Identity] := (Print[]; Print["R = ",ScalarCurvature // format];);

PrintEinstein[format_:Identity] := Module[
  {i,j,z}, 
  z = False;
  Do[If[Einstein[i,j] =!= 0,
        z = True;
        Print[]; 
        Print["G",Subscript[i],Subscript[j]," = ",Einstein[i,j] // format]],
     {i,FirstIndex,LastIndex},{j,i,LastIndex}];
  If[z == False, Print[]; Print["G",Subscript["ij"]," = 0"]];];

PrintCovariantWeyl[format_:Identity] := Module[
  {i,j,k,l,z}, 
  z = False;
  Do[If[CovariantWeyl[i,j,k,l] =!= 0,
        z = True;
        Print[]; 
        Print["C",Subscript[i],Subscript[j],Subscript[k],Subscript[l]," = ",
              CovariantWeyl[i,j,k,l] // format]],
     {i,FirstIndex,LastIndex-1},{j,    i+1   ,LastIndex},
     {k,i         ,LastIndex-1},{l,Max[k+1,j],LastIndex  }];
  If[z == False, Print[]; Print["C",Subscript["ijkl"]," = 0"]];];

PrintWeyl[format_:Identity] := Module[
  {i,j,k,l,z}, 
  z = False;
  Do[If[Weyl[i,j,k,l] =!= 0,
        z = True;
        Print[]; 
        Print["C",Superscript[i],Subscript[j],Subscript[k],Subscript[l]," = ",
              Weyl[i,j,k,l] // format]],
     {i,FirstIndex,LastIndex-1},{j,    i+1   ,LastIndex},
     {k,i         ,LastIndex-1},{l,Max[k+1,j],LastIndex  }];
  If[z == False, Print[]; Print["C",Superscript["i"],Subscript["jkl"]," = 0"]];];

PrintRicciPackage[format_:Identity] := (
       PrintCoordinate[format];
       PrintMetricTensor[format];
(*     PrintInverseMetric[format];           *)
       PrintCovariantChristoffel[format];
(*     PrintChristoffel[format];             *)
       PrintCovariantRiemann[format];
(*     PrintRiemann[format];                 *)
       PrintRicci[format];
       PrintScalarCurvature[format];
       PrintEinstein[format];
       PrintCovariantWeyl[format];
(*     PrintWeyl[format];                    *)
     );

);                      

(* :Note:
  End of ClearRicciPackage and do it!
*)

ClearRicciPackage[];

End[]

EndPackage[]     

(* ::Package:: *)

(* ::Subsection:: *)
(*How to solve the unstable boundary condition system*)
(**)


(* ::Input:: *)
(*ClearAll[x,B,a,c,b,f,M]*)
(*(*The smallest example of a similar unstable system is*)*)
(**)
(*M  =\!\(\**)
(*TagBox[*)
(*RowBox[{"(", "", GridBox[{*)
(*{*)
(*RowBox[{*)
(*SubscriptBox["J", "n"], "[", *)
(*RowBox[{"a1", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}], *)
(*RowBox[{*)
(*SubscriptBox["H", "n"], "[", *)
(*RowBox[{"a1", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}]},*)
(*{*)
(*RowBox[{"-", *)
(*RowBox[{*)
(*SubscriptBox["J", "n"], "[", *)
(*RowBox[{"a2", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}]}], *)
(*RowBox[{"-", *)
(*RowBox[{*)
(*SubscriptBox["H", "n"], "[", *)
(*RowBox[{"a2", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}]}]}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], "", ")"}],*)
(*Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);*)
(*B= {0,fo[n] Subscript[J, n][ko a[1]]};*)
(*B= {-fo[n],fo[n] };*)
(**)
(*subN = {a2->2.`,a1->1.5`,ko->0.0001`,k[1]->2.`,Subscript[J, n][ko a[1]]->1.`,fo[n]->1,Subscript[H, n_][x_]->HankelH1[n,x],Subscript[J, n_][x_]->BesselJ[n,x],Derivative[1][Subscript[H, n_]][x_]->1/2 (HankelH1[-1+n,x]-HankelH1[1+n,x]),Derivative[1][Subscript[J, n_]][x_]->1/2 (BesselJ[-1+n,x]-BesselJ[1+n,x])};*)
(*(*subN = subN/.BesselJ ->HankelH2;*)*)
(**)


(* ::Input:: *)
(*dn = 3;*)
(*ns = Range[0,35,dn];*)
(**)
(*NM = M//.subN;*)
(*NB = B//.subN;*)
(*(*NB = NM.xo*)*)


(* ::Input:: *)
(*{Abs@Det[NM/.n->#],SingularValueList[NM/.n->#]}&/@ns*)


(* ::Input:: *)
(**)
(*(*which has the well posed solution*)*)
(*eqs = M . {a[n],c[n]} -NB;*)
(*\[Epsilon] = 10^-14;*)
(*subsol=Solve[eqs==RandomReal[\[Epsilon],2] +RandomReal[\[Epsilon],2] I,{a[n],c[n]}]//Simplify//Flatten;*)
(*subsol2=Solve[eqs==RandomReal[\[Epsilon],2] +RandomReal[\[Epsilon],2] I,{a[n],c[n]}]//Simplify//Flatten;*)
(**)
(*(*test the solution*)*)
(*subNsols= (subsol//.subN//Flatten)/.n->#&/@ns ;*)
(*subNsols2= (subsol2//.subN//Flatten)/.n->#&/@ns ;*)
(**)
(*sols = Transpose[#/.Rule-> List][[2]]&/@subNsols;*)
(*sols2 = Transpose[#/.Rule-> List][[2]]&/@subNsols2;*)
(**)
(*diffsols = sols-sols2;*)
(*Norm/@diffsols*)
(**)
(**)
(*#[[1]]&/@diffsols;*)
(*Abs[% *( BesselJ[#,a1 k[1]//.subN]&/@ns)]*)
(**)
(**)
(*#[[2]]&/@diffsols;*)
(*Abs[% *( HankelH1[#,a1 k[1]//.subN]&/@ns)]*)


(* ::Input:: *)
(*invM = Inverse[M];*)
(*Nxs = Table[*)
(**)
(*invNM =invM//.subN/.n->n1;*)
(*sol =invNM . (NB +RandomReal[\[Epsilon],2] +RandomReal[\[Epsilon],2] I//.subN/.n->n1) *)
(*,{n1,ns[[1]],ns//Last,dn}];*)
(**)


(* ::Input:: *)
(*sublinsols =Flatten[ Thread[({a[n],c[n]}/.subN/.n->#)->LinearSolve[NM/.n->#,NB+ RandomReal[{-10^-15,10^-15},2]/.n->#]]&/@ns]//Quiet;*)
(**)
(**)


(* ::Subsubsection:: *)
(*Try to condition Matrix *)


(* ::Input:: *)
(**)
(*(*The smallest example of a similar unstable system is*)*)
(**)
(*M  =\!\(\**)
(*TagBox[*)
(*RowBox[{"(", "", GridBox[{*)
(*{"1", "1"},*)
(*{*)
(*RowBox[{*)
(*RowBox[{"-", *)
(*RowBox[{*)
(*SubscriptBox["J", "n"], "[", *)
(*RowBox[{"a2", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}]}], "/", *)
(*RowBox[{"(", *)
(*RowBox[{*)
(*RowBox[{*)
(*SubscriptBox["J", "n"], "[", *)
(*RowBox[{"a2", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}], "+", *)
(*RowBox[{*)
(*SubscriptBox["J", "n"], "[", *)
(*RowBox[{"a1", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}]}], ")"}]}], *)
(*RowBox[{*)
(*RowBox[{"-", *)
(*RowBox[{*)
(*SubscriptBox["H", "n"], "[", *)
(*RowBox[{"a2", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}]}], "/", *)
(*RowBox[{"(", *)
(*RowBox[{*)
(*RowBox[{*)
(*SubscriptBox["H", "n"], "[", *)
(*RowBox[{"a2", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}], "+", *)
(*RowBox[{*)
(*SubscriptBox["H", "n"], "[", *)
(*RowBox[{"a1", " ", *)
(*RowBox[{"k", "[", "1", "]"}]}], "]"}]}], ")"}]}]}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], "", ")"}],*)
(*Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);*)
(*B= {0,fo[n] Subscript[J, n][ko a[1]]};*)
(*B= {-fo[n],fo[n] };*)
(**)
(*subN = {a2->2.`,a1->1.5`,ko->0.0001`,k[1]->2.`,Subscript[J, n][ko a[1]]->1.`,fo[n]->1,Subscript[H, n_][x_]->HankelH1[n,x],Subscript[J, n_][x_]->BesselJ[n,x],Derivative[1][Subscript[H, n_]][x_]->1/2 (HankelH1[-1+n,x]-HankelH1[1+n,x]),Derivative[1][Subscript[J, n_]][x_]->1/2 (BesselJ[-1+n,x]-BesselJ[1+n,x])};*)
(*(*subN = subN/.BesselJ ->HankelH2;*)*)
(**)


(* ::Input:: *)
(*dn = 3;*)
(*ns = Range[0,35,dn];*)
(**)
(*NM = M//.subN;*)
(*NB = B//.subN;*)
(*(*NB = NM.xo*)*)


(* ::Input:: *)
(*{Abs@Det[NM/.n->#],SingularValueList[NM/.n->#]}&/@ns*)


(* ::Input:: *)
(**)
(*(*which has the well posed solution*)*)
(*eqs = M . {a[n],c[n]} -NB;*)
(*\[Epsilon] = 10^-14;*)
(*subsol=Solve[eqs==RandomReal[\[Epsilon],2] +RandomReal[\[Epsilon],2] I,{a[n],c[n]}]//Simplify//Flatten;*)
(*subsol2=Solve[eqs==RandomReal[\[Epsilon],2] +RandomReal[\[Epsilon],2] I,{a[n],c[n]}]//Simplify//Flatten;*)
(**)
(*(*test the solution*)*)
(*subNsols= (subsol//.subN//Flatten)/.n->#&/@ns ;*)
(*subNsols2= (subsol2//.subN//Flatten)/.n->#&/@ns ;*)
(**)
(*sols = Transpose[#/.Rule-> List][[2]]&/@subNsols;*)
(*sols2 = Transpose[#/.Rule-> List][[2]]&/@subNsols2;*)
(**)
(*diffsols = sols-sols2;*)
(*Norm/@diffsols*)
(**)
(*#[[1]]&/@diffsols*)
(**)
(*#[[2]]&/@diffsols*)
(**)


(* ::Subsection:: *)
(*4x4 matrix*)


(* ::Input:: *)
(*ClearAll[b,a,c,b,d];*)
(*(*The smallest example of a similar unstable system is*)*)
(*M = {*)
(*	{Subscript[J, n][a1 k[1]], Subscript[H, n][a1 k[1]],Subscript[J, n][a1 k[2]],Subscript[H, n][a1 k[2]]},*)
(*{Derivative[1][Subscript[J, n]][a1 k[1]], Derivative[1][Subscript[H, n]][a1 k[1]],Derivative[1][Subscript[J, n]][a1 k[2]],Derivative[1][Subscript[H, n]][a1 k[2]]},*)
(*	{Subscript[J, n][a2 k[1]], Subscript[H, n][a2 k[1]],Subscript[J, n][a2 k[2]],Subscript[H, n][a2 k[2]]},*)
(*         {Derivative[1][Subscript[J, n]][a2 k[1]], Derivative[1][Subscript[H, n]][a2 k[1]],Derivative[1][Subscript[J, n]][a2 k[2]],Derivative[1][Subscript[H, n]][a2 k[2]]}*)
(*};*)
(**)
(*no = 30;*)
(**)
(*x =RandomReal[1,4] + RandomReal[1,4] I;*)
(*B = M . x //.subN/.n->no*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)
(*B = {0,fo[n] Subscript[J, n][ko a1],0, fo[n] Derivative[1][Subscript[J, n]][ko a1]};*)
(*B = {0,fo[n],0, fo[n]};*)
(**)
(*subN = {a2->2.0,a1->1.5,ko-> 0.1,k[1] -> 2.0,k[2] -> 3.2,Subscript[J, n][ko a[1]]->1.0,fo[n]->1, Subscript[H, n_][x_] ->HankelH1[n,x], Subscript[J, n_][x_] ->BesselJ[n,x], Derivative[1][Subscript[H, n_]][x_] ->D[HankelH1[n,x],x],Derivative[1][Subscript[J, n_]][x_] ->D[BesselJ[n,x],x]};*)
(**)
(*NM = M//.subN;*)
(*NB = B//.subN;*)
(**)
(*{Abs@Det[NM/.n->#],SingularValueList[NM/.n->#]}&/@ns //Quiet*)


(* ::Input:: *)
(**)
(*(*which has the well posed solution*)*)
(*vars = {a[n],b[n],c[n],d[n]};*)
(*eqs = M . vars - B;*)
(**)
(*subsol=Solve[eqs==0,Evaluate@vars]//Simplify//Flatten;*)
(**)
(*(*test the solution*)*)
(*subNsols= (subsol//.subN//Flatten)/.n->#&/@ns //Flatten;*)
(**)
(*Neqs = eqs //.subN/.n->#&/@ns ;*)
(*errors=Norm/@(%/.subNsols)*)
(**)


(* ::Input:: *)
(*sublinsols =Flatten[ Thread[(vars/.subN/.n->#)->LinearSolve[NM/.n->#,NB+ RandomReal[{-10^-15,10^-15},4]/.n->#]]&/@ns]//Quiet;*)
(*Neqs/.sublinsols;*)
(*errors=Norm/@(%/.subNsols)*)


(* ::Input:: *)
(**)
(**)
(**)
(*subNsol=subsol//.subN/.n->1;*)
(*eqs//.subN/.n->1/.subNsol*)


(* ::Input:: *)
(*ClearAll[\[Psi]o,\[Psi],J,H]*)
(*\[Psi]o = fo[n] Subscript[J, n][ko r];*)
(*\[Psi][s] = Ao[n]Subscript[H, n][ko r];*)
(*\[Psi][1] = f[1,n] Subscript[J, n][k[1] r] + A[1,n]Subscript[H, n][k[1] r];*)
(*eqs = {\[Psi][1]  /.r-> a[0], \[Psi][s] + \[Psi]o  - \[Psi][1]/.r-> a[1],D[\[Psi][1],r]/\[Rho][1]  - D[\[Psi][s] + \[Psi]o ,r]/\[Rho]o/.r-> a[1]};*)
(**)
(*subsol = Solve[Thread[eqs==0],{Ao[n],f[1,n],A[1,n]}]  //Simplify;*)


(* ::Input:: *)
(*subN = {a[1]->2.0,a[0]->1.5,\[Rho]o->0.00001,\[Rho][1]->2.0,ko-> 0.0001,k[1] -> 2.0,Subscript[J, n][ko a[1]]->1.0,fo[n]->1, Subscript[H, n_][x_] ->HankelH1[n,x], Subscript[J, n_][x_] ->BesselJ[n,x], Derivative[1][Subscript[H, n_]][x_] ->D[HankelH1[n,x],x],Derivative[1][Subscript[J, n_]][x_] ->D[BesselJ[n,x],x]};*)
(*subsol//.subN//Flatten;*)


(* ::Input:: *)
(*%/.n->#&/@Range[1,20,3]*)
(**)


(* ::Input:: *)
(*(*Capsule with given forcing on boundary*)*)
(*ClearAll[\[Psi]o,\[Psi],J,H]*)
(*\[Psi]o = fo[n] Subscript[J, n][ko r];*)
(*\[Psi][1] = f[1,n] Subscript[J, n][k[1] r] + A[1,n]Subscript[H, n][k[1] r];*)
(*eqs = {\[Psi][1]  /.r-> a[0],  \[Psi]o  - \[Psi][1]/.r-> a[1]};*)
(*subsol = Solve[Thread[eqs==0],{f[1,n],A[1,n]}] /.{\[Rho]o -> qo ko,\[Rho][1] -> q[1] k[1],\[Rho][0] -> q[0] k[0]} //Simplify*)
(**)
(**)


(* ::Input:: *)
(*vars = {f[1,n],A[1,n]};*)
(*M = Coefficient[eqs,#]&/@vars;*)
(*M = Transpose@M;*)
(*b =  eqs - M . vars //Simplify;*)
(*M//MatrixForm*)


(* ::Input:: *)
(*subN = {a[1]->2.0,a[0]->1.5,k[1] -> 2.0,Subscript[J, n][ko a[1]]->1.0,fo[n]->1, Subscript[H, n_][x_] ->HankelH1[n,x], Subscript[J, n_][x_] ->BesselJ[n,x]};*)
(*ns = Range[1,40,3];*)
(*subsol//.subN//Flatten;*)
(*subNsol=Flatten[%/.n->#&/@ns];*)
(**)
(*eqs//.subN/.n->#&/@ns;*)
(*Norm/@(%/.subNsol)*)


(* ::Input:: *)
(*{Abs@f[1,n]Abs@BesselJ[n,a[1] k[1]],Abs@A[1,n]Abs@HankelH1[n,a[1] k[1]]}//.subN//Flatten;*)
(*%/.n->#&/@ns;*)
(*%/.Flatten@subNsol*)
(**)


(* ::Input:: *)
(*(*Is the matrix system ill-posed?*)*)
(*NM = M //.subN;*)
(*Nb = b //.subN;*)
(*{Abs@Det[NM/.n->#],SingularValueList[NM/.n->#]}&/@ns*)
(**)


(* ::Input:: *)
(*(*If we numerically solve the system:*)*)
(*sols =Flatten[ Thread[(vars/.subN/.n->#)->LinearSolve[NM/.n->#,-Nb+ RandomReal[{-10^-14,10^-14},2]/.n->#]]&/@ns]//Quiet;*)
(*eqs//.subN/.n->#&/@ns;*)
(*Norm/@(%/.sols)*)
(**)


(* ::Input:: *)
(*invM = Inverse[M]*)
(*invM . b*)
(*NinvM =invM //.subN;*)
(*Nb = b //.subN;*)


(* ::Input:: *)
(*invsols =Flatten[ Thread[(vars/.subN/.n->#)->(-NinvM . (Nb + RandomReal[{10^-18,10^-15}])/.n->#)]&/@ns];*)
(*Norm/@(eqs//.subN/.n->#&/@ns/.invsols)*)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*eqs//.subN/.n->#&/@ns/.invsols*)
(**)


(* ::Input:: *)
(*Nx = -NinvM . Nb/.n-> 1*)
(*M . Nx//.subN/.n-> 1*)
(*b//.subN/.n-> 1*)
(**)
(**)


(* ::Input:: *)
(*#->#&/@ns*)


(* ::Input:: *)
(*(* For one particle *)*)
(*(*Clear[dq0,dq]*)
(*denom = Y[a[1] k[1],a[0] k[1]] Derivative[1][Subscript[J, n]][a[0] k[0]]-dq0 Subscript[J, n][a[0] k[0]] Yd[a[1] k[1],a[0] k[1]];*)
(*force = (Ao[n] Subscript[H, n][ko a[1]]+fo[n] Subscript[J, n][ko a[1]])/denom;*)
(**)
(*numer = - dq0  Yd[a[0] k[1],a[0] k[1]]*)
(*force numer  == f[0,n]/.subsol [[1]]/.{dq0 \[Rule] q[0]/q[1]}//Simplify*)
(*numer = dq0 Subscript[J, n][a[0] k[0]] Derivative[1][Subscript[H, n]][a[0] k[1]]- Subscript[H, n][a[0] k[1]] Derivative[1][Subscript[J, n]][a[0] k[0]];*)
(*force numer  ==f[1,n]/.subsol [[1]]/.{dq0 \[Rule] q[0]/q[1]}//Simplify*)
(**)
(*numer =  Subscript[J, n][a[0] k[1]] Derivative[1][Subscript[J, n]][a[0] k[0]]-dq0 Subscript[J, n][a[0] k[0]] Derivative[1][Subscript[J, n]][a[0] k[1]]*)
(*force  numer== A[1,n]/.subsol [[1]]/.{dq0 \[Rule] q[0]/q[1]}//Simplify*)
(**)*)


(* ::Input:: *)
(*D[Y[x,y],x,y]*)
(**)
(*Yd [y,x]+Yd [x,y]//Simplify*)
(**)
(*capsule = tmp2;*)


Test subtracting the frequency response of the capsule . 
I'll assume the background is water, that there is Sunflower oil
in the tube, and the tube itself is silica;




subN = {Subscript[a, 1] -> 1, Subscript[a, 0] -> 3/4., qo-> 1, Subscript[q, 1]-> 6, Subscript[q, 0]-> 0.9, Subscript[k, 1]-> ko/3., Subscript[k, 0] -> 1.1 ko, fo[_] -> 1};
subHJ = { 
\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{
RowBox[{"Derivative", "[", "1", "]"}], "[", 
RowBox[{"Subscript", "[", 
RowBox[{"H", ",", "n_"}], "]"}], "]"}], "[", "k_", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) -> D[HankelH1[n,\[FormalK]],\[FormalK]]/.\[FormalK]->k, 
\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{
RowBox[{"Derivative", "[", "1", "]"}], "[", 
RowBox[{"Subscript", "[", 
RowBox[{"J", ",", "n_"}], "]"}], "]"}], "[", "k_", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) ->  D[BesselJ[n,\[FormalK]],\[FormalK]]/.\[FormalK]->k, 
\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{"Subscript", "[", 
RowBox[{"H", ",", "n_"}], "]"}], "[", "k_", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) :> HankelH1[n,k], 
\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{"Subscript", "[", 
RowBox[{"J", ",", "n_"}], "]"}], "[", "k_", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) :> BesselJ[n,k]
};

tube = capsule/.{Subscript[k, 0]-> ko, Subscript[q, 0]-> qo} //FullSimplify;
core = capsule/.{Subscript[k, 1]-> ko, Subscript[q, 1]-> qo} //FullSimplify;

tubeSubtract = capsule - tube  /.subHJ/.subN //Simplify;
Ncore = core/.subHJ/.subN //Simplify;


subN = {Subscript[a, 1] -> 1, Subscript[a, 0] -> 2., qo-> 4., 
    Subscript[q, 1]-> 1, Subscript[q, 0]-> 9., ko->0.5,
   Subscript[k, 1]-> 1., Subscript[k, 0] -> 1/3., fo[_] -> 1
};

Activate@capsule/.subHJ/.subN/.{n->1}




Activate@Inactive[Y][1\!\(\*
TagBox["*",
"InactiveToken",
BaseStyle->"Inactive",
SyntaxForm->"*"]\)1.`,2.`\!\(\*
TagBox["*",
"InactiveToken",
BaseStyle->"Inactive",
SyntaxForm->"*"]\)1.`]/.subHJ/.n->1
Activate@Ydd[1\!\(\*
TagBox["*",
"InactiveToken",
BaseStyle->"Inactive",
SyntaxForm->"*"]\)1.`,2.`\!\(\*
TagBox["*",
"InactiveToken",
BaseStyle->"Inactive",
SyntaxForm->"*"]\)1.`]/.subHJ/.n->1


Test subtracting the frequency response of the capsule . 
I'll assume the background is water, that there is Sunflower oil
in the tube, and the tube itself is silica;

Plot[{Abs@Ncore/.ko-> k/.n-> 0,  Abs@tubeSubtract/.ko-> k/.n-> 0},{k,0,1.}, 
PlotLegends->{"Just core", "Capsule - Tube"}, AxesLabel->{"k x tube radius", "response"} ]
Legended[Graphics[{{{}, {}, Annotation[{RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["


(* Rationale:  

Let u[c,t] be the reponse of core + tube (the capsule), u[0,t] be the reponse of just the tube,
and u[c,0] be the reponse of just the core, so that u[0,0] =0 
(u[0,0] the response with the tube and core being the same as the background material.

Then 
(1) u[c,t] - u[0,t] = \!\(
\*SubscriptBox[\(\[PartialD]\), \(c\)]\(u[0, t]\)\)c + O[c^2] = \!\(
\*SubscriptBox[\(\[PartialD]\), \(c\)]\(u[0, 0]\)\)c + \!\(
\*SubscriptBox[\(\[PartialD]\), \(ct\)]\(u[0, 0]\)\)c t + O[c t] + O[c^2]
(2) u[c,0] = u[0,0] + \!\(
\*SubscriptBox[\(\[PartialD]\), \(c\)]\(u[0, 0]\)\)c + O[c^2] = \!\(
\*SubscriptBox[\(\[PartialD]\), \(c\)]\(u[0, 0]\)\)c + O[c^2]
using (1) and (2) implies that 
u[c,t] - u[0,t] = O[c t] + O[c^2]
So the two are similar if the material in the tube is similar to 
the background material.  
    

















(* ::Section:: *)
(*Addition Theorems*)


(* ::Input:: *)
(*ClearAll[Graf,R,\[CapitalTheta]]*)
(*R[l_] := Norm[{x,y}-{x[l],y[l]}];*)
(*\[CapitalTheta][l_] := ArcTan@@({x,y}-{x[l],y[l]});*)
(*R[l_,j_] := Norm[{x[j],y[j]}-{x[l],y[l]}];*)
(*\[CapitalTheta][l_,j_] := ArcTan@@({x[j],y[j]}-{x[l],y[l]});*)
(**)
(*Graf[H_,M_]:= H[n,R[l]]E^(I n \[CapitalTheta][l])  - Sum[H[n-m,R[l,j]]E^(I (n-m) \[CapitalTheta][l,j]) BesselJ[m, R[j]]E^(I m \[CapitalTheta][j]),{m,-M,M}];*)
(*error =Abs[ Graf[HankelH1,50]/HankelH1[n,R[l]] ]/.{x[l]-> 0,y[l]-> 0, x[j]-> 2,y[j]-> 0, y->0};*)
(*Plot[error/.{n->3},{x,0.1,1},PlotRange-> All]*)
(**)
(*error =Abs[ Graf[HankelH1,3]/HankelH1[n,R[l]] ]/.{x[l]-> 0,y[l]-> 0, x[j]-> 2,y[j]-> 0, y->0};*)
(*Plot[error/.{n->3},{x,1.,2},PlotRange-> All]*)
(**)


(* ::Input:: *)
(*error =Abs[ Graf[BesselJ,20]/BesselJ[n,R[l]] ]/.{x[l]-> 0,y[l]-> 0, x[j]-> 2,y[j]-> 0, y->0};*)
(*Plot[error/.{n->3},{x,0.3,2},PlotRange-> All]*)
(**)
(**)

(* ::Package:: *)

(* :Name: MITBagModel` *)

(* :Title: Functions for MITBagModel. *)

(* :Author: SerialCore *)

(* :Copyright: (c) 2020, SerialCore. All rights reserved. *)

(* :Mathematica Version: 12 *)

(* :Package Version: 0.1 *)

(* :Summary: Basic Calculation for MITBagModel. *)



BeginPackage["MITBagModel`"]


Unprotect @@ Names["MITBagModel`*"];
ClearAll @@ Names["MITBagModel`*"];

MassVector::usage="MassVector[] returns mass vactor"
ChargeVector::usage="ChargeVector[] returns charge vactor"
BagConstant::usage="BagConstant[] return bag constant"
ZPEConstant::usage="ZPEConstant[] return zpe constant"

Hadron::usage="Hadron[structure, colorspin, correction] returns mass, radius, x"
HadronFitting::usage="HadronFitting[structure, colorspin, m, B, Z, correction] returns mass, radius, x"

CouplingConstant::usage="CouplingConstant[radius] returns coupling cosntant"
Charge::usage="Charge[quark and antiquark structure] returns charge"
ChargeRadius::usage="ChargeRadius[quark and antiquark structure, hadron] returns charge radius"
MagneticMoment::usage="MagneticMoment[quark and antiquark structure, hadron] returns magnetic moment"

Begin["`Private`"]


MassVector[]:=Return[{5.093,1.641,0.279,0,0}]

ChargeVector[]:=Return[{-1/3,2/3,-1/3,2/3,-1/3}]

BagConstant[]:=Return[0.145^4]

ZPEConstant[]:=Return[1.83]

Hadron[structure_,colorspin_,correction_]:=HadronFitting[structure,colorspin,MassVector[],BagConstant[],ZPEConstant[],correction]

HadronFitting[structure_,colorspin_,mass_,bagconstant_,zpe_,correction_]:=
	Module[{n=structure,CMI=colorspin,m=mass,B=bagconstant,Z=zpe,cor=correction},
		Clear[x,y,R,\[Omega],\[Mu],A,Sij,M];
		
		(* n: number of b, c, s, n quarks or antiquarks *)
		(* m: mass of b, c, s, n quarks or antiquarks *)
		(* B: bag constant *)
		(* Z: zpe *)
		
		(* momentum of b, c, s, q quarks *)
		x={3.08,2.94,2.48,2.04};
		
		(* functions *)
		\[Omega][mi_,xi_,r_]:=(mi^2+xi^2/r^2)^(1/2);
		\[Mu][mi_,xi_,r_]:=r/6 (4\[Omega][mi,xi,r]r+2mi*r-3)/(2\[Omega][mi,xi,r]r(\[Omega][mi,xi,r]r-1)+mi*r);
		y[xi_]:=xi-Sin[xi]Cos[xi];
		A[xi_,xj_]:=1+(xi*Sin[xi]^2-3/2 y[xi])^-1 (xj*Sin[xj]^2-3/2 y[xj])^-1 (-(3/2)y[xi]*y[xj]-2xi*xj*Sin[xi]^2 Sin[xj]^2+1/2 xi*xj(2xi*SinIntegral[2xi]+2xj*SinIntegral[2xj]-(xi+xj)SinIntegral[2(xi+xj)]-(xi-xj)SinIntegral[2(xi-xj)]));
		M[r_]:=Sum[n[[q]]*\[Omega][m[[q]],x[[q]],r],{q,1,4}]+(4\[Pi]*r^3)/3 B-Z/r-3*CouplingConstant[r]*Sum[CMI[[p,q]]*(\[Mu][m[[p]],x[[p]],r]*\[Mu][m[[q]],x[[q]],r])/r^3 A[x[[p]],x[[q]]],{p,1,4},{q,1,4}];
		
		(* use iteration to calculate the minimal mass with precise momentum x *)
		For[loop=0,loop<5,loop++,
			(* perform variation *)
			minimal = Minimize[{M[R],0<R<7},R];
			x[[1]]=X /.FindRoot[Tan[X]==X/(1-m[[1]]*(R/.Last[minimal])-(m[[1]]^2 (R/.Last[minimal])^2+X^2)^(1/2)),{X,2.04,3.14}];
			x[[2]]=X /.FindRoot[Tan[X]==X/(1-m[[2]]*(R/.Last[minimal])-(m[[2]]^2 (R/.Last[minimal])^2+X^2)^(1/2)),{X,2.04,3.14}];
			x[[3]]=X /.FindRoot[Tan[X]==X/(1-m[[3]]*(R/.Last[minimal])-(m[[3]]^2 (R/.Last[minimal])^2+X^2)^(1/2)),{X,2.04,3.14}];
			x[[4]]=X /.FindRoot[Tan[X]==X/(1-m[[4]]*(R/.Last[minimal])-(m[[4]]^2 (R/.Last[minimal])^2+X^2)^(1/2)),{X,2.04,3.14}]
		];

		(* the final result of mass, x, radius, and coupling constant *)
		{First[minimal]+cor,R/.Last[minimal],x}
	]
	
CouplingConstant[radius_]:=
	Module[{r=radius},
		Return[0.296/Log[1+1/(r*0.281)]]
	]
	
Charge[qstructure_,antistructure_]:=
	Module[{qn=qstructure,an=antistructure},
		Q=ChargeVector[];
		Sum[qn[[q]]*Q[[q]]-an[[q]]*Q[[q]],{q,1,5}]
	]
	
ChargeRadius[qstructure_,antistructure_,hadron_]:=
	Module[{qn=qstructure,an=antistructure,hd=hadron},
		m=MassVector[];
		Q=ChargeVector[];
		R=hd[[2]];
		x={hd[[3,1]],hd[[3,2]],hd[[3,3]],hd[[3,4]],hd[[3,4]]};
		
		\[Omega][mi_,xi_,r_]=(mi^2+xi^2/r^2)^(1/2);
		rsquare[mi_,xi_,r_]=r^2 ((\[Omega][mi,xi,r]r(2xi^2 (\[Omega][mi,xi,r]r-1)+4\[Omega][mi,xi,r]r+2mi*r-3)-3/2 mi*r(4\[Omega][mi,xi,r]r+2mi*r-2xi^2-3))/(3xi^2 (2\[Omega][mi,xi,r]r(\[Omega][mi,xi,r]r-1)+mi*r)));
		Abs[Sum[qn[[q]]*Q[[q]]*rsquare[m[[q]],x[[q]],R]-an[[q]]*Q[[q]]*rsquare[m[[q]],x[[q]],R],{q,1,5}]]^(1/2)/5.067
	]
	
MagneticMoment[qstructure_,antistructure_,hadron_]:=
	Module[{qn=qstructure,an=antistructure,hd=hadron},
		m=MassVector[];
		Q=ChargeVector[];
		R=hd[[2]];
		x={hd[[3,1]],hd[[3,2]],hd[[3,3]],hd[[3,4]],hd[[3,4]]};
		
		\[Omega][mi_,xi_,r_]=(mi^2+xi^2/r^2)^(1/2);
		\[Mu][mi_,xi_,r_]=r/6 (4\[Omega][mi,xi,r]r+2mi*r-3)/(2\[Omega][mi,xi,r]r(\[Omega][mi,xi,r]r-1)+mi*r);
		Sum[qn[[q]]*Q[[q]]*\[Mu][m[[q]],x[[q]],R]-an[[q]]*Q[[q]]*\[Mu][m[[q]],x[[q]],R],{q,1,5}]
	]



End[]


EndPackage[]

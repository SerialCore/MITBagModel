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

MassVector::usage="MassVector[] return mass vactor"
ChargeVector::usage="ChargeVector[] return charge vactor"
Flavour::usage="Flavour[structure] returns only flavour configuration"

Hadron::usage="Hadron[structure, colorspin, correction] returns mass, radius"
HadronFitting::usage="HadronFitting[structure, colorspin, m, B, Z, correction] returns mass, radius"

CouplingConstant::usage="CouplingConstant[radius] returns coupling cosntant"
ChargeRadius::usage="ChargeRadius[quark and antiquark structure, hadron] return charge radius"
MagneticMoment::usage="MagneticMoment[quark and antiquark structure, hadron] return magnetic moment"

Begin["`Private`"]


(* ::InheritFromParent:: *)
(*"ChargeRadius[quark and antiquark structure, x, hadron] return charge radius"*)


(* ::InheritFromParent:: *)
(*"MagneticMoment[quark and antiquark structure, hadron] return magnetic moment"*)


MassVector[]:=Return[{5.090,1.640,0.279,0,0}]

ChargeVector[]:=Return[{-1/3,2/3,-1/3,2/3,-1/3}]

Flavour[structure_]:=
	Module[{n=structure},
		f={0,0,0,1};
		For[i=1,i<=3,i++,
			If[n[[i]]==0,f[[i]]=0,f[[i]]=1]
		];
		Return[f]
	]

Hadron[structure_,colorspin_,correct_]:=HadronFitting[structure,colorspin,MassVector[],{0,0,0,0.145^4},{0,0,0,1.83},correct]

HadronFitting[structure_,colorspin_,quarks_,bag_,zpe_,correction_]:=
	Module[{n=structure,CMI=colorspin,m=quarks,B=bag,Z=zpe,cor=correction},
		Clear[x,R];
		
		(* n: number of b, c, s, q quarks *)
		(* m: mass of b, c, s, q quarks *)
		(* B: bag constant of b, c, s, q quarks *)
		(* Z: zpe of b, c, s, q quarks *)
		
		(* flavour configuration of quarks *)
		f=Flavour[n];
		
		(* momentum of b, c, s, q quarks *)
		x={3.08,2.94,2.48,2.04};
		
		(* use iteration to calculate the minimal mass with precise momentum x *)
		For[i=0,i<3,i++,
			(* define M(R) with new x every time *)
			\[Omega][mi_,xi_,R_]=(mi^2+xi^2/R^2)^(1/2);
			\[Mu][mi_,xi_,R_]=R/6 (4\[Omega][mi,xi,R]*R+2mi*R-3)/(2\[Omega][mi,xi,R]*R(\[Omega][mi,xi,R]*R-1)+mi*R);
			y[xi_]=xi-Sin[xi]Cos[xi];
			A[xi_,xj_,R_]=1+(xi*Sin[xi]^2-3/2 y[xi])^-1 (xj*Sin[xj]^2-3/2 y[xj])^-1 (-(3/2)y[xi]*y[xj]-2xi*xj*Sin[xi]^2 Sin[xj]^2+1/2 xi*xj(2xi*SinIntegral[2xi]+2xj*SinIntegral[2xj]-(xi+xj)SinIntegral[2(xi+xj)]-(xi-xj)SinIntegral[2(xi-xj)]));
			M[r_]=Sum[Abs[n[[q]]]*\[Omega][m[[q]],x[[q]],r]+(4\[Pi]*r^3)/3 f[[q]]*B[[q]]-(f[[q]]*Z[[q]])/r,{q,1,4}]-3*CouplingConstant[r]*Sum[CMI[[p,q]]*(\[Mu][m[[p]],x[[p]],r]*\[Mu][m[[q]],x[[q]],r])/r^3 A[x[[p]],x[[q]],r],{p,1,4},{q,1,4}]+cor;
			
			(* perform variation *)
			minimal = Minimize[{M[R],0<R<6},R];
			x[[1]]=X /.FindRoot[Tan[X]==X/(1-m[[1]]*(R/.Last[minimal])-(m[[1]]^2 (R/.Last[minimal])^2+X^2)^(1/2)),{X,2.04,3.14}];
			x[[2]]=X /.FindRoot[Tan[X]==X/(1-m[[2]]*(R/.Last[minimal])-(m[[2]]^2 (R/.Last[minimal])^2+X^2)^(1/2)),{X,2.04,3.14}];
			x[[3]]=X /.FindRoot[Tan[X]==X/(1-m[[3]]*(R/.Last[minimal])-(m[[3]]^2 (R/.Last[minimal])^2+X^2)^(1/2)),{X,2.04,3.14}];
			x[[4]]=X /.FindRoot[Tan[X]==X/(1-m[[4]]*(R/.Last[minimal])-(m[[4]]^2 (R/.Last[minimal])^2+X^2)^(1/2)),{X,2.04,3.14}]
		];

		(* the final result of mass, x, radius, and coupling constant *)
		{First[minimal],R/. Last[minimal],x}
	]
	
CouplingConstant[radius_]:=
	Module[{R=radius},
		Return[0.296/Log[1+1/(R*0.281)]]
	]
	
ChargeRadius[qstructure_,antistructure_,hadron_]:=
	Module[{qn=qstructure,an=antistructure,hd=hadron},
		m=MassVector[];
		Q=ChargeVector[];
		R=hd[[2]];
		x={hd[[3,1]],hd[[3,2]],hd[[3,3]],hd[[3,4]],hd[[3,4]]};
		
		\[Omega][mi_,xi_,R_]=(mi^2+xi^2/R^2)^(1/2);
		r2[mi_,xi_,R_]=R^2 ((\[Omega][mi,xi,R]R(2xi^2 (\[Omega][mi,xi,R]R-1)+4\[Omega][mi,xi,R]R+2mi*R-3)-3/2 mi*R(4\[Omega][mi,xi,R]R+2mi*R-2xi^2-3))/(3xi^2 (2\[Omega][mi,xi,R]R(\[Omega][mi,xi,R]R-1)+mi*R)));
		Abs[Sum[qn[[q]]*Q[[q]]*r2[m[[q]],x[[q]],R]-an[[q]]*Q[[q]]*r2[m[[q]],x[[q]],R],{q,1,5}]]^(1/2)/5.067
	]
	
MagneticMoment[qstructure_,antistructure_,hadron_]:=
	Module[{qn=qstructure,an=antistructure,hd=hadron},
		m=MassVector[];
		Q=ChargeVector[];
		R=hd[[2]];
		x={hd[[3,1]],hd[[3,2]],hd[[3,3]],hd[[3,4]],hd[[3,4]]};
		
		\[Omega][mi_,xi_,R_]=(mi^2+xi^2/R^2)^(1/2);
		\[Mu][mi_,xi_,R_]=R/6 (4\[Omega][mi,xi,R]R+2mi*R-3)/(2\[Omega][mi,xi,R]R(\[Omega][mi,xi,R]R-1)+mi*R);
		Sum[qn[[q]]*Q[[q]]*\[Mu][m[[q]],x[[q]],R]-an[[q]]*Q[[q]]*\[Mu][m[[q]],x[[q]],R],{q,1,5}]
	]



End[]


EndPackage[]

(* ::Package:: *)

BeginPackage["FourPointSpectrum`"]

intervals::usage = 
" Generate intervals for a given beta. Exact values of beta lead to exact solutions. 

Also accepts optional argument epsilon (cutoff for digonal values, default value:10^-4) "

image::usage = 
" Plot intervals. "

points::usage =
" Plot the endpoints of the intervals. "

Begin["`Private`"]
SetSystemOptions["CompileOptions" -> "CompileReportExternal"->True]

(* Compiled helper functions. *)
cC1 = Compile[{{x, _Real, 1}, {b, _Real}}, b^(x+1.0)/(1.0-b),
CompilationOptions -> {"InlineCompiledFunctions" -> True, "ExpressionOptimization" -> True, "InlineExternalDefinitions" -> True},
CompilationTarget -> "C", RuntimeOptions -> "Speed"]; 
cC2= Compile[{{x, _Real, 1}, {b, _Real}}, x+(b^(x+1.0)-b)/(1.0-b),  (* sum to the left at position x>0 *)
CompilationOptions -> {"InlineCompiledFunctions" -> True, "ExpressionOptimization" -> True, "InlineExternalDefinitions" -> True},
CompilationTarget -> "C", RuntimeOptions -> "Speed"];
cI1 = Compile[{{a, _Real, 1}, {b, _Real}}, Ceiling[Log[a]/Log[b]]-1, 
CompilationOptions -> {"InlineCompiledFunctions" -> True, "ExpressionOptimization" -> True, "InlineExternalDefinitions" -> True},
CompilationTarget -> "C", RuntimeOptions -> "Speed"]; 
gcd = Compile[{{aa, _Integer}, {bb, _Integer}}, Module[{t, b=bb, a=aa}, While [b != 0, t = b; b = Mod[a, b]; a = t;]; a],  
CompilationOptions -> {"InlineCompiledFunctions" -> True, "ExpressionOptimization" -> True, "InlineExternalDefinitions" -> True},
CompilationTarget -> "C", RuntimeOptions -> "Speed"];

(* Symbolic helper functions. *)
intervals[b_?MachineNumberQ, eps_: 0.0001] := ints[b, eps]

C1[x_, c_] := c^(x+1)/(1-c);  
C2[x_, c_] := x+(c^(x+1)-c)/(1-c);  
I1[a_, c_] := Ceiling[Log[a]/Log[c]]-1;

(* The main symbolic module. *)
intervals[b_, eps_: 10^-4] := 
	Module[{t, endpoints, vI1, vI2, max, min, lend, cend, mid, vCc, vMc, i, k, M, A, a},
		(* endpoints of intervals on which Cc is constant *)
		M = I1[N[eps], N[b]];
		t = Table[b^k, {k, 1, M}];
		AppendTo[t, eps];
		(* Sort must be done according to numerical value. *)
		endpoints = SortBy[Join[1-t, t], N];
		(* values of Cc and Mc on the intervals *)
		(* just check midpoints, since constant inside *)
		mid = (endpoints[[;;-2]] + endpoints[[2;;]])/2;
		vI1 = I1[#, b]& /@ mid;
		vI2 = Reverse[vI1];
		vCc = C1[vI1, b] + C2[vI2, b];
		vMc = vI2-vI1;
		(* max values of x and min values of y, given a fixed slope *)
		i = 1;
		k = Length[endpoints]-1;
		max = Reap[While[k>0,
				If[vCc[[k]]/(i+vMc[[k]]) > endpoints[[k]],
					Sow[vCc[[k]]/(i+vMc[[k]])]; i++,
					k--
				]
			  ]]//Last//Last;
		M = Length[max];
		min = 1-max;
		(* lengths of intervals *)
		(* from the side *)
		lend = Table[g=GCD[N2, m]; A=m-N2 min[[m]]; 
			If[A<=0, Unevaluated@Sequence[],
			Table[
			If[
			    GCD[g, N1] > 1,
			    Unevaluated@Sequence[],
			    a = Min[max[[N1+N2-m]], m/(N1+N2), A/N1]; {{0, m/N2}, {a, (m-N1 a)/N2}}
			], {N1, 1, M-N2}]],
			{N2, 1, M-1}, {m, 1, N2-1}] // Flatten[#, 2]&;
		(* from vertex *)
		cend = Table[
			If[
			    GCD[N1, N2] > 1,
			    Unevaluated@Sequence[],
			    a = Min[max[[N1]], N2/(N1+N2), N2 max[[N2]]/N1]; {{0, 1},{a, 1-N1 a/N2}}
			],
			{N1, M}, {N2, M-N1}] // Flatten[#, 1]&;
		(* starting points vertical, horizontal and vertex *)
		lend = Join[lend, 1-lend, cend];
 		(* rend is a reflection of lend about y=x *)
		lend = Join[lend, RotateLeft[lend, {0,0,1}]];
		lend
];

(* The main numerical, compiled module. *)
ints = With[{I1=cI1, C1=cC1, C2=cC2, gcd=gcd},
        Compile[{{b, _Real}, {eps, _Real}}, 
    	    Module[
		{t, endpoints, vI1, vI2, max, bag, bag2, min, lend, cend, s, mid, vCc, vMc, i, k, M, temp, N1, N2, a=0.0, m, g, lst, A=0.0},
		(* endpoints of intervals on which Cc is constant *)
		M = I1[{eps}, b][[1]];
		t = Table[b^k, {k, 1, M}];
		AppendTo[t, eps];
		endpoints = Sort[Join[1.0-t, t]];
		(* values of Cc and Mc on the intervals *)
		(* just check midpoints, since constant inside *)
		mid = (endpoints[[1;;-2]] + endpoints[[2;;-1]])/2.0;
		vI1 = I1[mid, b];
		vI2 = Reverse[vI1];
		vCc = C1[vI1, b] + C2[vI2, b];
		vMc = vI2-vI1;
		(* max values of x and min values of y, given a fixed slope *)
		i = 1;
		k = Length[endpoints]-1;
		bag = Internal`Bag[];
		While[k>0,
		    	temp = vCc[[k]]/(i+vMc[[k]]);
			If[temp > endpoints[[k]],
			    Internal`StuffBag[bag, temp]; i++,
			    k--
		    	]
		];
		max = Internal`BagPart[bag, All];
		M = Length[max];
		min = 1-max;
		(* lengths of intervals *)
		(* from the side *)
		bag = Internal`Bag[]; 
		bag2 = Internal`Bag[];
		Do[
    		  Do[
    		    g = gcd[N2, m];
		    A = m-N2 min[[m]];
		    If[A>0,
		      Do[
			If[g==1 || gcd[N1, g]==1,
			    a = Min[max[[N1+N2-m]], N[m]/(N1+N2), A/N1];
			    Internal`StuffBag[bag, {N[m]/N2, a, (m-N1 a)/N2}, 1];
			],
		      {N1, 1, M-N2}]],
		  {m, 1, N2-1}];
		    Do[
		    If[gcd[N1, N2]==1, 
			a = Min[max[[N1]], N[N2]/(N1+N2), (N2 max[[N2]])/N1];
			Internal`StuffBag[bag2, {a, 1.0 - N1 a/N2}, 1]];,
		    {N1, 1, M-N2}],
		{N2, 1, M}];     
		lst = Internal`BagPart[bag, All];
		bag = Internal`Bag[{0.0}];
		If[Length[lst]>0,
			lend = {{0.0, #[[1]]}, {#[[2]], #[[3]]}}& /@ Partition[lst, 3],
			lend = Most[{{{0.0,0.0},{0.0,0.0}}}]
		];
		lst = Internal`BagPart[bag2, All];
		bag2 = Internal`Bag[{0.0}];
		If[Length[lst]>0,
			cend = {{0.0, 1.0}, #}& /@ Partition[lst, 2],
			cend = Most[{{{0.0,0.0},{0.0,0.0}}}]
		];
		lend = Join[lend, 1.0-lend, cend];
		lend = Join[lend, RotateLeft[lend, {0,0,1}]];
		lend
],
CompilationTarget -> "C", RuntimeOptions -> "Speed",
CompilationOptions -> {"InlineCompiledFunctions" -> True, "ExpressionOptimization" -> True, "InlineExternalDefinitions" -> True}
]];

diag = Plot[x, {x,0,1}, PlotStyle->Dashed, PlotRange->{{0,1},{0,1}}, AspectRatio->1, Frame -> True]; 
image[b_, eps_:0.0001] := Labeled[Show[diag, 
			 Graphics[{AbsoluteThickness[1-b], Line[ints[N[b], N[eps]]]}, 
			 	  PlotRange->{{0,1},{0,1}}, AspectRatio->1, AntiAliasing->False], 
			 ImageSize->400], b]
points[b_, eps_:0.0001] := Labeled[Show[diag, 
			 Graphics[{Point[ints[N[b], N[eps]][[All,2]]]}, 
			 	  PlotRange->{{0,1},{0,1}}, AspectRatio->1], 
			 ImageSize->400], b]

End[]

EndPackage[]


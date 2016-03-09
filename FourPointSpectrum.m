(* ::Package:: *)

BeginPackage["FourPointSpectrum`"]

image::usage = 
"A few ways to use the function image[beta]:

- Generate an image.
b = 95/100;
im = image[b] // AbsoluteTiming;
Export[\" image.png \", im[[2, 1]], ImageSize -> 2000]

- Export a movie. Saving avi takes longer than generating the data:
B = 0.80;
animate=Monitor[Table[image[b],
{b,0.01,B,0.002}],Row[{ProgressIndicator[b,{1/100,B}],\"b=\",N[b]}]] // AbsoluteTiming;
animate[[1]]
Export[\"~/movie.avi\",animate[[2,All,1]],ImageSize->800]

- Manipulate beta almost in real time:
B=0.8
Manipulate[image[b],{b,0.01,0.8}]"

Begin["`Private`"]


Print["Use image[beta] to generate the set of possible eigenvalues E1 and E2 for a class of self-adjoint operators (parametrized by beta), with the point spectrum equal to 0=E0<E1<E2<E3=1. "]

(* Clean version of line endpoints with arbitrary beta *) 
(* all multiplicities considered, based on epsilon *)
eps=10^(-4); (* accuracy of computations *) 

(* PART 1: functions C, D,  m for symmetric geomtric sequence *)
C1[x_, b_] := b^(x+1)/(1-b);  (* sum to the left at position x<0 *)
C2[x_, b_] := x+(b^(x+1)-b)/(1-b);  (* sum to the left at position x>0 *)
I1[a_, b_] := Ceiling[Log[a]/Log[b]]-1; (* nearest position to a, a<1/2 *)
I2[a_, b_] := Ceiling[Log[1-a]/Log[b]]-1; (* nearest position to a, a>1/2 *)
Cc[a_, b_] := C1[I1[a, b], b] + C2[I2[a, b], b];(* Sum C(a) *) 
Mc[a_, b_] := I2[a, b] - I1[a, b]; (* Counting elements function m(a) *)


(* Remove all N[] calls to get rational results. The performance remains very good. *)
diag = Plot[x, {x,0,1}, PlotStyle->Dashed, PlotRange->{{0,1},{0,1}}, AspectRatio->1, Frame -> True]; 
image[b_] := 
	Block[{t, endpoints, int, vI1, vI2, max, min, lend, cend, s, time, mid, vCc, vMc, i, k, M, a, tbl},
		(*time = AbsoluteTime[];*)
		(* endpoints of intervals on which Cc is constant *)
		M = Quiet[FullSimplify[I1[eps, b]]];
		t = Table[b^k, {k, 1, M}];
		AppendTo[t, eps];
		endpoints = Sort[Join[1-t, t]];
		(* values of Cc and Mc on the intervals *)
		(* just check midpoints, since constant inside *)
		mid = (endpoints[[;;-2]] + endpoints[[2;;]])/2;
		vI1 = Map[I1[#, b]&, mid];
		vI2 = Reverse[vI1];
		vCc = C1[vI1, b] + C2[vI2, b];
		vMc = vI2-vI1;
		(* max values of x and min values of y, given a fixed slope *)
		i = 1;
		k = Length[endpoints]-1;
		max = Reap[While[k>0,
				If[vCc[[k]]/(i+vMc[[k]]) > endpoints[[k]],
					Sow[vCc[[k]]/(i+vMc[[k]])]; i+=1,
					k-=1
				]
			  ]][[2,1]]//N;
		M = Length[max];
		min = 1-max;
		(*maxtime=AbsoluteTime[]- time;
		time = AbsoluteTime[];*)
		(* lengths of intervals *)
		(* parallel is not always faster *)
		tbl = If[M > 25, DistributeDefinitions[M, max, min, eps]; ParallelTable, Table];
		(* from the side *)
		lend = tbl[
			Which[
			    GCD[N1,N2,m] > 1,
			    Unevaluated@Sequence[],
			    a = Min[max[[N1+N2-m]], N[m]/(N1+N2), (m-N2 min[[m]])/N1]; a>=eps,
			    {{0., N[m]/N2}, {a, (m-N1 a)/N2}},
			    True,
			    Unevaluated@Sequence[]
			],
			{N1, M}, {N2, M-N1}, {m, 1, N2-1}] // Flatten[#, 2]&;
		tbl = If[M > 100, ParallelTable, Table];
		(* from vertex *)
		cend = tbl[
			Which[
			    GCD[N1,N2] > 1,
			    Unevaluated@Sequence[],
			    a = Min[max[[N1]], N[N2]/(N1+N2), N2 max[[N2]]/N1]; a>=eps,
			    {{0.,1.},{a,1.-N1 a/N2}},
			    True,
			    Unevaluated@Sequence[]
			],
			{N1, M}, {N2, M-N1}] // Flatten[#, 1]&;
		(* starting points vertical, horizontal and vertex *)
		lend = Join[lend, 1-lend, cend];
 		(* rend is a reflection of lend about y=x *)
		lend = Join[lend, RotateLeft[lend, {0,0,1}]];
		(*lendtime=AbsoluteTime[]- time;*)
		(* graph of line segements *)
		s = Graphics[{AbsoluteThickness[1-b], Line[lend]}, PlotRange->{{0,1},{0,1}}, AspectRatio->1];
		Labeled[Show[diag, s, ImageSize->400], N[b]]
];

End[]

EndPackage[]


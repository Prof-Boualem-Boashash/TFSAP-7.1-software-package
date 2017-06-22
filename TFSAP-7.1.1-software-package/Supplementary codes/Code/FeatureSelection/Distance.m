function Res = Distance(P,Q)


	Res = Divergence(P,Q,0) + Divergence(Q,P,0);

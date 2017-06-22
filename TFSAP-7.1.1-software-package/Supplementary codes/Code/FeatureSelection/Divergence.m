function Res = Divergence(P,Q,meth)
	

	
	if meth == 0
		P = P+eps;
		Q = Q+eps;

	elseif meth == 1
		indice_P = find(P==0);
		indice_Q = find(Q==0);
		indices = union(indice_P,indice_Q);
		P(indices) = [];
		Q(indices) = [];
	end
	
		P = P/sum(P);
		Q = Q/sum(Q);
	
	Res = sum(P.*(log10(P)-log10(Q)));

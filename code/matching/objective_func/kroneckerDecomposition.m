function [B,C, bVec, cVec] = kroneckerDecomposition(A, Br, Bc, doFull)

	if ( ~exist('doFull', 'var') )
		doFull = 0;
	end
	if ( ~doFull )
		% factorise A = kron(B,C) in the least-squares sense, where size(B) = [Br,Bc]
		Atilde = kron2vec(A, Br, Bc);
		[u,s,v] = svds(Atilde,1);
		
		bVec = u*sqrt(s);
		cVec = v*sqrt(s);
		
        
		B = reshape(bVec, Br, Bc);
		
		Ar = size(A,1)/Br;
		Ac = size(A,2)/Bc;
		
		C = reshape(cVec, Ar, Ac);
	else
		% factorise A = sum_i kron(B_i,C_i) exactly, where size(B) = [Br,Bc]
		Atilde = kron2vec(A, Br, Bc);
        bb = full(Atilde);
		[u,s,v] = svd(full(Atilde));
		
		epsilon = 1e-10;
		sigmaMax = s(1,1);
		idx = find(diag(s)./sigmaMax>epsilon);
		
		rankA = numel(idx);
		
		B = {};
		C = {};
		bVec = {};
		cVec = {};
		for c=1:rankA
			bVec{c} = u(:,c)*sqrt(s(c,c));
			cVec{c} = v(:,c)*sqrt(s(c,c));
			
			B{c} = reshape(bVec{c}, Br, Bc);
			
			Ar = size(A,1)/Br;
			Ac = size(A,2)/Bc;
			
			C{c} = reshape(cVec{c}, Ar, Ac);
		end
		
	end
end
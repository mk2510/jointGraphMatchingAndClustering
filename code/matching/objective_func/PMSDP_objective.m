function [ob obj2 obj3] = PMSDP_objective(K)
	K = -K; % now we can formulate minimisation problems
	
    % added here the * for matrix multiplication
    
    %TODO: SDP notation, the P should probably be a SPD value matrix
	obj = @(P) full(P(:)' *  K * P(:));
	 
	% factorise K
	X = rand(N,N);
	
	
	
	offDiagMask = ones(N,N)-eye(N);
	
	switchAB = 0;
	
	diagAc = {};
	diagBc = {};

	% full kronecker decomposition
	[AcCell,BcCell] = kroneckerDecomposition(K, N, N, 1);
	C = numel(AcCell);
	
	% move out diagonal elements of Ac and Bc into diagAc/diagBc
	AcSum = zeros(N);
	BcSum = zeros(N);
    
% 	figure;
	for c=1:C
		diagAc{c} = diag(AcCell{c});
		diagBc{c} = diag(BcCell{c});
		
		AcCell{c} = offDiagMask.*AcCell{c};
		BcCell{c} = offDiagMask.*BcCell{c};
% 		imagesc(AcCell{c});
% 		drawnow;
% 		pause(0.5);
		AcSum = AcSum + AcCell{c};
		BcSum = BcSum + BcCell{c};
	end

	obj2 = 0;
	obj3 = 0;
	for c=1:C
		obj2 = obj2 + trace(AcCell{c}*X'*BcCell{c}*X) + ...
			diagAc{c}'*X*diagBc{c};
	
		obj3 = obj3 + trace(BcCell{c}*X'*AcCell{c}*X) + ...
			diagAc{c}'*X*diagBc{c};
% 		obj3 = obj3 + trace(AcCell{c}*(X'*BcPlus{c}*X - X'*BcMinus{c}*X));	
    end
    
    ob = obj(X);
	[ob obj2 obj3]
end


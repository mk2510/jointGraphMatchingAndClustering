function Aout = kron2vec(A, Br, Bc)
% [Br,Bc] = size(B);
% for A = kron(B,C): [kron2vec(A) = B(:)*C(:)']
% 	D = size(A,1)/2;
% 	N = D;
	
    % I think this function is wrong, index error?
	% reshape
	%Aout = [];
	%for ii=1:Br
	%	iIdx = (ii-1)*Br+1:ii*Br;
	%	for jj=1:Bc
	%		jIdx = (jj-1)*Bc+1:jj*Bc;
	%		newRow = A(jIdx, iIdx);
	%		Aout = [Aout; newRow(:)'];
	%	end
    %end
	
    Aout = [];
	for ii=1:Bc
		iIdx = (ii-1)*Br+1:ii*Br;
		for jj=1:Br
			jIdx = (jj-1)*Bc+1:jj*Bc;
			newRow = A(jIdx, iIdx);
			Aout = [Aout; newRow(:)'];
		end
	end
    
    
% 	D = size(A,1)/2;
% 	N = D;
% 	
% 	% reshape
% 	Aout = [];
% 	for ii=1:D
% 		iIdx = (ii-1)*N+1:ii*N;
% 		for jj=1:D
% 			jIdx = (jj-1)*N+1:jj*N;
% 			newRow = A(jIdx, iIdx);
% 			Aout = [Aout; newRow(:)'];
% 		end
% 	end
end
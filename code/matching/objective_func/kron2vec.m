% 
% original code by: Florian Bernard (f.bernardpi@gmail.com)
% 

function Aout = kron2vec(A, Br, Bc)
	
    Aout = [];
	for ii=1:Bc
		iIdx = (ii-1)*Br+1:ii*Br;
		for jj=1:Br
			jIdx = (jj-1)*Bc+1:jj*Bc;
			newRow = A(jIdx, iIdx);
			Aout = [Aout; newRow(:)'];
		end
	end
   
end

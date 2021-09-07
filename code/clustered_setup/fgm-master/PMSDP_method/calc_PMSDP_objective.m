function [objY] = calc_PMSDP_objective(Y,probDim,W,normPSquared, normQSquared, X, k, n)
% calculates the objective used in the PMSDP discretisation    
    objY = 0;
    for ii = 1:k
        Y{ii} = sdpvar( n, probDim ^ 2, 'full' );
        objY = objY + 2 * trace( W(:,(ii-1)*n + 1 : ii*n) * Y{ii} );
    end
    objY = objY - normPSquared - normQSquared' * sum(X,2);

end


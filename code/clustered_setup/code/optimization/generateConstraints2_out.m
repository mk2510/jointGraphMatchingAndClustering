function [F,objY,problem] = generateConstraints2_out( problem,params, number, numberOfGroups, groupLs)
%===============================================================
% module:
% ------
% generateConstraints2_out.m
%
% Description:
% -----------
% generates PM-SDP constraints

%===============================================================

%--------------------------------------------
% Initialization
%--------------------------------------------
n = params.n;
k = params.k;
probDim = params.probDim{number};
F_blocks = cell( n );
objY = 0;
% utilize R
utilizeRFlag = params.utilizeRFlag;
Rtol = params.Rtol;
% partial flag
partialNFlag = params.n ~= params.k;
% frobenious normQSquared
normQSquared = diag(problem.Q{number}'*problem.Q{number});
normPSquared = problem.normPSquared{number};
UtilizeXflag = params.UtilizeXflag;
permConstraint = params.permConstraint;
X = problem.X;
R = problem.R{number};
Y = problem.Y{number};
B = problem.B{number};
A = problem.A{number};
W = problem.W{number};
%============================================

%--------------------------------------------
% Prepare rows and cols to remove
%--------------------------------------------
% true means that we want to remove this entry
RrowsColsOffset = 1 + n;
if utilizeRFlag
    RlogicMat = true( probDim );
    RzeroRIdx = (tril(RlogicMat,Rtol) == 0) | (triu(RlogicMat,-Rtol) == 0);
    %RlogicMat = true(probDim);
    %for ki = 1:numberOfGroups
    %    RdiagMat{ki} = true(groupLs(ki));
    %end
    %Rblkdiag = blkdiag(RdiagMat{:});
    %RzeroRIdx = RlogicMat - Rblkdiag;
    %RzeroRIdx = logical(RzeroRIdx);
else
    RzeroRIdx = false( probDim);
end

RzeroRIdxWithOffset = [false(RrowsColsOffset,1); RzeroRIdx(:)];
XIrowsColsOffset = 1;
XzeroIdxWithOffsetCells=cell(1,k);

for ii=1:k
    if UtilizeXflag
        XIzeroRIdx  = ~permConstraint(:,ii);
    else
        XIzeroRIdx = false( n,1 );
    end
    XzeroIdxWithOffsetCells{ii} = [false(XIrowsColsOffset,1); XIzeroRIdx(:); false(probDim^2,1)];
end
%============================================

%--------------------------------------------
% Generate constraints
%--------------------------------------------
% get constraints of ds and rotation
[ A_ds, b_ds ] = getDoublyStochasticConstraints( n, k );
[ A_rot, b_rot ] = getLiftedRotationConstraints( probDim );
% set elements equal zero constraints according to assumptions on R

F_zerosR = (R(RzeroRIdx) == 0);
F_zerosB = (B(:,RzeroRIdx) == 0);

for ii = 1 : k
    %ii
    % current variables
    Y{number, ii} = sdpvar( n, probDim ^ 2, 'full' );
    A{number, ii} = diag(X(:,ii));
    % objective, keep in mind that we are maximizing
    objY = objY + 2 * trace( W(:,(ii-1)*n + 1 : ii*n) * Y{number, ii} );
    % block constraints
    constraintMat = [1 [X(:,ii)'  colStack(R)']; [X(:,ii); colStack(R)]...
        [A{number, ii} Y{number, ii};Y{number, ii}' B]];
    % logical vector for row removel - combine R and XI vectors
    zeroIdxWithOffset = XzeroIdxWithOffsetCells{ii} | RzeroRIdxWithOffset;
    % zero cols from Y
    F_blocks{ii} = (constraintMat(~zeroIdxWithOffset,~zeroIdxWithOffset) >= 0)...
        + (ones(1,n) * Y{number, ii} == colStack(R)') + (Y{number, ii}(:,RzeroRIdx) == 0);
    
end

% general constraints- ds and rotation in the lift
% partial
if partialNFlag
    F_general = (A_ds(1:n,:) * X(:) <= b_ds(1:n)) + (A_ds((n + 1) : end,:) * X(:) == b_ds((n + 1) : end)) + (X(:) >= 0) ...
        + (A_rot * colStack(B) == b_rot);
else
    F_general = (A_ds * X(:) == b_ds) + (X(:) >= 0) + (A_rot * colStack(B) == b_rot);
end

% use AGD information for reducing problem size
if UtilizeXflag
    % set elements equal zero constraints according to assumptions on X using
    % AGD
    F_Xi_AGD=[];
    F_Ai_AGD=[];
    F_Yi_AGD=[];
    
    for ii = 1:k
        % constraints on X
        Xi=X(:,ii);
        badIdx=~permConstraint(:,ii);
        if sum(badIdx)>0
            F_Xi_AGD = [F_Xi_AGD Xi(badIdx)==0];
            % constraints on Ai
            F_Ai_AGD = [F_Ai_AGD A{number, ii}(badIdx,:)==0];
            F_Ai_AGD = [F_Ai_AGD A{number, ii}(:,badIdx)==0];
            % constraints on Yi
            F_Yi_AGD = [F_Yi_AGD Y{number, ii}(badIdx,:)==0];
        end
    end
    % add constraints to general constraints
    F_general=F_general + F_Xi_AGD + F_Ai_AGD + F_Yi_AGD;
end

% gather all constraints
F = [F_blocks{:}] + F_general + F_zerosR + F_zerosB;
%============================================

%--------------------------------------------
% Final objective
%--------------------------------------------
objY = objY - normPSquared - normQSquared' * sum(X,2);
%============================================
problem.Y = Y;
problem.A = A;

end


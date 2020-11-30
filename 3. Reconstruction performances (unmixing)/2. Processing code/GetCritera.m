function [M,x,MSE,SAD,aSAD,aMapsMSE,Maps_est] = GetCritera(Xhat,Y,I)

    load('../../synth_data/synthetic_data.mat','X','Spectra','Maps')
    load('../../reproducible_results/synthetic_analysis.mat','FourPixels')

    [n,m,l] = size(Xhat);

    % Resized maps
    x0 = transpose(reshape(cat(3,Maps{:}),[n*m 4]));
    M0 = Spectra;

    %% MSE
    MSE = norm(X(:)-Xhat(:),'fro')^2/(norm(X(:),'fro')^2);


    %% Spectra estimation

    % Estimates the unmixed spectra using always the same 4 pixels and using
    % acquired data
    Y = transpose(reshape(Y,[n*m l]));
    init = Y(:,FourPixels); 

    M = sisal(Y(:,I),4,'verbose',0,'M0',init); 

    % to be sure that the spectra are +
    M = M.*repmat(1-2*(sum(M>0, 1) < 10), [size(M,1) 1]);

    %% Maps Estimation

    x = sunsal(M,transpose(reshape(Xhat,[n*m l])),'POSITIVITY', 'yes', 'ADDONE', 'yes');

    %% Identification between original and estimated maps and spectra

    % Get matching permutation vector
    tab_M = GetMatching(x0,x);

    % Permute
    M = M(:,tab_M); x = x(tab_M,:);
                
    %% Compute metrics            
    SAD = zeros(1,4);
    MapsMSE = zeros(1,4);

    %SAD
    SAD = acos( diag(M0'*M)./( sqrt(diag(M0'*M0)).*sqrt(diag(M'*M)) ) );
    aSAD = mean(SAD);

    % Map MSE
    aMapsMSE = norm(x-x0,'fro')/norm(x0,'fro');

    %% Output maps
    Maps_est = cell(1,4);
    for ind = 1:4
        Maps_est{ind} = reshape(x(ind,:),[n m]);
    end

end




function [ tab_M ] = GetMatching( M1,M2 )
    %Compute the cross correlation beetween two (M,N) matrices, i.e. the
    %(M,M) matrix whose (i,j) element is Corr(M1(i,:),M2(j,:)). A M-length permutation
    %vector P is then deduced such as M1 matches with M2(P,:).

    M = size(M1,1);

    % Correlation operator
    corrx = @(x,y) abs(x'*y)/(norm(x)*norm(y));

    % Computation
    Corr = zeros(M);
    for ind_x = 1:M
        for ind_y = 1:M
            Corr(ind_x,ind_y) = corrx(transpose(M1(ind_x,:)),transpose(M2(ind_y,:)));
        end
    end

    % Find matching beetween M1 and M2
    Corr2 = Corr;
    tab_M = zeros(1,M);
    for ind = 1:M
        [Max,ind_Max] = max2(Corr2);
        tab_M(ind_Max(1)) = ind_Max(2);
        Corr2(ind_Max(1),:) = 0; Corr2(:,ind_Max(2)) = 0;
    end

end


function [ m , indices ] = max2( Mat )
    %Returns the maximum of a n dimensionnal matrix and give its position.
    %
    % [ m , indices ] = minn( Mat );
    %
    %Input: 
    %   Mat: The n-dimensionnal matrix to study
    %
    %Output:
    %   m: the maximum of the whole Mat matrix
    %   indices: the position of the maximum value in Mat
    %   

    if nargin ~= 1
        error('Invalid number of argument : one matrix required')
    end

    % Dim and init
    N = ndims(Mat);    
    M = cell(1,N+1); I = cell(1,N);
    M{N+1} = Mat;

    % Makes all min
    for ind = N:-1:1
        [M{ind} I{ind}] = max(M{ind+1},[],ind);
    end

    % Find outputs
    m = M{1};

    indices = zeros(1,N);
    indices(1) = I{1};

    for ind = 2:N
        indices(ind) = I{ind}(I{ind-1});
    end

end



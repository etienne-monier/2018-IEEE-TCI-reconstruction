function [Spectra, Maps, x] = perform_unmixing(Xhat,Y,I,x0)

    FirstFlag = ~exist('x0','var');
    Norm = @(x) (x-min(x(:)))/(max(x(:))-min(x(:)));

    FourPixels = [11 304 688 1333 2272];

    [n,m,l] = size(Y);

    %% Spectra estimation

    % Estimates the unmixed spectra using always the same 4 pixels and using
    % acquired data
    Y = transpose(reshape(Y,[n*m l]));

    init = Y(:,FourPixels); 
    M = sisal(Y(:,I),5,'verbose',0,'M0',init); 

    % To be sure spectra are +
    M = M.*repmat(1-2*(sum(M>0,1)<size(M,1)*0.5),[size(M,1) 1]);

    %% Maps Estimation
    x = sunsal(M,transpose(reshape(Xhat,[n*m l])),'POSITIVITY', 'yes', 'ADDONE', 'yes');

    %% Matching between first data and these estimated maps and spectra

    if (FirstFlag) % first study, sort the results according to their power
        Pow = mean(x.^2,2);
        [Pow2, tab_M] = sort(Pow);
    else
        tab_M = GetMatching(Norm(x0),Norm(x));
    end
    M = M(:,tab_M); x = x(tab_M,:);

    %% Output maps
    Maps = cell(1,5);
    for ind = 1:5, Maps{ind} = reshape(x(ind,:),[n m]); end
    Spectra = M;

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

    if nargin ~=1
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


function [ H,S ] = PcaDimensionReduction( Y, M, I )
    %Computes the PCA of Y (whose size is (m,n) and keeps only the M first
    %components. The output are a (m,M) matrix of principal components vectors
    %and a (M,n) abundance matrix.
    %   Input:
    %       Y: (m,n) data matrix
    %       M: desired output data dimension
    %       I: sampled pixel table
    %   Output:
    %       H: (m,M) matrix of principal components vectors
    %       S: (M,n) of principal components abundance maps
    [m n] = size(Y);

    N = length(I);
    if(M>(N-1) || M>m)
        error(['Wrong choice of dimension to keep. It should be less than ' num2str(min(m,N-1))])
    end

    [V,D] = eig(cov(transpose(Y(:,I))));
    V = fliplr(V); d = flipud(diag(D));

    H = V(:,1:M); S = H'*Y;

end

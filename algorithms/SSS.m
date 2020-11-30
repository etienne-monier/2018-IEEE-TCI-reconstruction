function [ Xhat ] = SSS( Y, I, lambda, varargin )
% Reconstructs a full multiband image using the Smoothed SubSpace (3S) 
%  algorithm. For further details, see "Reconstruction of partially sampled
%  multi-band images -- Application to EELS microscopy,"
%  (http://monier.perso.enseeiht.fr/publications.php)
%
%   Input:
%       Y: the acquired data of size (n, m, l) where n and m are the spatial
%           dimensions and l is the number of channels. Y should be filled with 0 at
%           non-sampled pixels,
%       I: Sampling path. The elements of I are the indices of the acquired pixels in
%           column-major order. 
%       lambda: spectral regularization parameter,
%
%   Optional:
%       'init': sets the initial point for FISTA. Should be a (l,n*m)
%           matrix,
%       'verbose': set the display option by choosing 'y' whereas 'n'
%           option displays nothing.
%
%   Output:
%       Xhat: the reconstructed image under a (n,m,l) matrix form.
%
% -------------------------------------------------------
%
% This code is distributed under the terms of the GNU General Public License 3.0.
%
% ----------------------------------------------------------------------
% 

%% Getting input

defaultInitPoint = randn(size(Y,3),size(Y,1)*size(Y,2));
defaultverbose = 'y';
expectedShapes = {'y','n'};

p = inputParser;
addRequired(p,'Y');
addRequired(p,'I');
addRequired(p,'lambda');

addParameter(p,'init',defaultInitPoint, @(x) isnumeric(x)); % size(x) == size(Y) );
addParameter(p,'verbose',defaultverbose,  @(x) any(validatestring(x,expectedShapes)) );

parse(p, Y, I, lambda, varargin{:});

if (p.Results.verbose == 'y')
    verb = true;
else
    verb = false;
end

InitPoint   = p.Results.init;


%% Size and path

Ns = length(I);
[nr,nc,Nb] = size(Y);
Np = nr*nc;

%% Remove mean

Y = transpose(reshape(Y,[Np Nb]));

A = zeros(Nb,Np); A(:,I)=1;

Ym = repmat(mean(Y(:,I),2),[1 Np]);
Y = A.*(Y-Ym);

%% Performing PCA
if (verb)
    disp(sprintf('\n----------- PCA Study -----------\n'))
    disp('')
end

[V,D] = eig(cov(Y(:,I)'));
d = flipud(diag(D));
V = fliplr(V);

if (Ns<=Nb)
    d = d(1:(Ns-1));
    V = V(:,1:(Ns-1));
end

[ dout, sigma, RhoMax ] = EigenEstimate( d, Ns );

if (verb)
    disp(['Chosen sigma^2 is ' num2str(sigma^2)])
    disp(['Chosen signal subspace dimension is ' num2str(RhoMax)])
end

%% Sparse operators

temp = mat2cell( repmat( subarray( -eye(nr) + diag(1*ones(1,nr-1),-1),1:nr,1:(nr-1)) ,[1 nc] ) ,nr,repmat(nr-1,[1 nc]));
NablaX = spdiags([ones(nr*(nc-1),1), -ones(nr*(nc-1),1)],[-nr 0],nr*nc,nr*(nc-1));
NablaY = sparse(blkdiag( temp{:} ));
Nabla = [NablaY NablaX];
Delta = -Nabla*Nabla';

%% Reconstruction

if (verb)
    disp('')
    disp(sprintf('\n----------- Reconstruction -----------\n'))
    disp('')
end

step = @(lambda,W,rho) 1/(8/rho + lambda*max(diag(W)) );
Gradient = @(Sy, lambda, W, rho) -Sy*Delta/rho + lambda*W*Sy;
cost = @(S,lambda,W, rho) (1/(2*rho)) * ( norm( S*Nabla , 'fro')^2 + lambda*norm(W.^(1/2)*S , 'fro')^2   );

%-- Stopping conditions
ConvE = @(E,step,n) abs(E(n)-E(n-1))/(E(n-1)*step);

rho = RhoMax;
           
% Knowledge
H = V(:,1:rho); 
Yu = H'*Y;

%-- The initial point
S0 = InitPoint(1:Nb,1:Np);
Nit_MAX = 2500;

% Weight computation needs rho
epsilon=5e-6;
W = diag( sigma^2./max(d(1:rho)-sigma^2 , epsilon) );     

% descent step
tau = step(lambda,W,rho); 

% Descent curves data
E = zeros(1,Nit_MAX); ConvRule = zeros(1,Nit_MAX);

% Initialisation
S = S0(1:rho,:);
Sm1 = S0(1:rho,:);
Sy = S0(1:rho,:);
t = 1;

% Projected Gradient Descent
n = 0;
r = 2;

tic
while ((n<2 || ConvE(E,tau,n) > 1e-3) && n<Nit_MAX )
    
    % Gradient descent
    S = Sy - tau*Gradient(Sy,lambda, W, rho);

    % Projection
    for samp = 1:Ns
        nDiff = norm(Yu(:,I(samp)) - S(:,I(samp)));
        if nDiff > sqrt(rho)*sigma*r
            u = (S(:,I(samp))-Yu(:,I(samp)))/nDiff;
            S(:,I(samp)) = Yu(:,I(samp))+ sqrt(rho)*sigma*r*u;
        end 
    end

    % t update
    tp1 = (1+sqrt(1+4*t^2))/2;


    % Sy update
    Sy = S + ((t-1)/tp1)*(S-Sm1);

    % Step update
    Sm1 = S;
    t = tp1;

    % Cost update
    n = n+1;
    % keyboard
    E(n) = cost(S,lambda,W,rho);
    if (n>1)
        ConvRule(n-1) = ConvE(E,tau,n);
    end
    
    % Display It #
    if(verb)
        if n<3
            fprintf('E is %.5f and n = %d\n',E(n),n);
        else    
            fprintf('E is %.5f and n = %d and diff = %.5f (goal is 1e-3)\n',E(n),n,ConvE(E,tau,n));
        end
    end 
    
end

t =toc;

if verb
    if (n < Nit_MAX)
        fprintf('nit fin %d \nElapsed time is %8f s\n',n,t)
    else
        fprintf('nit max reached %d \nElapsed time is %8f s\n',n,t)
    end
end

Xhat = reshape( transpose(H*S + repmat(Ym(:,1), [1 Np])), [nr nc Nb]);
                

end


function array_out=subarray(array_in, varargin)
% Author: Peter H. Mao, Caltech
%
% simple program to access an element of an array when the array is
% the output of another function
%
% usage array_out=subarray(array_in,dim1_indices,dim2_indices,....);
% to use the ':' symbol or the 'end' keyword, pass the range as a string
%
% intent: I call this function when I have a function that returns an array
% and I want to access some subset or element of that array.  this saves me from
% having a temporary variable in my workspace.

  dimlength = size(array_in);
  NonSingletonDimensions = find(dimlength > 1);
  n_dimensions = length(NonSingletonDimensions);
  if n_dimensions < length(varargin)
    warning('you may have too many index ranges specified');
  end

  kk = NonSingletonDimensions(1);
  for jj=1:length(varargin)
    if ischar(varargin{jj}) & ~strcmp(varargin{jj},':')
      range = regexprep(varargin{jj},'end', 'dimlength(kk)');
      varargin{jj} = eval(range);
    end
    kk=kk+1;
  end
  
  array_out = array_in(varargin{:});
  return;
end


function [ lout, Sigma, Rho ] = EigenEstimate( l, Ns )
% Computes an estimate of the covariance eigenvalues given the sample
% covariance eigenvalues. The Stein estimator coupled with isotonic
% regression has been used here.
%
% Input:
%   l: vector containing the sample eigenvalues
%   Ns: Number of observations
% Output:
%   lout: vector containing the estimated covariance matrix eigenvalues

% Check l
if (~isvector(l))
    error('Error. \nInput must be a vector, not a matrix');
end
turn = false;
if (~iscolumn(l))
    l = transpose(l);
    turn = true;
end

M = length(l); % data dimension

%% Stein estimator (cf MESTRE, Xavier. Improved estimation of eigenvalues and eigenvectors of covariance matrices using their sample estimates. IEEE Transactions on Information Theory, 2008, vol. 54, no 11, p. 5113-5129.)

% Stein Est
ListT = [l zeros(size(l))];

for ind = 1:M
    ListT(ind,2) = ( 1+(1/M)* sum((l(ind) + l(~ismember(1:M,ind))) ./(l(ind) - l(~ismember(1:M,ind)) ))  );
end



%% Isotonic regression

% Procedure 1st step
CellT = mat2cell(1:size(ListT,1),1,ones(1,size(ListT,1)));

while (sum(ListT(:,2)<0)>0)  % there are <0 alphai
    ind = size(ListT,1);
    while (ListT(ind,2)>0)
        ind = ind-1;
    end
    
    if (ind == 2)
        ListT = [sum(ListT( (ind-1):ind , :),1); ListT((ind+1):end,:)];
        CellT = [  {[CellT{1} CellT{2}]} CellT(3:end) ];
    elseif (ind == size(ListT,1))
        ListT = [ListT(1:(ind-2),:); sum(ListT( (ind-1):ind , :),1)];
        CellT = [ CellT(1:(ind-2)) {[CellT{ind-1} CellT{ind}]} ];
    else
        ListT = [ListT(1:(ind-2),:); sum(ListT( (ind-1):ind , :),1); ListT((ind+1):end,:)];
        CellT = [ CellT(1:(ind-2)) {[CellT{ind-1} CellT{ind}]} CellT((ind+1):end) ];
    end
end

% Procedure 2nd step

ListT(:,3) = ListT(:,1)./ListT(:,2);

ind = size(ListT,1)-1;
while (ind>1)
    
    while (ListT(ind+1,3)<ListT(ind,3) && ind>1)
        ind = ind-1;
    end
    if (ind>1)        
        if (ind == (size(ListT,1)-1))
            ListT = [ListT(1:(ind-1),:); sum(ListT( ind:(ind+1) , :),1)]; ListT(ind,3) = ListT(ind,1)/ListT(ind,2);
            CellT = [ CellT(1:(ind-1)) {[CellT{ind} CellT{ind+1}]} ];
            ind = ind-1;
        else
            ListT = [ListT(1:(ind-1),:); sum(ListT( ind:(ind+1) , :),1); ListT((ind+2):end,:)];  ListT(ind,3) = ListT(ind,1)/ListT(ind,2);
            CellT = [ CellT(1:(ind-1)) {[CellT{ind} CellT{ind+1}]} CellT((ind+2):end) ];
        end
    end
end

% Output

lout = zeros(M,1);

for ind = 1:length(ListT)
    lout(CellT{ind})=ListT(ind,3)*ones(1,length(CellT{ind}));
end

% if (Ns<=M)
%     % Because an outlier can ocure when not enough samples
%     ind = CellT{end-1};
%     Sigma = sqrt(lout(ind(1)));
% else
    Sigma = sqrt(lout(end));
% end


eps = 0.1;
temp = (lout>(1+eps)*Sigma^2); 
Rho = find(temp(1:(end-1)) - temp(2:end) );


%% Output
if (turn)
    lout = transpose(lout);
end

end

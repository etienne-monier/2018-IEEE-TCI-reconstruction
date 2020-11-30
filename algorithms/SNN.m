function [ Xhat, varargout ] = SNN( Y, I, varargin )
%Reconstructs a full multiband image using the Smoothed Nuclear Norm (SNN). 
%algorithm. 
%   Input:
%       Y: the acquired data of size (n,m,l) where n and m are the spatial
%        dimensions and l is the band number. Y should be filled with 0 at
%        non-sampled pixels,
%       I: vector of acquired pixel indexes (see sub2ind function).
%
%   Optional input:
%       'parameters': a 2-length positive vector to define optimization
%        problem parameters,
%       'init': sets the initial point for FISTA. Should be a (l,n*m)
%        matrix,
%       'verbose': set the display option by choosing 'y' whereas 'n'
%        option displays nothing.
%
%   Output:
%       Xhat: the reconstructed image under a (n,m,l) matrix form.
%
%   Optional output:
%       [Xhat, lmabda, mu]: gives back the OP parameters.    

%% Options

if(nargout~=1 && nargout ~=3)
    error('One or three outputs allowed.');
end

% default entries
defaultInitPoint = randn(size(Y,3),size(Y,1)*size(Y,2));
defaultverbose = 'y';
expectedShapes = {'y','n'};

% Parser parameters
p = inputParser;
addRequired(p,'Y',@(x) assert(isnumeric(x),'Value must be numeric.'));
addRequired(p,'I',@(x) assert(isnumeric(x) && isvector(x),'Value must be a numeric vector.'));

addParameter(p,'parameters',[], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2 && sum(x>=0)==2,'Value must be a 2-length positive and numeric vector.'));
addParameter(p,'init',defaultInitPoint, @(x) isnumeric(x)); % size(x) == size(Y) );
addParameter(p,'verbose',defaultverbose,  @(x) assert(any(validatestring(x,expectedShapes)),'Input must be y or n') );

% Parsing
parse(p,Y,I,varargin{:});

% post
verb = (p.Results.verbose == 'y');

paramSearch = isempty(p.Results.parameters);

if (~paramSearch) lambda = p.Results.parameters(1); mu = p.Results.parameters(2); end

InitPoint   = p.Results.init;


%% Parameter choice

if (paramSearch)
    [lambda, mu, Xhat] = SNNParameterSearch(Y, I, InitPoint, verb);
else
    [ Xhat ] = Reconstruction( lambda, mu, Y, I, InitPoint, verb );
end

if(nargout == 3)
    varargout{1} = lambda;
    varargout{2} = mu;
end



end

function [ Xhat ] = Reconstruction( lambda, mu, Y, I, InitPoint, verb )
% Here, the solution of 
%       Xhat = argMin_{X\in R^{Nb\times Ns}} \frac{1}{2}||Y-X_I||_F^2 +
%       \frac{\lambda}{2}||XD||_F^2
%
% is computed.
% 
% Input : 

%% Size
Ns = length(I);
[nr,nc,Nb] = size(Y);
Np = nr*nc;

Y = transpose(reshape(Y,[Np Nb]));

%% Sparse operators

temp = mat2cell( repmat( subarray( -eye(nr) + diag(1*ones(1,nr-1),-1),1:nr,1:(nr-1)) ,[1 nc] ) ,nr,repmat(nr-1,[1 nc]));
NablaX = spdiags([ones(nr*(nc-1),1), -ones(nr*(nc-1),1)],[-nr 0],nr*nc,nr*(nc-1));
NablaY = sparse(blkdiag( temp{:} ));
Nabla = [NablaY NablaX];
Delta = -Nabla*Nabla';


%% Parameters
% Mask
A = zeros(Nb,nr*nc); A(:,I) = 1;

% Stopping conditions
ConvE = @(E,step,n) abs(E(n)-E(n-1))/(E(n-1)*step);     % The chosen one

% Regularisation operators definition
step = @(lambda) 1/(1 + 8*lambda);
Gradient = @(lambda, Zn, Y) (Zn-Y).*A - lambda*Zn*Delta;
cost = @(Xn, Y, lambda, mu,s) (1/2)*norm(A.*(Xn-Y), 'fro')^2  + (lambda/2) *  norm( Xn*Nabla , 'fro')^2 + mu * s;

sThresholding = @(x,t) sign(x).*(abs(x)-t).*(abs(x)>t);

% Parameters
Nit_MAX = 2000;
lim = 1e-3;

%% Optimal X search

% Init matrix
X0 = InitPoint(1:Nb,1:nr*nc);

% descent step
tau = step(lambda);

% Descent curves data
E = zeros(1,Nit_MAX);

% Initialisation
Xn = X0;
Zn = X0;
tn = 1;

% iterations
n = 0;

nMsg = 0;
while (n<Nit_MAX &&(n<3 || ConvE(E,tau,n) > lim))
    
    % Display It #
    if(verb)
        fprintf(repmat('\b',1,nMsg));
        if n<3
            msg = ['n = ' num2str(n)];
        else    
            msg = ['n = ' num2str(n) ' and diff = ' num2str(ConvE(E,tau,n)) ' (goal is ' num2str(lim) ')' ];
        end
        fprintf(msg);
        nMsg=numel(msg);
    end
    
    % Gradient update
    Yn = Zn - tau*Gradient(lambda, Zn, Y);
    
    % Projection
    [U,S,V] = svd(Yn,'econ');
    Sbar = diag( sThresholding(diag(S),tau*mu));
    Xnp1 = U*Sbar*V'; s = sum(diag(Sbar));
    
    % mean parameter
    tnp1 = (1/2)*(1+sqrt(4*tn^2+1));
    lambdan = (tn-1)/(tnp1);

    Zn = Xn + lambdan*(Xnp1-Xn);

    % update
    tn=tnp1;
    Xn = Xnp1;
    n = n+1;

    % Error
    E(n) = cost(Xn, Y, lambda, mu,s);
    
end

% Display
if(verb)
    fprintf(repmat('\b',1,nMsg));
    fprintf(['nit fin ' num2str(n) '\n'])
end

Xhat = reshape(Xn',[nr nc Nb]);
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




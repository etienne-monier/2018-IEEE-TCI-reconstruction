function [lambda, mu, Xhat] = SNNParameterSearch(Y, I, InitPoint, verb)

%% Size and data
[nr, nc, Nb] = size(Y);
Np = nr*nc;
Ns = length(I);

%% Parameters

% Max number of it to find opt. lambda or mu or scale
Nmax = 1000;    
% Sets the optimal scale between upper and lower bound of lambda / mu
lim = 0.05;

% Limits for parameters
logLambdaMin(1) = -2; logLambdaMax(1) = 2;
logMuMin(1) = -2; logMuMax(1) = 2;

% Scale limits
scaleMin = 0.01; scaleMax = 1;

%% Estimation of sigma

Ytemp = transpose(reshape(Y,[Np Nb]));
[ sigma ] = madNoiseEst( Ytemp(:,I) );

%% Sparse operators

temp = mat2cell( repmat( subarray( -eye(nr) + diag(1*ones(1,nr-1),-1),1:nr,1:(nr-1)) ,[1 nc] ) ,nr,repmat(nr-1,[1 nc]));
NablaX = spdiags([ones(nr*(nc-1),1), -ones(nr*(nc-1),1)],[-nr 0],nr*nc,nr*(nc-1));
NablaY = sparse(blkdiag( temp{:} ));
Nabla = [NablaY NablaX];
Delta = -Nabla*Nabla';

%% First search : DataFidelity + lambda sobolev

if (verb) fprintf('\n#--------------------------------------#\nFirst Search : Lambda with spatial regularization\n'), end

ind_l = 1;
logLambda(1) = 0.5*(logLambdaMax(1) + logLambdaMin(1));

while ( 10^logLambdaMax(ind_l)/(10^logLambdaMin(ind_l)) - 1 > lim && ind_l < Nmax)
    
    % Compute Xhat for central Lambda
    [ Xhat ] = SearchOptimalLambda( 10^logLambda(ind_l), Ytemp, I, Ns, Nb, nr, nc, InitPoint, Delta, Nabla, verb );
    
    % Looks if the function ||Y-X_I||_F^2-NsNb*sigmaHat is positive or
    % negative at central point
    if ( norm(Ytemp(:,I) - Xhat(:,I)  ,'fro')^2-Ns*Nb*sigma^2 > 0 )
        % Keep left half of the segment
        logLambdaMin(ind_l+1) = logLambdaMin(ind_l);
        logLambdaMax(ind_l+1) = logLambda(ind_l);
    else
        % Keep right half of the segment
        logLambdaMin(ind_l+1) = logLambda(ind_l);
        logLambdaMax(ind_l+1) = logLambdaMax(ind_l);
    end
    
    % Update center point
    logLambda(ind_l+1) = 0.5*(logLambdaMax(ind_l+1) + logLambdaMin(ind_l+1));
    
    ind_l = ind_l + 1;
end

if (verb) disp(['Chosen Lambda is ' num2str(10^logLambda(end))]), end

lambda0 = 10^logLambda(end);


%% Second search : DataFidelity + mu nuclear norm

if (verb) fprintf('\n#--------------------------------------#\nSecond Search : Mu with spectral regularization\n'), end

% Init
ind_m = 1;
logMu(1) = 0.5*(logMuMax(1) + logMuMin(1));

while ( 10^logMuMax(ind_m)/(10^logMuMin(ind_m)) - 1 > lim && ind_m < Nmax)
    
    % Compute Xhat for central Lambda
    [ Xhat ] = SearchOptimalMu( 10^logMu(ind_m), Ytemp, I, Ns, Nb, nr, nc, InitPoint, Delta, Nabla, verb );
    
    % Looks if the function ||Y-X_I||_F^2-NsNb*sigmaHat is positive or
    % negative at central point
    if ( norm(Ytemp(:,I) - Xhat(:,I)  ,'fro')^2-Ns*Nb*sigma^2 > 0 )
        % Keep left half of the segment
        logMuMin(ind_m+1) = logMuMin(ind_m);
        logMuMax(ind_m+1) = logMu(ind_m);
    else
        % Keep right half of the segment
        logMuMin(ind_m+1) = logMu(ind_m);
        logMuMax(ind_m+1) = logMuMax(ind_m);
    end
    
    % Update center point
    logMu(ind_m+1) = 0.5*(logMuMax(ind_m+1) + logMuMin(ind_m+1));
    
    ind_m = ind_m + 1;
end

if (verb) disp(['Chosen Mu is ' num2str(10^logMu(end))]), end

mu0 = 10^logMu(end);


%% Third search : DataFidelity + lambda spatial reg. + mu nuclear norm

if (verb) fprintf('\n#--------------------------------------#\nThird Search : Scale with full regularisation\n'), end

% Init
ind_s = 1;
scale(1) = 0.5*(scaleMax(1) + scaleMin(1));

while ( 10^scaleMax(ind_s)/(10^scaleMin(ind_s)) - 1 > lim && ind_s < Nmax)
    
    % Compute Xhat for central Lambda
    [ Xhat ] = SearchOptimalScale( scale(ind_s)*10^logLambda(end), scale(ind_s)*10^logMu(end), Ytemp, I, Ns, Nb, nr, nc, InitPoint, Delta, Nabla,verb );
    
    % Looks if the function ||Y-X_I||_F^2-NsNb*sigmaHat is positive or
    % negative at central point
    if ( norm(Ytemp(:,I) - Xhat(:,I)  ,'fro')^2-Ns*Nb*sigma^2 > 0 )
        % Keep left half of the segment
        scaleMin(ind_s+1) = scaleMin(ind_s);
        scaleMax(ind_s+1) = scale(ind_s);
    else
        % Keep right half of the segment
        scaleMin(ind_s+1) = scale(ind_s);
        scaleMax(ind_s+1) = scaleMax(ind_s);
    end
    
    % Update center point
    scale(ind_s+1) = 0.5*(scaleMax(ind_s+1) + scaleMin(ind_s+1));
    
    ind_s = ind_s + 1;
end

scale0 = scale(end);

if (verb) 
    disp(['Chosen scale is ' num2str(scale(end))])
    disp('')
    fprintf('Final choices of parameters are lambda = %d and mu = %d\n\n',scale(end)*10^logLambda(end),scale(end)*10^logMu(end))
end

%% End


lambda = scale0*lambda0;
mu = scale0*mu0;

Xhat = reshape(Xhat',[nr nc Nb]);

end




function [ Xhat ] = SearchOptimalLambda( lambda, Y, I, Ns, Nb, nr, nc, InitPoint, Delta, Nabla, verb )
% Here, the solution of 
%       Xhat = argMin_{X\in R^{Nb\times Ns}} \frac{1}{2}||Y-X_I||_F^2 +
%       \frac{\lambda}{2}||XD||_F^2
%
% is computed.
% 
% Input : 

% Mask
A = zeros(Nb,nr*nc); A(:,I) = 1;

% Stopping conditions
ConvE = @(E,step,n) abs(E(n)-E(n-1))/(E(n-1)*step);

% Regularisation operators definition
step = @(lambda) 1/(8*lambda + 1);
Gradient = @(lambda, Zn, Y) (Zn-Y).*A - lambda*Zn*Delta;
cost = @(Xn, Y, I, lambda) (1/2)*norm(Xn(:,I)-Y(:,I), 'fro')^2  + (lambda/2) *  norm( Xn*Nabla , 'fro')^2;

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
    Xnp1 = Zn - tau*Gradient(lambda,Zn,Y);

    % mean parameter
    tnp1 = (1/2)*(1+sqrt(4*tn^2+1));
    lambdan = 1+(tn-1)/(tnp1);

    Zn = Xn + lambdan*(Xnp1-Xn);

    % update
    tn=tnp1;
    Xn = Xnp1;
    n = n+1;

    % Error
    E(n) = cost(Xn,Y,I,lambda);

end

% Display
if (verb)
    fprintf(repmat('\b',1,nMsg));
    fprintf(['nit fin ' num2str(n) '\n'])
end

Xhat = Xn;


end


function [ Xhat ] = SearchOptimalMu( mu, Y, I, Ns, Nb, nr, nc, InitPoint, Delta, Nabla, verb )
% Here, the solution of 
%       Xhat = argMin_{X\in R^{Nb\times Ns}} \frac{1}{2}||Y-X_I||_F^2 +
%       \frac{\lambda}{2}||XD||_F^2
%
% is computed.
% 
% Input : 

% Mask
A = zeros(Nb,nr*nc); A(:,I) = 1;

% Stopping conditions
ConvE = @(E,step,n) abs(E(n)-E(n-1))/(E(n-1)*step);     % The chosen one

% Regularisation operators definition
Gradient = @(Zn, Y) (Zn-Y).*A;
cost = @(Xn, Y, mu,s) (1/2)*norm(A.*(Xn-Y), 'fro')^2  + mu * s;

sThresholding = @(x,t) sign(x).*(abs(x)-t).*(abs(x)>t);

% Parameters
Nit_MAX = 2000;
lim = 1e-3;

%% Optimal X search

% Init matrix
X0 = InitPoint(1:Nb,1:nr*nc);

% descent step
tau = 1;

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
    Yn = Zn - tau*Gradient(Zn,Y);
    
    % Projection
    [U,S,V] = svd(Yn,'econ');
    Sbar = diag( sThresholding(diag(S),tau*mu) );
    Xnp1 = U*Sbar*V';
    
    % mean parameter
    tnp1 = (1/2)*(1+sqrt(4*tn^2+1));
    lambdan = 1+(tn-1)/(tnp1);

    Zn = Xn + lambdan*(Xnp1-Xn);

    % update
    tn=tnp1;
    Xn = Xnp1;
    n = n+1;

    % Error
    E(n) = cost(Xn,Y,mu,sum(diag(Sbar)));

end

% Display
if (verb)
    fprintf(repmat('\b',1,nMsg));
    fprintf(['nit fin ' num2str(n) '\n'])
end

Xhat = Xn;

end


function [ Xhat ] = SearchOptimalScale( lambda, mu, Y, I, Ns, Nb, nr, nc, InitPoint, Delta, Nabla, verb )
% Here, the solution of 
%       Xhat = argMin_{X\in R^{Nb\times Ns}} \frac{1}{2}||Y-X_I||_F^2 +
%       \frac{\lambda}{2}||XD||_F^2
%
% is computed.
% 
% Input : 

% Mask
A = zeros(Nb,nr*nc); A(:,I) = 1;

% Stopping conditions
ConvE = @(E,step,n) abs(E(n)-E(n-1))/(E(n-1)*step);     % The chosen one

% Regularisation operators definition
step = @(lambda) 1/(8*lambda + 1);
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
    Sbar = diag( sThresholding(diag(S),tau*mu) );
    Xnp1 = U*Sbar*V';
    
    % mean parameter
    tnp1 = (1/2)*(1+sqrt(4*tn^2+1));
    lambdan = 1+(tn-1)/(tnp1);

    Zn = Xn + lambdan*(Xnp1-Xn);

    % update
    tn=tnp1;
    Xn = Xnp1;
    n = n+1;

    % Error
    E(n) = cost(Xn, Y, lambda, mu,sum(diag(Sbar)));
end

% Display
if(verb)
    fprintf(repmat('\b',1,nMsg));
    fprintf(['nit fin ' num2str(n) '\n'])
end

Xhat = Xn;
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
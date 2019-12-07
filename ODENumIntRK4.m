function xOut = ODENumIntRK4(rhs,time,x0,varargin)

% Implement the Euler method to numerically integrate ODEs.
% 
% Usage:
%   xOut = ODENumIntEuler(rhs,time,x0)
%
% Inputs:
%   rhs     Function handle for the "right hand side" dynamic equations.
%           Must be defined like this:  rhs = @(t,x) YOUR_RHS_FUNCTION(...)
%
%   time    Time vector                   (1,nt)  
%           The time step between each (i) and (i+1) set of elements in the
%           time vector is used as the time step for the integration.
%
%   x0      Initial state vector          (ns,1) or (1,ns)
%
% varargin  (Variable argument input) Takes additional parameter inputs 
%           and places them in a cell array in the order they were input.
%           to be used if the rhs function requires additional parameters
%           to solve.
%           If no inputs beyond the required 3 are input: varargin = 0x0
%
% Outputs:
%   xOut    Integrated state over time    (nt,ns)
%
params = varargin;  % Additional inputs for rhs
nt = length(time);  % number of time points
ns = length(x0);    % number of states

% here we initialize the xOut matrix. 
% we will fill in its contents below
xOut = zeros(nt,ns); 

% first, record the initial state in the first row of xOut
xOut(1,:) = x0;

% numerically integrate using the Euler method
for j=1:nt-1
  
  xj = xOut(j,:)';        
  % NOTE: after transposing with the ' operator, xj is a column vector
  h = time(j+1)-time(j);
  
  % slope at time point j, or you could say the slope is dx/dt
  K1 = rhs(time(j),xj, params{1:end});  
  K2 = rhs(time(j+1)/2, xj + K1 * h/2,params{1:end});
  K3 = rhs(time(j+1)/2, xj + K2 * h/2,params{1:end});
  K4 = rhs(time(j+1) , xj + K3 * h,params{1:end});
  
  K = (K1/6 + K2/3 + K3/3 + K4/6);  
  
  % The "rhs" method should work the same whether xj is a column or row vector.
  % However, note that we expect the output of rhs to always be a column vector.
      
  % estimated state value at time point j+1 
  xNext = xj + h*K;   
  % both xj and K are column vectors, so "xNext" is also a column vector

  % record it in xOut
  xOut(j+1,:) = xNext; 
  % This puts "xNext" (a column vector) into the j+1 row of "xOut".
  % Matlab automatically reshapes it from a column into a row here.
  if(norm(xOut(j,1:3))-3390 < 5)
      xOut = xOut(1:j,:);
      break
  end
end



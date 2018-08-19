function [L,M,X,S] = logfactorial(N,varargin)

% logfactorial implementation from MATLAB file exchange
% https://www.mathworks.com/matlabcentral/fileexchange/14920-log-of-factorial-of-large-numbers
% Copyright (c) 2007, Yvan Lengwiler 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% * Redistributions of source code must retain the above copyright 
% notice, this list of conditions and the following disclaimer. 
% * Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in 
% the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%
% LOGFACTORIAL computes log10 of the factorial.
%
% USAGE:
%   logfactorial(N)
%   logfactorial(N,method)
%   logfactorial(N,method,fmt)
%   logfactorial(N,'sum',fmt,lengthOfSequence)
%   [M,X,L,S] = logfactorial(N,...);
%
% INPUT:
%   N is a positive integer or an array of positive integers.
%
%   method is either 'gamma' (the default) or 'sum'. It selects the method
%   that is used for computation.
%
%   fmt is the format in which the mantissa in S is printed. This argument
%   is optional; default format is '%g'.
%
%   lengthOfSequence is the length of a vector that MATLAB is able to
%   assign. One can set this to Inf, but this entails the risk of running
%   out of memory. The default is 1E7. If MATLAB runs out of memory, choose
%   a smaller number. This option is only relevant when using the 'sum'
%   method.
%
% OUTPUT:
%   L is an array of doubles containing log10 of the factorials.
%   M is an array of doubles containing the mantissas of the results.
%   X is an array of integers containing the exponents of the results.
%   S is an array of cellstrings that represents the results as strings in
%   the form 'Me+X'.
%
%   The factorial of N is M .* 10.^X, or equivalently 10.^L. Consequently,
%   10.^logfactorial(N) produces the same result as MATLAB's standard
%   function factorial(N).
%
% EXAMPLE 1
%   logfactorial(1E6)
%   produces 5.5657e+006, which is log10(1E6!).
%
%   [L,M,X,S] = logfactorial(1E6)
%   produces
%   L = 5.5657e+006         % log10(N!)
%   M = 8.2639              % mantissa
%   X = 5565708             % exponent
%   S = '8.26393e+5565708'  % N! as a string
%
% EXAMPLE 2
%   [L,M,X,S] = logfactorial(21:25)
%   returns
%   L =  19.7083   21.0508   22.4125   23.7927   25.1906
%   M =   5.1091    1.1240    2.5852    6.2045    1.5511
%   X =  19        21        22        23        25
%   S =  '5.10909e+19' '1.124e+21' '2.5852e+22' '6.20448e+23' '1.55112e+25'
%
% REMARK: The 'gamma' method computes the log (wrt base 10) of the
% factorial by means of the gamma function. The 'sum' method computes the
% log of the factorial as sum(log10(i)), where i = 1 ... n. The
% exponential of the factorial is just the integer part of log10 of the
% factorial; the mantissa is 10 to the fractional part.
% The 'gamma' method is generally faster and more accurate, and should
% therefore always be used. The 'sum' method is provided here only for
% comparison, and because the algorithm is quite transparent.
%
% Author: Yvan Lengwiler (yvan.lengwiler@unibas.ch)
% Date: May 9, 2007
%
% Acknowledgement: Urs Schwarz and, in particular, John D'Errico (reviewers
% at MathCentral) for helpful comments.
%
% See also: FACTORIAL, LARGEFACTORIAL (MathCentral file id 14816)

% Compare the two methods with factorial.m for N = 1 ... 170 as follows:
%   dgamma = @(x) log10(factorial(x)) - logfactorial(x,'gamma');
%   dsum = @(x) log10(factorial(x)) - logfactorial(x,'sum');
%   hold all
%   plot(dgamma(1:170))
%   plot(dsum(1:170))
%   legend('gamma method','sum method')
% One can see the differences in precision already at these small values of
% N. The loss of precision gets worse for larger arguments.

    % --- check validity of arguments
    
    Nsort = N(:);   % This is not sorted yet, but will be sorted later. The
                    % same variable name is used to save memory.

    if any(fix(Nsort) ~= Nsort) || any(Nsort < 0) || ...
            ~isa(Nsort,'double') || ~isreal(Nsort)
        error('N must be a matrix of non-negative integers.')
    end
    
    if length(varargin) > 3
        error('Too many input arguments.')
    end
    
    % --- assign default arguments if necessary
    % (this trick is taken from Nabeel Azar, see
    % http://www.ee.columbia.edu/~marios/matlab/Programming%20Patterns%20Think%20Globally,%20Act%20Locally.pdf)
    
    defaultValues                 = {'gamma','%g',1E7};
    nonemptyIdx                   = ~cellfun('isempty',varargin);
    defaultValues(nonemptyIdx)    = varargin(nonemptyIdx);
    [method fmt lengthOfSequence] = deal(defaultValues{:});
    
    % --- perform computation
    
    switch lower(method)
        
        case 'gamma'	% ... using gammaln

            L = gammaln(N+1)/log(10);   % John D'Errico one line solution
                                        % of the problem
            
        case 'sum'      % ... using direct summation

            % We compute the smallest factorial first, ...
            [Nsort,map] = sort(Nsort);
            Nsort = [0;Nsort];
            % ... and then compute only the increments.
            L = arrayfun(@sumlogIncrement, 1:numel(N));
            L = cumsum(L);  % aggregate increments
            L(map) = L;     % put results back into original ordering
                            % ("un-sort")
            L = reshape(L,size(N));     % bring output into correct shape

        otherwise

            error(['Method ''', method, ''' is unknown.'])

    end

    if nargout > 1
        % 10^L is the factorial we are after, so ...
        X = fix(L);      % ... X is the exponential of the factorial ...
        M = 10.^(L-X);   % ... and M is the mantissa.
    end
    
    if nargout > 3
        S = arrayfun(@makeString, 1:numel(N));
        S = reshape(S,size(N));     % bring output into correct shape
    end

    % --- nested function -------------------------------------------------
    function s = sumlogIncrement(j)
    % sumlogIncrement(j) computes sum(log10(i)) for i = Nsort(j)+1 ...
    % Nsort(j+1); this is used for the 'sum' method only.
        lengthInterval = Nsort(j+1)-Nsort(j);
        a = Nsort(j);
        if lengthInterval <= lengthOfSequence
            s = sum(log10(a+(1:lengthInterval)));
        else
            % MATLAB sometimes runs out of memory when allocating long
            % vectors. It is assumed that MATLAB can always allocate a
            % vector of size 'lengthOfSequence' (chosen by user, or set to
            % 1E7 by default).
            seq = 1:lengthOfSequence;   % "short" sequence
            rounds = floor(lengthInterval/lengthOfSequence);  % # of rounds
            partialSum = zeros(1,rounds);
            for r = 1:rounds
                partialSum(r) = sum(log10(seq+(r-1)*lengthOfSequence+a));
            end
            s = sum(partialSum) + ...
                sum(log10(((rounds*lengthOfSequence+1):lengthInterval)+a));
        end
    end

    function str = makeString(j)
    % makeString(j) creates a string in the form 'Me+X'.
        str = cellstr(strcat(num2str(M(j),fmt), 'e+', int2str(X(j))));
    end

end
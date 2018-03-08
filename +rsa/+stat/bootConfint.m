function [ci_lo, ci_up, p] = bootConfint(x, B, alternative, user_opts)

% bootConfint
%   Compute confidence intervals based on bootstrap distribution using
%   various methods and use them to derive p-values for hypothesis testing.
% 
%   [ci_lo, ci_up, p] = bootConfint(x, B, alternative, user_opts)
%
%   Currently two methods are available, "percentile" and "bca". P-values
%   can be computed for both one-tailed and two-tailed testing.
%
%   INPUT
%   -----
%   x           : vector, original sample of length #subjects
%   B           : vector, bootstrap distribution of the statistic of 
%                 interest of length #bootstrap_iterations
%   alternative : string, either 'two-tailed' or 'greater'
%   user_opts   : struct, user options used by RSA Toolbox. Options are ...
%                 user_opts.bootCImethod:
%                 - "percentile":  bootstrap percentile method (default)
%                 - "bca": bias-corrected and accelerated CI
%                 user_opts.bootCIalpha: type 1 error probability to use
%                 for CI definition, default is .05.
%
%   OUTPUT
%   ------
%   ci_lo, ci_up: lower and upper bounds of CI, respectively
%   p           : p value of hypothesis test
%
%   Literature   
%   ----------
%   "An Introduction to the Bootstrap" (Efron & Tibshirani, 1993)
%
%   2018-01-27 michael.bannert@tuebingen.mpg.de

user_opts = rsa.util.setIfUnset(user_opts, 'bootCImethod', 'percentile');
user_opts = rsa.util.setIfUnset(user_opts, 'bootCIalpha', .05);

method = user_opts.bootCImethod;
ci_alpha = user_opts.bootCIalpha;

n_B = length(B);
n_samples = length(x);
B_sorted = sort(B);

switch lower(method)
    case 'bca'        
        % Bias-correction, eq. 14.14
        z_hat0 = norminv(sum(B < mean(x)) / n_B);

        % Compute acceleration using jackknife
        theta_i = nan(n_samples, 1);
        for n = 1:n_samples
            x_i = x(setdiff(1:n_samples, n));
            theta_i(n) = mean(x_i);
        end
        theta_pt = mean(theta_i);

        % Eq. 14.15
        a_hat = sum((theta_pt * ones(n_samples, 1) - theta_i).^3) / ...
            6 * (sum((theta_pt * ones(n_samples, 1) - theta_i).^2)) ^(3/2);

        z_alpha = norminv(ci_alpha / 2);
        z_one_minus_alpha = norminv(1 - ci_alpha / 2);

        % Eq. 14.10
        alpha1 = normcdf( ...
            z_hat0 + (z_hat0 + z_alpha) / ...
            (1 - a_hat * (z_hat0 + z_alpha)) ...
            );

        alpha2 = normcdf( ...
            z_hat0 + (z_hat0 + z_one_minus_alpha) / ...
            (1 - a_hat * (z_hat0 + z_one_minus_alpha)) ...
            );

        % Find indices of bootstrap distribution corresponding to desired 
        % alpha
        alpha1_n = max(floor(alpha1 * n_B), 1); % At least 1
        alpha2_n = ceil(alpha2 * n_B);

        % Eq. 14.9
        ci_lo = B_sorted(alpha1_n);
        ci_up = B_sorted(alpha2_n);
        
        % Calculate p-value by solving eq. 14.10 for alpha. P-value should 
        % be one-tailed for the null hypothesis that the true value is <= 0.
        % Find the z_alpha for which the threshold 0 equals the lower CI 
        % bound.

        % 1. Precompute inverse cumulative Gaussian of lower confidence 
        % bound that would just include 0. (Used several times in equation 
        % below.)
        
        % Index must be at least 1. Increase number of bootstrap iterations
        % for more precise p-values
        alpha0_n = find(B_sorted > 0, 1) - 1;
        alpha1_norminv = norminv(alpha0_n / n_B);        
   
        % 2. Find CI alpha for which the lower bound would just enclose 
        % 0. 
        % This is eq. 14.10 solved for z_alpha
        z_alpha0 = ...
            (alpha1_norminv * (1 - a_hat * z_hat0) + ...
            z_hat0 * (z_hat0 * a_hat - 2)) / ...
            (a_hat * (alpha1_norminv - z_hat0) - 1);
        
        % Super non-significant, i.e., no positive values in bootstrap
        % distribution
        if isempty(z_alpha0)
            z_alpha0 = -inf;
        end
        
        % Super significant, i.e., only positive values in bootstrap
        % distribution. P-value only limited by number of bootstrap
        % iterations. Increasing number of iterations would possibly help.
        if isnan(z_alpha0)
            z_alpha0 = norminv(1 - 1/n_B);
        end
            
        switch lower(alternative)
            case 'greater'
                % Convert z_alpha to p-value
                p = 1 - normcdf(z_alpha0);
                
            case 'two-tailed'
                % For which z_alpha does the upper CI bound equal 0? See
                % one-tailed case for explanation
                
                % Index must be at least n_B - 1. Increase number of
                % bootstraps for more precise p-values
                alpha0_n_up = min(n_B-1, find(B_sorted < 0, 1, 'last') + 1);
                if isempty(alpha0_n_up)
                    alpha0_n_up = n_B-1;
                end
                alpha1_norminv_up = norminv(alpha0_n_up / n_B);               
                
                z_alpha0_up = ...
                    (alpha1_norminv_up * (1 - a_hat * z_hat0) + ...
                    z_hat0 * (z_hat0 * a_hat - 2)) / ...
                    (a_hat * (alpha1_norminv_up - z_hat0) - 1);

                % Compare the two one-tailed p-values and return the
                % smaller one multiplied by two.
                p = min(1 - normcdf(z_alpha0), normcdf(z_alpha0_up)) * 2;
                if isempty(p) || p > 1
                    keyboard;
                end
        end
                    
    case 'percentile'
        % Efron & Tibshirani (1993, p. 170)
        ci_lo = B_sorted(max(1, floor((ci_alpha / 2) * n_B)));
        ci_up = B_sorted(ceil((1 - ci_alpha / 2) * n_B));
        
        switch lower(alternative)
            case 'greater'
                % One-tailed p-value for null hypothesis that true value is
                % below or equal to 0.
                p = max(sum(B < 0) / n_B, 1/n_B);
                
            case 'two-tailed'
                % Two-tailed p-value for null hypothesis that true value
                % lies within confidence interval
                p = max((min(sum(B < 0), sum(B >= 0)) / n_B) * 2, 1/n_B);
        end
        
    otherwise
        error(['Valid method ''%s'' to compute bootstrap CIs\nValid ' ...
            'methods are: ''bca'', ''percentile'''], method);
end

assert(p <= 1 && p >= 0, '''p'' is not a probability.');

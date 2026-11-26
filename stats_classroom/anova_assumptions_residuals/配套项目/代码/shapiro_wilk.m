function [W, p] = shapiro_wilk(x)
%SHAPIRO_WILK 完整独立样本的Shapiro-Wilk正态性检验。
%   [W,P] = SHAPIRO_WILK(X) returns the Shapiro-Wilk statistic W and its
%   p-value P for the finite, non-constant real vector X. The implementation
%   accepts 3 <= numel(X) <= 5000 and does not remove missing values.
%
%   This is an uncensored MATLAB port of the SciPy 1.17.1 swilk routine,
%   which follows Royston's AS R94 algorithm. See the accompanying
%   LICENSE_shapiro_wilk.txt and the source copy in
%   https://github.com/scipy/scipy/blob/v1.17.1/scipy/stats/_ansari_swilk_statistics.pyx. Shapiro-Francia is
%   not used here.

    validateattributes(x, {'numeric'}, {'real', 'vector', 'finite'}, ...
        mfilename, 'x', 1);
    n = numel(x);
    if n < 3 || n > 5000
        error('shapiro_wilk:SampleSize', 'x must contain between 3 and 5000 values.');
    end
    if all(x(:) == x(1))
        error('shapiro_wilk:ConstantSample', 'x must not be constant.');
    end
    x = sort(double(x(:)));

    % Constants are the coefficients in SciPy swilk.pyx. They are supplied
    % in ascending-power order; fliplr makes them suitable for polyval.
    c1 = [0, 0.221157, -0.147981, -0.207119e1, 0.4434685e1, -0.2706056e1];
    c2 = [0, 0.42981e-1, -0.293762, -0.1752461e1, 0.5682633e1, -0.3582633e1];
    c3 = [0.5440, -0.39978, 0.25054e-1, -0.6714e-3];
    c4 = [0.13822e1, -0.77857, 0.62767e-1, -0.20322e-2];
    c5 = [-0.15861e1, -0.31082, -0.83751e-1, 0.38915e-2];
    c6 = [-0.4803, -0.82676e-1, 0.30302e-2];
    g = [-0.2273e1, 0.459];

    if n == 3
        a = 1 / sqrt(2);
    else
        n2 = floor(n / 2);
        an = double(n);
        an25 = an + 0.25;
        m = norminv(((1:n2) - 0.375) / an25);
        m = m(:);
        summ2 = 2 * sum(m.^2);
        ssumm2 = sqrt(summ2);
        rsn = 1 / sqrt(an);
        A1 = polyval(fliplr(c1), rsn) - m(1) / ssumm2;
        if n > 5
            A2 = -m(2) / ssumm2 + polyval(fliplr(c2), rsn);
            fac = sqrt((summ2 - 2*m(1)^2 - 2*m(2)^2) / ...
                (1 - 2*A1^2 - 2*A2^2));
            a = m / fac;
            a(2) = A2;
            i1 = 3;
        else
            fac = sqrt((summ2 - 2*m(1)^2) / (1 - 2*A1^2));
            a = m / fac;
            i1 = 2;
        end
        a(1) = A1;
        a(i1:end) = -a(i1:end);
    end

    % Pair the positive half-coefficients with high minus low order values.
    d = x(end:-1:end-numel(a)+1) - x(1:numel(a));
    S = sum(a .* d);
    Sxx = sum((x - mean(x)).^2);
    W = (S^2) / Sxx;
    W = min(1, max(0, W));

    if n == 3
        p = 1 - (6 / pi) * acos(sqrt(W));
        p = min(1, max(0, p));
        return;
    end

    w1 = max(realmin, 1 - W);
    y = log(w1);
    if n <= 11
        gamma = polyval(fliplr(g), n);
        if y >= gamma
            p = 1e-19;
            return;
        end
        y = -log(gamma - y);
        mu = polyval(fliplr(c3), n);
        sigma = exp(polyval(fliplr(c4), n));
    else
        ln = log(n);
        mu = polyval(fliplr(c5), ln);
        sigma = exp(polyval(fliplr(c6), ln));
    end
    p = 0.5 * erfc((y - mu) / (sigma * sqrt(2)));
    p = min(1, max(0, p));
end

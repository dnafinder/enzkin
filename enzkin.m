function enzkinout = enzkin(S, v)
%ENZKIN Estimate Michaelis-Menten parameters (Km, Vmax) with multiple regressions
%   ENZKINOUT = ENZKIN(S, V) estimates Km and Vmax from initial-rate
%   enzyme kinetics data using six different regression approaches and
%   generates diagnostic plots.
%
%   INPUT
%     S - Row vector of substrate concentrations [S].
%         Must be numeric, real, finite, nonnegative, increasing.
%
%     V - Row vector of initial velocities v corresponding to S.
%         Must be numeric, real, finite, nonnegative, increasing.
%
%   OUTPUT
%     ENZKINOUT - Structure with fields:
%         ENZKINOUT.KM   - 6 x 4 matrix with Km estimates:
%                          [value, standard error, lower CI, upper CI]
%         ENZKINOUT.VMAX - 6 x 4 matrix with Vmax estimates:
%                          [value, standard error, lower CI, upper CI]
%
%     Rows correspond to:
%       1) Lineweaver-Burk   (linear: 1/v vs 1/[S])
%       2) Hanes-Woolf       (linear: [S]/v vs [S])
%       3) Eadie-Hofstee     (linear: v vs v/[S])
%       4) Scatchard         (linear: v/[S] vs v)
%       5) Logarithmic non-linear Michaelis-Menten
%       6) Hyperbolic (standard) non-linear Michaelis-Menten
%
%   Besides returning ENZKINOUT (if requested), the function:
%     - Prints detailed regression summaries to the Command Window
%     - Displays tables of Km and Vmax estimates (with CI)
%     - Creates a figure with 6 diagnostic subplots:
%         * Michaelis-Menten hyperbolic plot
%         * Lineweaver-Burk plot
%         * Hanes-Woolf plot
%         * log(Michaelis-Menten) plot
%         * Eadie-Hofstee plot
%         * Scatchard plot
%
%   REQUIREMENTS
%     - Custom regression function MYREGR by G. Cardillo:
%         https://it.mathworks.com/matlabcentral/fileexchange/15473-myregression
%     - Curve Fitting Toolbox (FIT, FITTYPE) for non-linear regressions.
%
%   EXAMPLE
%     S = [0.5 1 2 4 8];
%     v = [0.12 0.22 0.39 0.55 0.65];
%     out = enzkin(S, v);
%
%   ------------------------------------------------------------------
%   Author and citation:
%   ------------------------------------------------------------------
%   Created by:  Giuseppe Cardillo
%   E-mail:      giuseppe.cardillo.75@gmail.com
%
%   To cite this file:
%   Cardillo G. (2010). enzkin.m – A tool to estimate Michaelis-Menten
%   kinetic parameters using multiple linear and non-linear regressions.
%
%   GitHub repository:
%   https://github.com/dnafinder/enzkin
%   ------------------------------------------------------------------

%% Input error handling
p = inputParser;
addRequired(p, 'S', @(x) validateattributes(x, {'numeric'}, ...
    {'vector','real','finite','nonnan','nonnegative','row','increasing'}));
addRequired(p, 'v', @(x) validateattributes(x, {'numeric'}, ...
    {'vector','real','finite','nonnan','nonnegative','row','increasing'}));
parse(p, S, v);
S = p.Results.S;
v = p.Results.v;
assert(length(S) == length(v), 'S and v must have the same length.');
clear p

%% Check for MYREGR dependency
assert(exist('myregr.m', 'file') ~= 0, ...
    ['You must download the myregr function from ' ...
     'https://it.mathworks.com/matlabcentral/fileexchange/15473-myregression ' ...
     'and add it to your MATLAB path.']);

%% Basic constants and containers
n  = length(S);           % number of data points
vc = tinv(0.95, n - 2);   % critical t-value for 95% confidence intervals
tr = repmat('-', 1, 100); % text divider line

% Labels for the plots (used by kingraph)
txtlbl = { ...
    'Michaelis & Menten non linear fit',      '[S]',   'v'; ...
    sprintf('Lineweaver-Burk\n(x=1/S; y=1/v)'), '1/[S]','1/V'; ...
    sprintf('Hanes-Woolf\n(x=[S]; y=[S]/v)'), '[S]',   '[S]/v'; ...
    sprintf('Logarithmic non linear fit'),    'log(S)','log(v)'; ...
    sprintf('Eadie-Hofstee\n(x=v/[S]; y=v)'), 'v/[S]','v'; ...
    sprintf('Scatchard\n(x=v; y=v/[S])'),     'v',     'v/[S]'}; 

% Matrices to store Km and Vmax estimates:
%   columns: [value, standard error, lower CI, upper CI]
%   rows   : 1..6 methods
KM   = NaN(6, 4);
VMAX = NaN(6, 4);

%% 1) Lineweaver-Burk linearization
disp(tr)
fprintf('Lineweaver-Burk linearization (x=1/S; y=1/v => Vmax=1/q; Km=m/q)...\n');
disp(tr)
[x, Idx] = sort(1 ./ S); 
y        = 1 ./ v(Idx);
[slope, intercept, stat] = myregr(x, y, 0);

VMAX(1, 1:2) = [1 / intercept.value, intercept.se / intercept.value^2];
VMAX(1, 3:4) = VMAX(1, 1) + [-1 1] .* vc .* VMAX(1, 2);

KM(1, 1:2) = [slope.value / intercept.value, ...
    realsqrt((intercept.value * slope.se)^2 + ...
             (slope.value * intercept.se)^2) / intercept.value^2];
KM(1, 3:4) = KM(1, 1) + [-1 1] .* vc .* KM(1, 2);

disppar(1)
clear x y slope intercept stat
disp(' ')

%% 2) Hanes-Woolf linearization
disp(tr)
fprintf('Hanes-Woolf linearization (x=[S]; y=[S]/v => Vmax=1/m; Km=q/m)...\n');
disp(tr)
[slope, intercept, stat] = myregr(S, S ./ v, 0);

VMAX(2, 1:2) = [1 / slope.value, slope.se / slope.value^2];
VMAX(2, 3:4) = VMAX(2, 1) + [-1 1] .* vc .* VMAX(2, 2);

KM(2, 1:2) = [intercept.value / slope.value, ...
    realsqrt((slope.value * intercept.se)^2 + ...
             (intercept.value * slope.se)^2) / slope.value^2];
KM(2, 3:4) = KM(2, 1) + [-1 1] .* vc .* KM(2, 2);

disppar(2)
clear slope intercept stat
disp(' ')

%% 3) Eadie-Hofstee linearization
disp(tr)
[x, Idx] = sort(v ./ S); 
y        = v(Idx);
fprintf('Eadie-Hofstee linearization (x=v/[S]; y=v => Vmax=q; Km=-m)...\n');
disp(tr)
[slope, intercept, stat] = myregr(x, y, 0);

VMAX(3, :)   = [intercept.value, intercept.se, intercept.lv, intercept.uv];
KM(3, 1:2)   = [-slope.value, slope.se];
KM(3, 3:4)   = KM(3, 1) + [-1 1] .* vc .* KM(3, 2);

disppar(3)
clear x y slope intercept stat
disp(' ')

%% 4) Scatchard linearization
disp(tr)
[x, Idx] = sort(v); 
y        = v(Idx) ./ S(Idx);
fprintf('Scatchard linearization (x=v; y=v/[S] => Vmax=-q/m; Km=-1/m)...\n');
disp(tr)
[slope, intercept, stat] = myregr(x, y, 0);

KM(4, 1:2) = [-1 / slope.value, slope.se / slope.value^2];
KM(4, 3:4) = KM(4, 1) + [-1 1] .* vc .* KM(4, 2);

VMAX(4, 1:2) = [-intercept.value / slope.value, ...
    realsqrt((-intercept.se / slope.value)^2 + ...
             (-intercept.value / slope.value^2 * slope.se)^2)];
VMAX(4, 3:4) = VMAX(4, 1) + [-1 1] .* vc .* VMAX(4, 2);

disppar(4)
clear x y slope intercept stat
disp(' ')

%% Non-linear log(Michaelis-Menten) fit
xfit  = S(:);
yfit  = v(:);
lyfit = log(v(:));

% Check finite values
ok = isfinite(log(xfit)) & isfinite(lyfit);
if ~all(ok)
    warning('enzkin:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data for log(Michaelis-Menten) fit.');
end

% Robust, simple initial guesses from raw data (NO dependence on KM/VMAX)
Vmax0 = max(yfit);
posS  = xfit(xfit > 0);
if isempty(posS)
    Km0 = max(xfit);
else
    Km0 = median(posS);
end

Km0   = max(Km0,   eps);
Vmax0 = max(Vmax0, eps);

% Fitting options: parameters Km and Vmax in linear scale
fo = fitoptions('Method', 'NonlinearLeastSquares', ...
                'Robust', 'On', ...
                'Lower', [0 0], ...
                'StartPoint', [Km0 Vmax0]);

% log(Michaelis-Menten) model: log(v) = log(Vmax) + log(S) - log(Km + S)
ft = fittype('log(Vmax)+log(x)-log(Km+x)', ...
    'dependent',   {'y'}, ...
    'independent', {'x'}, ...
    'coefficients', {'Km', 'Vmax'});

[cfl, goodness] = fit(xfit(ok), lyfit(ok), ft, fo);

disp(tr)
fprintf('Logarithmic nonlinear fit...\n');
disp(tr)
disp(cfl)
fprintf('\tR = \t%0.4f\n', realsqrt(goodness.rsquare));
disp(tr)
disp(' ')

p            = coeffvalues(cfl);
KM(5, 1)     = p(1);
VMAX(5, 1)   = p(2);
p            = confint(cfl);
KM(5, 3:4)   = p(1, :)';
VMAX(5, 3:4) = p(2, :)';
KM(5, 2)     = (KM(5, 4)   - KM(5, 1))   / vc;
VMAX(5, 2)   = (VMAX(5, 4) - VMAX(5, 1)) / vc;
clear p fo ft cfl goodness

%% Non-linear Michaelis-Menten hyperbolic fit
ok = isfinite(xfit) & isfinite(yfit);
if ~all(ok)
    warning('enzkin:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data for hyperbolic fit.');
end

% Simple initial guesses from raw data (independent of previous methods)
Vmax0 = max(yfit(ok));
posS  = xfit(ok);
posS  = posS(posS > 0);
if isempty(posS)
    Km0_h = max(xfit(ok));
else
    Km0_h = median(posS);
end

Km0_h   = max(Km0_h,   eps);
Vmax0_h = max(Vmax0,   eps);

% Fitting options for standard Michaelis-Menten
fo = fitoptions('Method', 'NonlinearLeastSquares', ...
                'Robust', 'On', ...
                'Lower', [0 0], ...
                'StartPoint', [Km0_h Vmax0_h], ...
                'DiffMaxChange', 9.9999999999999995475e-07);

% Standard Michaelis-Menten: v = (Vmax * S) / (Km + S)
ft = fittype('(Vmax*x)/(Km+x)', ...
    'dependent',   {'y'}, ...
    'independent', {'x'}, ...
    'coefficients', {'Km', 'Vmax'});

[cf, goodness] = fit(xfit(ok), yfit(ok), ft, fo);


disp(tr)
fprintf('Michaelis & Menten nonlinear fit...\n');
disp(tr)
disp(cf)
fprintf('\tR = \t%0.4f\n', sqrt(goodness.rsquare));
disp(tr)
disp(' ')

p            = coeffvalues(cf);
KM(6, 1)     = p(1);
VMAX(6, 1)   = p(2);
p            = confint(cf);
KM(6, 3:4)   = p(1, :)';
VMAX(6, 3:4) = p(2, :)';
KM(6, 2)     = (KM(6, 4)   - KM(6, 1))   / vc;
VMAX(6, 2)   = (VMAX(6, 4) - VMAX(6, 1)) / vc;
clear p fo ft cf goodness

%% Display summary tables for Km and Vmax
disp(array2table(KM, ...
    'VariableNames', {'KM', 'Std_Err', 'Lower_bound', 'Upper_bound'}, ...
    'RowNames', {'Lineweaver_Burk', 'Hanes_Woolf', 'Eadie_Hofstee', ...
                 'Scatchard', 'Log_fit', 'Hyperbolic_fit'}));

disp(array2table(VMAX, ...
    'VariableNames', {'VMAX', 'Std_Err', 'Lower_bound', 'Upper_bound'}, ...
    'RowNames', {'Lineweaver_Burk', 'Hanes_Woolf', 'Eadie_Hofstee', ...
                 'Scatchard', 'Log_fit', 'Hyperbolic_fit'}));

%% Output struct
if nargout
    enzkinout.KM   = KM;
    enzkinout.VMAX = VMAX;
end

%% Plot figure with 6 subplots
figure('Color', [1 1 1], 'OuterPosition', get(groot, 'ScreenSize'));

% Lineweaver-Burk plot 
m = KM(1, 1) / VMAX(1, 1);
q = 1 / VMAX(1, 1);
x = [-1 / KM(1, 1)  0  1 ./ S];
y = [0              q  1 ./ v];
stringa = { ...
    ['Vmax = ' num2str(VMAX(1, 1)) '    Km = ' num2str(KM(1, 1))]; ...
    'y = (Km/Vmax) * x + (1/Vmax)'; ...
    ['x=0; y=1/Vmax= ' num2str(y(2))]; ...
    ['y=0; x=-1/Km= ' num2str(x(1))]};
kingraph(2, m, q, stringa)

% Hanes-Woolf plot 
m = 1 / VMAX(2, 1);
q = KM(2, 1) / VMAX(2, 1);
x = [-KM(2, 1)  0  S];
y = [0          q  S ./ v];
stringa = { ...
    ['Vmax = ' num2str(VMAX(2, 1)) '    Km = ' num2str(KM(2, 1))]; ...
    'y = (1/Vmax) * x + (Km/Vmax)'; ...
    ['x=0; y=Km/Vmax=' num2str(y(2))]; ...
    ['y=0; x=-Km=' num2str(x(1))]};
kingraph(3, m, q, stringa)

% Eadie-Hofstee plot
m = -KM(3, 1);
q = VMAX(3, 1);
x = [VMAX(3, 1) / KM(3, 1)  0  v ./ S];
y = [0                      VMAX(3, 1)  v];
stringa = { ...
    ['Vmax = ' num2str(VMAX(3, 1)) '    Km = ' num2str(KM(3, 1))]; ...
    'y = -Km * x + Vmax'; ...
    ['x=0; y=Vmax=' num2str(y(2))]; ...
    ['y=0; x=Vmax/Km=' num2str(x(1))]};
kingraph(5, m, q, stringa)

% Scatchard plot
m = -1 / KM(4, 1);
q = VMAX(4, 1) / KM(4, 1);
x = [VMAX(4, 1)  0  v];
y = [0           VMAX(4, 1) / KM(4, 1)  v ./ S];
stringa = { ...
    ['Vmax = ' num2str(VMAX(4, 1)) '    Km = ' num2str(KM(4, 1))]; ...
    'y = -1/Km * x + Vmax/Km'; ...
    ['x=0; y=Vmax/Km=' num2str(y(2))]; ...
    ['y=0; x=Vmax=' num2str(x(1))]};
kingraph(6, m, q, stringa)

% Michaelis & Menten plot 
x = S; 
y = v;
[~, xtick, ytick] = kingraph(1, VMAX(6, 1), KM(6, 1));
H = text(5 * xtick, VMAX(6, 1) + 3 * ytick, ...
    ['Vmax = ' num2str(VMAX(6, 1))]); 
set(H, 'Color', 'm')
H = text(KM(6, 1) + 5 * xtick, VMAX(6, 1) / 2, ...
    ['Km = ' num2str(KM(6, 1))]); 
set(H, 'Color', 'm')

% log(Michaelis & Menten) plot 
[as, xtick, ytick] = kingraph(4, VMAX(5, 1), KM(5, 1));
H = text(as(1) + 5 * xtick, log(VMAX(5, 1)) + 2 * ytick, ...
    'log(v)=log(Vmax)+log(S)-log(Km+S)'); 
set(H, 'Color', 'm')
H = text(as(1) + 5 * xtick, log(VMAX(5, 1)) + 5 * ytick, ...
    ['log(Vmax) = ' num2str(log(VMAX(5, 1))) '; Vmax = ' num2str(VMAX(5, 1))]); 
set(H, 'Color', 'm')
H = text(log(KM(5, 1)) + 5 * xtick, log(VMAX(5, 1) / 2), ...
    ['log(Km) = ' num2str(log(KM(5, 1))) '; Km = ' num2str(KM(5, 1))]); 
set(H, 'Color', 'm')

%% Nested helper functions

    function disppar(I)
        %DISPPAR Display regression parameters for method I
        z = [struct2array_local(slope), ...
                (1 - tcdf(abs(slope.value / slope.se), n - 2)) * 2; ...
             struct2array_local(intercept), ...
                (1 - tcdf(abs(intercept.value / intercept.se), n - 2)) * 2; ...
             stat.r(1:4), ...
                (1 - tcdf(abs(stat.r(1) / stat.r(2)), n - 2)) * 2; ...
             VMAX(I, :), ...
                (1 - tcdf(abs(VMAX(I, 1) / VMAX(I, 2)), n - 2)) * 2; ...
             KM(I, :), ...
                (1 - tcdf(abs(KM(I, 1) / KM(I, 2)), n - 2)) * 2; ...
            ];
        disp(array2table(z, ...
            'RowNames', {'Slope', 'Intercept', 'R', 'Km', 'Vmax'}, ...
            'VariableNames', {'Value', 'SE', 'Lower_bound', 'Upper_bound', 'p_value'}));
        disp(tr)
    end

    function a = struct2array_local(s)
        %STRUCT2ARRAY_LOCAL Convert scalar structure to numeric row
        c = struct2cell(s);
        a = [c{:}];
    end

    function [out1, out2, out3] = kingraph(I, m, q, stringa)
        %KINGRAPH Plot helper for the 6 diagnostic panels.
        % Uses outer variables x, y, txtlbl.

        if I == 4
            xtick = range(log(x)) / 50;
            ytick = range(log(y)) / 50;
        else
            xtick = range(x) / 50;
            ytick = range(y) / 50;
        end

        subplot(2, 3, I)
        hold on

        if I == 1 || I == 4
            startIdx = 1;
        else
            startIdx = 3;
        end

        % Data points (blue circles)
        if I == 4
            plot(log(x(startIdx:end)), log(y(startIdx:end)), ...
                'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
                'MarkerSize', 4, 'LineStyle', 'none');
        else
            plot(x(startIdx:end), y(startIdx:end), ...
                'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
                'MarkerSize', 4, 'LineStyle', 'none');
        end

        % Axis intercepts (red circles) for linear plots
        if I ~= 1 && I ~= 4
            plot(x(1:2), y(1:2), ...
                'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', ...
                'MarkerSize', 4, 'LineStyle', 'none');
        end

        % Regression line (orange)
        if I == 1
            xs = linspace(0, max(x), 500);
            ys = (m .* xs) ./ (q + xs);
            plot(xs, ys, '-', 'Color', [0.87 0.49 0]);
        elseif I == 4
            xs = linspace(min(x), max(x), 500);
            ys = log(m) + log(xs) - log(q + xs);
            plot(log(xs), ys, '-', 'Color', [0.87 0.49 0]);
        else
            xs = linspace(min(x), max(x), 500);
            plot(xs, polyval([m q], xs), '-', 'Color', [0.87 0.49 0]);
        end

        as = axis;
        if nargout
            out1 = as;
            out2 = xtick;
            out3 = ytick;
        end

        if I == 1
            % Vmax asymptote and Km coordinates
            plot(as(1:2), [m m], 'k--')
            plot([0 q], [m m]/2, 'k--', [q q], [0 m/2], 'k--')
        elseif I == 4
            as(4) = as(4) + 10 * ytick;
            axis(as)
            plot(as(1:2), log([m m]), 'k--')
            ax = axis;
            plot([as(1) log(q)], log([m m]/2), 'k--', ...
                 log([q q]), [log(m/2) ax(3)], 'k--')
        else
            axis([as(1)-10*xtick as(2)+10*xtick as(3)-10*ytick as(4)+10*ytick])
            xk1 = [0 0];
            yk1 = as(3:4) + [-10 10] .* ytick;
            xk2 = as(1:2) + [-10 10] .* xtick;
            yk2 = [0 0];
            plot(xk1, yk1, 'k', xk2, yk2, 'k')
        end

        hold off
        axis square

        % Titles and axis labels
        title(txtlbl(I, 1))
        xlabel(txtlbl(I, 2))
        ylabel(txtlbl(I, 3))

        % Annotations for linear (non-MM) plots
        if I ~= 1 && I ~= 4
            H = text(as(1) - 5 * xtick, as(4) + 5 * ytick, stringa(1));
            set(H, 'Color', 'm')
            H = text(as(1) - 5 * xtick, as(4) - ytick, stringa(2));
            set(H, 'Color', 'm')
            H = text(5 * xtick, y(2), stringa(3));
            set(H, 'Color', 'm')

            if I == 2 || I == 3
                H = text(x(1) + 0.5 * xtick, -5 * ytick, stringa(4));
            elseif I == 5 || I == 6
                H = text(x(1) - 12 * xtick, -5 * ytick, stringa(4));
            end
            set(H, 'Color', 'm')
        end
    end

end

function app = rs_explorer_app(P, edge_factor)
%RS_EXPLORER_APP  Interactive slider app: drag the square rotation angle.
%
%   rs_explorer_app            uses params_default
%   rs_explorer_app(P)         your own parameter struct
%   rs_explorer_app(P, ef)     with the Stage 4 edge factor (S4.edge_factor)
%
%   Left: the 5x5 lattice at the chosen angle over the dashed as-printed
%   outline; ligament colour = stress status. Right: live numbers and the
%   strain-angle curve with a moving marker. Base MATLAB (uifigure).
here = fileparts(mfilename('fullpath'));
addpath(here);  addpath(genpath(fullfile(here, 'src')));
if nargin < 1 || isempty(P), P = params_default(); end
if nargin < 2, edge_factor = 1; end

g   = P.geom;
th0 = deg2rad(g.theta0_deg);
E   = hinge_modulus(P);
NL  = rs_nonlinear_response(g.a, g.Lh, g.w, g.t, E, th0, ...
        P.model.root_length_factor, P.stage2.lock_margin);
th_max = rad2deg(NL.theta_lock);
N   = P.anim.N;
L0  = rs_layout(th0, g.a, g.Lh, th0, N, P.anim.lig_points);
Lmx = rs_layout(NL.theta_lock, g.a, g.Lh, th0, N, P.anim.lig_points);
half = 0.5*(N-1)*Lmx.p + g.a;
th_curve = linspace(0, th_max, 200);
eps_curve = arrayfun(@(t) rs_state_at_theta(P, deg2rad(t)).eps, th_curve);

% ---------------- layout ------------------------------------------------------
app.fig = uifigure('Name', 'Rotating-squares explorer', 'Position', [100 100 1050 600], ...
                   'Color', 'w');
if exist('theme', 'file') || exist('theme', 'builtin'), theme(app.fig, 'light'); end
gl = uigridlayout(app.fig, [2 2], 'RowHeight', {'1x', 70}, 'ColumnWidth', {'1.2x', '1x'});
app.ax = uiaxes(gl);  app.ax.Layout.Row = 1;  app.ax.Layout.Column = 1;
right = uigridlayout(gl, [2 1], 'RowHeight', {170, '1x'});
right.Layout.Row = 1;  right.Layout.Column = 2;
app.txt = uitextarea(right, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 13);
app.ax2 = uiaxes(right);
plot(app.ax2, th_curve, 100*eps_curve, 'Color', P.plot.c_series1, 'LineWidth', 2);
hold(app.ax2, 'on');
app.marker = plot(app.ax2, g.theta0_deg, 0, 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', P.plot.c_ink, 'MarkerFaceColor', P.plot.c_good);
xlabel(app.ax2, '\theta (deg)');  ylabel(app.ax2, 'Strain (%)');  grid(app.ax2, 'on');

ctl = uigridlayout(gl, [1 2], 'ColumnWidth', {150, '1x'});
ctl.Layout.Row = 2;  ctl.Layout.Column = [1 2];
uilabel(ctl, 'Text', 'Square rotation θ (deg)', 'FontSize', 13);
app.slider = uislider(ctl, 'Limits', [0 th_max], 'Value', g.theta0_deg, ...
    'MajorTicks', 0:5:floor(th_max));
app.slider.ValueChangingFcn = @(~, evt) update(evt.Value);
app.slider.ValueChangedFcn  = @(src, ~) update(src.Value);

update(g.theta0_deg);

    function update(theta_deg)
        th = deg2rad(theta_deg);
        St = rs_state_at_theta(P, th, edge_factor);
        L  = rs_layout(th, g.a, g.Lh, th0, N, P.anim.lig_points);
        draw_lattice(app.ax, P, L, L0, St.rgb);
        xlim(app.ax, [-half half]);  ylim(app.ax, [-half half]);
        app.marker.XData = theta_deg;  app.marker.YData = 100*St.eps;
        app.marker.MarkerFaceColor = St.rgb;
        app.txt.Value = {
            sprintf('theta            %6.2f deg  (printed %g)', theta_deg, g.theta0_deg)
            sprintf('strain x = y     %+6.2f %%', 100*St.eps)
            sprintf('porosity         %6.1f %%', 100*St.porosity)
            sprintf('nu ideal / beam  %6.3f / %6.4f', St.nu_ideal, St.nu_compliant)
            sprintf('lattice stress   %+6.2f MPa', St.sigma)
            sprintf('peak hinge stress %s', St.sig_text)
            sprintf('status           %s', St.status)};
    end
end

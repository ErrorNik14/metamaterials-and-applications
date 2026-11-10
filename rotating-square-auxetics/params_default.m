function P = params_default()
%PARAMS_DEFAULT  Every input to the rotating-squares project lives here.
%
%   P = params_default() returns one struct P. No other file contains a
%   physical number: change a value here and every stage picks it up.
%
%   UNITS (consistent set, used everywhere):
%       length  -> mm
%       force   -> N
%       stress  -> MPa  (= N/mm^2)
%       angle   -> degrees in this file, converted to radians in the code
%
%   To try a variant without editing this file, do e.g.
%       P = params_default();  P.geom.a = 12;  S1 = stage1_run(P);

%% ===================== 1. GEOMETRY OF ONE CELL (mm) =====================
P.geom.a          = 10.0;   % side length of each square
P.geom.Lh         = 1.0;    % hinge ligament length (gap bridged between corners)
P.geom.w          = 0.8;    % hinge ligament width (in-plane, bending direction)
P.geom.t          = 5.0;    % thickness out of plane (extrusion depth)
P.geom.theta0_deg = 15;     % as-printed rotation angle of the squares
                            %  0  = squares touching (fully closed)
                            %  45 = fully open (diamond holes, kinematic lock)

%% ===================== 2. MATERIAL (PLA, MPa) ===========================
P.mat.name           = 'PLA';
P.mat.E              = 3500;   % Young's modulus, datasheet [MPa]
P.mat.sigma_y        = 50;     % yield / tensile strength [MPa]
P.mat.nu             = 0.36;   % Poisson's ratio of the bulk PLA
P.mat.print_knockdown = 1.0;   % multiply E and sigma_y by this for FDM parts
                               % (try 0.7-0.85 for printed PLA; 1.0 = datasheet)

%% ===================== 3. 3D PRINTER ====================================
P.print.nozzle       = 0.4;    % nozzle diameter [mm]
P.print.layer_height = 0.2;    % layer height [mm]
P.print.min_perims   = 2;      % minimum extrusion widths across a thin wall

%% ===================== 3b. HINGE MODEL CHOICES (Stage 2) ================
P.model.hinge_E_mode = 'plane_stress';  % 'plane_stress' -> use E
                                        % 'wide_beam'    -> use E/(1-nu^2)
                                        %   (ligament depth t >> width w)
P.model.root_length_factor = 0;  % effective ligament length = Lh + factor*w.
                                 % 0 = pure beam theory. Raise it (e.g. 0.5)
                                 % only to calibrate against Ansys/test data.

%% ===================== 3c. CAD / EXPORT (Stage 6) =======================
P.cad.N            = 5;      % squares per side of the printed part
P.cad.lig_overlap  = 0.5;    % ligament runs this many widths w into each square
                             % so the union is one solid (not a point contact)
P.cad.min_gap_factor = 2;    % warn if a pore/slot is narrower than this x nozzle
P.export.strain_grid = -0.02:0.0005:0.02;   % strains in the results CSV

%% ===================== 4. FINITE-ELEMENT CHECK (Stage 3) ================
P.fe.N_cells       = 5;        % array is N x N squares
P.fe.applied_strain = 0.01;    % nominal |strain| applied in tension/compression
P.fe.n_el_lig      = 4;        % beam elements along each ligament
P.fe.square_E_factor = 1e3;    % squares modelled as beams this many times
                               % stiffer than PLA (=> effectively rigid)
P.fe.N_convergence = [3 5 7 9 11 13];        % array sizes for edge-effect study
P.fe.w_compare     = [0.4 0.8 1.2 1.6 2.0 2.5 3.0];  % widths for the comparison
P.fe.plot_scale    = 5;        % displacement magnification in the deformed plot
P.fe.tol_nu_vs_S2  = 0.01;     % accepted model-to-model difference (FE vs Stage 2)

%% ===================== 4b. FAILURE ANALYSIS (Stage 4) ===================
P.fail.safety_factor = 1.5;    % on stress (yield) and on the instability strain
P.fail.Kt            = 1.0;    % stress-concentration factor at the ligament root
                               % 1.0 = ideal beam theory. A sharp re-entrant
                               % corner raises it; fillets lower it. Take a
                               % value from Ansys (peak / nominal stress).
P.fail.buckle_ref_strain = -0.01;  % compressive pre-strain for the buckling solve
P.fail.n_modes       = 3;          % buckling modes to report
P.fail.Lh_design     = [0.5 4];    % ligament-length range of the design chart [mm]
P.fail.w_design      = [0.8 1.2 1.6];  % widths drawn in the design chart [mm]

%% ===================== 5. ANALYSIS RANGES ===============================
P.stage1.theta_range_deg = [0 60]; % plotted range (beyond 45 deg the cell re-closes)
P.stage1.n_theta         = 241;    % number of angle samples
P.stage1.snapshot_deg    = [0 15 30 45];  % angles shown in the geometry figure
P.stage1.n_vis           = 4;      % squares per side in the geometry figure
P.stage1.n_check         = 6;      % squares per side used by the numerical check

P.stage2.w_range   = [0.2 3.0];  % ligament width sweep [mm]
P.stage2.Lh_range  = [0.4 4.0];  % ligament length sweep [mm]
P.stage2.a_range   = [5 20];     % square side sweep [mm]
P.stage2.n_sweep   = 120;        % points per sweep
P.stage2.a_family  = [5 10 20];  % curves drawn in the nu-vs-w figure [mm]
P.stage2.lock_margin = 0.05;     % stop the nonlinear curve when the lever
                                 % arm dp/dtheta falls below 5 % of its start

%% ===================== 6. NUMERICAL TOLERANCES ==========================
P.tol.geom = 1e-9;   % hinge gap allowed [mm] (should be ~1e-14 in practice)
P.tol.nu   = 1e-6;   % allowed |nu - (-1)| for the ideal model
P.tol.rel  = 1e-6;   % relative agreement between two routes to the same number

%% ===================== 7. PLOT STYLE (slide-ready) ======================
P.plot.font_name   = 'Arial';
P.plot.font_size   = 11;       % axis tick labels [pt]
P.plot.title_size  = 12;
P.plot.line_width  = 2.0;
P.plot.marker_size = 6;
P.plot.fig_width_cm  = 16;     % one full-width slide figure
P.plot.fig_height_cm = 9;
P.plot.dpi         = 300;      % PNG export resolution
P.plot.save_fig    = true;     % also save editable .fig files
P.plot.save_vector = true;     % also save a vector PDF (sharpest in slides)
% Colour roles (fixed order: same quantity = same colour in every figure)
P.plot.c_series1 = hex2rgb_local('#2a78d6');  % blue   - main quantity (x)
P.plot.c_series2 = hex2rgb_local('#eb6834');  % orange - second quantity (y)
P.plot.c_series3 = hex2rgb_local('#1baf7a');  % aqua   - numerical check
P.plot.c_ink     = hex2rgb_local('#0b0b0b');  % text
P.plot.c_muted   = hex2rgb_local('#52514e');  % secondary text, reference lines
P.plot.c_square  = hex2rgb_local('#86b6ef');  % fill of the squares
P.plot.c_ghost   = hex2rgb_local('#b4b2a9');  % undeformed outline
P.plot.c_lock    = hex2rgb_local('#d03b3b');  % locking limits (status: critical)
P.plot.c_good    = hex2rgb_local('#0ca30c');  % status: inside safe range
P.plot.c_warn    = hex2rgb_local('#fab219');  % status: past safe, before yield

%% ===================== 7b. ANIMATION & APP (Stage 5) ====================
P.anim.N          = 5;       % squares per side in the animation
P.anim.n_frames   = 120;     % frames in one open-close cycle
P.anim.fps        = 24;      % frames per second (GIF and MP4)
P.anim.open_frac  = 0.97;    % go this fraction of the way to the open lock
P.anim.close_deg  = 0;       % most-closed angle shown (deg)
P.anim.size_cm    = [16 10]; % frame size
P.anim.lig_points = 12;      % points drawn along each (curved) ligament

%% ===================== 8. OUTPUT FOLDERS ================================
root = fileparts(mfilename('fullpath'));
P.paths.root    = root;
P.paths.results = fullfile(root, 'results');
P.paths.figs    = fullfile(root, 'results', 'figures');
P.paths.data    = fullfile(root, 'results', 'data');
P.paths.anim    = fullfile(root, 'results', 'animation');
end

function rgb = hex2rgb_local(hex)
% '#RRGGBB' -> [r g b] in 0..1 (base MATLAB, no toolbox).
hex = strrep(hex, '#', '');
rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))] / 255;
end

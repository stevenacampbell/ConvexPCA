%% Code for a runtime comparison of [C & W (2022)] with [Cazelles et al. (2018)]

% All of the necessary MATLAB code is taken from the repository:
%   -"https://github.com/ecazelles/2017-GPCA-vs-LogPCA-Wasserstein/tree/master"
% NOTE: This code file will NOT run unless their code is downloaded and
% placed in the same directory.

%% Set Reference Directory
addpath('toolbox')
clear all
close all

%% Load Data
% NOTE: The code file "DataPrep_Cazellesetal_Comparison.R" must be run in
% order to produce this data.
Omega = csvread('omega_grid.csv', 1, 0);
mu = csvread('histograms.csv', 1, 0);

% Reformat and rescale for compatibility purposes
Omega = Omega.';
domega = max(diff(Omega));
Omega = Omega /domega;
mu = mu * domega;

% Plot data
figure(1)
plot(Omega,mu)
title('Data');
%% Barycenter

method = 'pchip';
n_inv = 10000; 
[Bs, FBs] = wasserstein_barycenter_1D_smooth(mu,Omega,method,n_inv);
f = [Bs(1) Bs]; % Smooth histogram of the barycenter

%% Log Maps

OmegaExt = [Omega(1)-1, Omega];
V = zeros(size(mu,1),length(Omega)+1);
for i = 1:size(mu,1)
    V(i,:) = logMap(mu(i,:),FBs,Omega,method); % log maps of the data at the barycenter
end

%% Compute PCA on logmap (log-PCA)

Vc = mean(V,1);
Vp = bsxfun(@minus,V,Vc)*diag(sqrt(f));
C = Vp'*Vp/size(Vp,1);
[eigV, eigVals] = eig(C);
nonzero_ind = (f > 0);
eigV(nonzero_ind,:) = diag(1./sqrt(f(nonzero_ind))) * eigV(nonzero_ind,:);
eigV = eigV(:,end:-1:1);

%% GPCA - Iterative Geodesic Approach
 
% Run optimzation and output computational cost
tic
L = 2;
% Choose initialization
%V0 = rand(L,length(Omega)+1); % random
V0 = eigV(:,1:2).';
range_t0=[-1:0.01:1];
% Run optimization
[v_gpca_iter, t_gpca_iter,t0_iter,residual_iter,W_residual_iter] = algo_GPCA_1D_iter(V,OmegaExt,L,V0,f,range_t0);
toc

% Representation of the 1st, 2nd components of the iterative geodesic 
% approach of Cazelles et al.
% CREDIT: All plotting functionality is taken from their repository.
figure(2)
map1=autumn(size(mu,1));
map2=cool(size(mu,1));

% 1st component (Note: The main effect of perturbing the variance is captured)
subaxis(1,2,1,'SpacingVert',0.07,'SpacingHoriz',0.07,'ML',0.15);
h_iter = zeros(size(mu,1),length(OmegaExt));
[~, I] = sort(t_gpca_iter(:,1));
for i=1:size(mu,1)
    T_iter = OmegaExt + t_gpca_iter(I(i),1) * v_gpca_iter(1,:);
    h_iter(i,:) = pushforward_density(T_iter, f, OmegaExt);
    plot(OmegaExt, h_iter(i,:), 'Color', map1(i,:))
    axis([Omega(1) Omega(end) -0.002 0.15])
    hold on
end
pl = plot(OmegaExt, f, '-k', 'linewidth', 1);
legend(pl, 'Wasserstein barycenter');
ylabel('Iterative Geodesic approach', 'FontSize', 10, 'FontWeight', 'bold');
title('First PG', 'FontSize', 10, 'FontWeight', 'bold');

% 2nd component (Note: The tail effects represented by the second PC
% are not visualized very well when using histograms. 
% The influence is very small in the center of the distribution and so most
% of the projections onto the second component overlap.
subaxis(1,2,2,'SpacingVert',0.07,'SpacingHoriz',0.07,'ML',0.15);
[~, I] = sort(t_gpca_iter(:,2));
for i=1:size(mu,1)
    T_iter = OmegaExt + t_gpca_iter(I(i),2) * v_gpca_iter(2,:);
    h_iter(i,:) = pushforward_density(T_iter, f, OmegaExt);
    plot(OmegaExt, h_iter(i,:), 'Color', map2(i,:))
    axis([Omega(1) Omega(end) -0.002 0.15])
    hold on
end
pl = plot(OmegaExt, f, '-k', 'linewidth', 1);
legend(pl, 'Wasserstein barycenter');
title('Second PG', 'FontSize', 10, 'FontWeight', 'bold');

drawnow;

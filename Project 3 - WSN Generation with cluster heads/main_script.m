clear; clc; close all;
rng(24); %set the random seed!

%% WSN Parameters
N = 100; %The number of sensors
area_width = 100;  area_height = 100;
eta_short = 10;     %nJ/bit/m2;
eta_long = 0.0013;  %nJ/bit/m4;
d0 = sqrt(eta_short/eta_long); 
E_initial = 2e9; % in nJ (the starting energy)
E_elec = 50; E_agg = 50;
k = 625*8; %number_of_bits per cycle
R = 30;  %threshold for one-hop/dual-hop transmission

func_tx_energy = @(E,eta,exp,d) max(0, E - (k*(E_elec + E_agg) + k*eta*d.^exp));
func_rx_energy = @(E, n) max(0, E - (n*k*E_elec));

%% Generating the location of sensors
locs = cat(1, randperm(area_width), randperm(area_height))';
% assert no node existance @ the network center
for i=1:size(locs,2)
    if locs(i,:) == [50 50]
        locs(i,:) = locs(i,:) + randi([1,20], 1 , 2);
    end
end
% Add the sink node @ the network center
f1 = figure('Name', 'Network Topology - - sink centered');
scatter(locs(:,1),locs(:,2));
hold on
scatter(50,50,'filled');
locs = cat(1, locs, [50, 50]);
legend('source', 'sink');
title('Network Topology - sink centered');
xlabel('distance on x-axis');
ylabel('distance on y-axis');
% saveas(f1, [pwd '/Figures/network_topology_Sink_Centered']);
%% Calculate the distances
dists = sqrt(sum((locs-[50 50]).^2,2));
%% Sort Nodes according to Distances
[dists, ids] = sort(dists);
locs = locs(ids,:);
%% remove sink to end of vectors
dists = [dists(2:end); dists(1)];
locs = [locs(2:end,:); locs(1,:)];

%% GO on cycles!
    % Initiate the energy
energy = E_initial*ones(1,N+1);
active_nodes =(N);
[assigned_heads,result_energies]  = CH_election(5, 5, energy(1:end-1), locs(1:end-1,:), dists(1:end-1), d0);
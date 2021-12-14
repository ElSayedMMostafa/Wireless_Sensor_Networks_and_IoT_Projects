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
R = 10; %threshold for one-hop/dual-hop transmission
func_tx_energy = @(E,eta,exp,d) max(0, E - (k*(E_elec + E_agg) + k*eta*d^exp));
func_rx_energy = @(E, n) max(0, E - (n*k*E_elec));

%% Generating the location of senssors
locs = cat(1, randperm(area_width), randperm(area_height))';
% assert no node existance @ the network center
for i=1:size(locs,2)
    if locs(i,:) == [50 50]
        locs(i,:) = locs(i,:) + randi([1,20], 1 , 2);
    end
end
% Add the sink node @ the network center
f1 = figure('Name', 'Network Topology');
scatter(locs(:,1),locs(:,2));
hold on
scatter(50,50,'filled');
locs = cat(1, locs, [50, 50]);
legend('source', 'sink');
title('Network Topology');
xlabel('distance on x-axis');
ylabel('distance on y-axis');
% saveas(f1, [pwd '/Figures/network_topology_Sink_Centered']);
%% Calculate the distances
dists = sqrt(sum((locs-[50 50]).^2,2));
%% Initiate the energy
energy = E_initial*ones(1,N+1);


%% GO on cycles!
active_nodes =(N);

while 1
    % Calculate the dissipated & remaining energy
for i=1:N
    dist = dists(i);
    if dist > R
        [sending_index, dist_to_node] = intermediate_selection(i, locs(1:end-1, :), dists(1:end-1));
        % Transmission between node and intermediate
        if dist_to_node <= d0 
            energy(i) = func_tx_energy(energy(i), eta_short, 2, dist_to_node);
        else
            energy(i) = func_tx_energy(energy(i), eta_long, 4, dist_to_node);
        end
        energy(sending_index) = func_rx_energy(energy(sending_index), 1);
        % Transmission between intermediate and sink
        dist = dists(sending_index);
        if dist <= d0 
            energy(sending_index) = func_tx_energy(energy(sending_index), eta_short, 2, dist);
        else
            energy(sending_index) = func_tx_energy(energy(sending_index), eta_long, 4, dist);
        end  
    % Otherwise, THE NORMAL CASE
    elseif dist <= d0
        energy(i) = func_tx_energy(energy(i), eta_short, 2, dist);
    else
        energy(i) = func_tx_energy(energy(i), eta_long, 4, dist);
    end
end
energy(N+1) = func_rx_energy(energy(N+1), active_nodes(end));
    % Get the no. of active nodes
tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
active_nodes = cat(2, active_nodes, sum(energy(1:N) >= tx_threshold));

rx_threshold = k*E_elec;
if active_nodes(end) == 0 || energy(N+1) < rx_threshold
    break;
end
end

% remove the first element in the active_nodes vector (intialization)
 active_nodes = active_nodes(2:end);
f2 = figure('Name', 'Curve of Active Nodes');
plot(linspace(1,length(active_nodes),length(active_nodes)), active_nodes);
grid on
xlabel('no. of cycle');
ylabel('no. of active nodes');
title('The number of active nodes at each cycle - Sink Centered');

%% TRIAL DEMO -- In case of dual hop transmission, select the intermediate node! 
min_index = intermediate_selection(1, locs(1:end-1, :), dists(1:end-1));
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
k = 625*1024; %number_of_bits per cycle
func_tx_energy = @(E,eta,exp,d) max(0, E - (k*(E_elec + E_agg) + k*eta*d^exp));
func_rx_energy = @(E) max(0, E - (k*(E_elec+E_agg)));

%% Generating the location of senssors
locs = cat(1, randperm(area_width), randperm(area_height))';
% assert no node existance @ the network center
for i=1:size(locs,2)
    if locs(i,:) == [50 50]
        locs(i,:) = locs(i,:) + randi([1,20], 1 , 2);
    end
end
% Add the sink node @ the network center
figure('Name', 'Network Topology')
scatter(locs(:,1),locs(:,2));
hold on
scatter(50,50,'filled');
locs = cat(1, locs, [50, 50]);
legend('source', 'sink');
%% Calculate the distances
dists = sqrt(sum((locs-[50 50]).^2,2));
%% Initiate the energy
energy = E_initial*ones(1,N+1);
%% GO on cycles!
active_nodes =[];

while 1
    % Calculate the dissipated & remaining energy
for i=1:N
    dist = dists(i);
    if dist <= d0
        energy(i) = func_tx_energy(energy(i), eta_short, 2, dist);
    else
        energy(i) = func_tx_energy(energy(i), eta_long, 4, dist);
    end
end
energy(N+1) = func_rx_energy(energy(N+1));
    % Get the no. of active nodes
    rx_threshold = k*(E_elec+E_agg);
    tx_threshold = rx_threshold + k*eta_short*0.5^2;
active_nodes = cat(2, active_nodes, sum(energy(1:N) >= tx_threshold));

if active_nodes(end) == 0 || energy(N+1) < rx_threshold
    break;
end
end
    
figure('Name', 'Curve of Active Nodes')
plot(linspace(1,length(active_nodes),length(active_nodes)), active_nodes);
grid on
xlabel('no. of cycle');
ylabel('no. of active nodes');
title('The number of active nodes at each cycle');


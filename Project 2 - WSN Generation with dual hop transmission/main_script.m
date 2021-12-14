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
R = 30; %threshold for one-hop/dual-hop transmission

func_tx_energy = @(E,eta,exp,d) max(0, E - (k*(E_elec + E_agg) + k*eta*d^exp));
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
%energy = [];
%% GO on cycles!
active_nodes =(N);
t = 1; %initial time, step
while 1
t = t+1; % time step
energy = cat(1, energy, E_initial*ones(1,N+1));
% Calculate the dissipated & remaining energy
for i=1:N
    dist = dists(i);
    if dist > R
        [sending_index, dist_to_node] = intermediate_selection(i, locs(1:end-1, :), dists(1:end-1));
        % Transmission between node and intermediate
        if dist_to_node <= d0 
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_short, 2, dist_to_node);
        else
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_long, 4, dist_to_node);
        end
        energy(t,sending_index) = func_rx_energy(energy(t-1,sending_index), 1);
        % Transmission between intermediate and sink
        dist = dists(sending_index);
        if dist <= d0 
            energy(t,sending_index) = func_tx_energy(energy(t,sending_index), eta_short, 2, dist);
        else
            energy(t,sending_index) = func_tx_energy(energy(t,sending_index), eta_long, 4, dist);
        end  
    % Otherwise, THE NORMAL CASE
    elseif dist <= d0
        energy(t,i) = func_tx_energy(energy(t-1,i), eta_short, 2, dist);
    else
        energy(t,i) = func_tx_energy(energy(t-1,i), eta_long, 4, dist);
    end
end
energy(t,N+1) = func_rx_energy(energy(t-1,N+1), active_nodes(end));
    % Get the no. of active nodes
tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
active_nodes = cat(2, active_nodes, sum(energy(t,1:N) >= tx_threshold));

rx_threshold = k*E_elec;
if active_nodes(end) == 0 || energy(t,N+1) < rx_threshold
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
title('The number of active nodes at each cycle');

%% Find the lifetime T1
T1 = find(active_nodes < N, 1);

f3= figure('Name', 'Curve of Active Nodes with 1st lifetime');
plot(linspace(1,length(active_nodes),length(active_nodes)), active_nodes);
grid on
hold on
plot(T1, active_nodes(T1),'r*');

xlabel('no. of cycle');
ylabel('no. of active nodes');
legend('no. active nodes','T1');
title('The number of active nodes at each cycle');
%saveas(f3, [pwd '/Figures/curve of active nodes with 1st lifetimes']);

%% Plotting the remaining energies at T1 cycles
energy_T1 = energy(T1+1, :);
f4 = figure('Name','Remaining energies of the N nodes after T1 cycles');
stem(linspace(1,N,N),energy_T1(1:end-1),"filled")
hold on
stem(N+1,energy_T1(end),"filled",'LineWidth',1,'Color','r')
hold off
grid on
xlabel('Node Index');
ylabel('Energy (nJ)');
title('Remaining energies of the N nodes after T1 cycles')
%saveas(f4, [pwd '/Figures/Remaining energies of the N nodes after T1 cycles'])

%% Relation between R and T1
% We're going to repeat the previous procedure with a range of R values
% Rmax = 25*sqrt(2) since that the area is a square 100*100
Rmax = 25*sqrt(2);
Rs = 10:10:Rmax;
T1s = [];  % a vector to save T1 values
for R = Rs
    % (1) Energy Initiation
energy = E_initial*ones(1,N+1);

    % (2) Go on Tx & Rx cycles
active_nodes =(N);
t = 1; %initial time, step
while 1
t = t+1; % time step
energy = cat(1, energy, E_initial*ones(1,N+1));
% Calculate the dissipated & remaining energy
for i=1:N
    dist = dists(i);
    if dist > R
        [sending_index, dist_to_node] = intermediate_selection(i, locs(1:end-1, :), dists(1:end-1));
        % Transmission between node and intermediate
        if dist_to_node <= d0 
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_short, 2, dist_to_node);
        else
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_long, 4, dist_to_node);
        end
        energy(t,sending_index) = func_rx_energy(energy(t-1,sending_index), 1);
        % Transmission between intermediate and sink
        dist = dists(sending_index);
        if dist <= d0 
            energy(t,sending_index) = func_tx_energy(energy(t,sending_index), eta_short, 2, dist);
        else
            energy(t,sending_index) = func_tx_energy(energy(t,sending_index), eta_long, 4, dist);
        end  
    % Otherwise, THE NORMAL CASE
    elseif dist <= d0
        energy(t,i) = func_tx_energy(energy(t-1,i), eta_short, 2, dist);
    else
        energy(t,i) = func_tx_energy(energy(t-1,i), eta_long, 4, dist);
    end
end
energy(t,N+1) = func_rx_energy(energy(t-1,N+1), active_nodes(end));
    % Get the no. of active nodes
tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
active_nodes = cat(2, active_nodes, sum(energy(t,1:N) >= tx_threshold));

rx_threshold = k*E_elec;
if active_nodes(end) == 0 || energy(t,N+1) < rx_threshold
    break;
end
end

    % (3) Get T1 for R
    T1 = find(active_nodes < N, 1);
    T1s = cat(2, T1s, T1);
end
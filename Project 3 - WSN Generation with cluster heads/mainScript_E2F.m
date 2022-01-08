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
R = 25;  % radius of the cluster heads

%% Definition of tx and Rx functions
% Tx for regular nodes
func_tx_energy = @(E,eta,exp,d) max(0, E - (k*(E_elec + E_agg) + k*eta*d.^exp));
% Rx from regular nodes
func_rx_energy = @(E, n) max(0, E - (n*k*E_elec));

%% Generating the location of sensors and cluster heads
% (1) Sensor Node
locs = cat(1, randperm(area_width), randperm(area_height))';
% (2) Cluster Heads
uniform_angles = linspace(0,2*pi-2*pi/5,5);
heads_locs = [R*cos(uniform_angles) + 50; R*sin(uniform_angles) + 50]';
% Add the sink node @ the network center
f1 = figure('Name', 'Network Topology with special heads');
scatter(locs(:,1),locs(:,2));
hold on
scatter(heads_locs(:,1), heads_locs(:,2),'filled', 'blue');
hold on
scatter(50,50,'filled');
locs = cat(1, locs, [50, 50]);
legend('source', 'heads','sink');
title('Network Topology with special heads');
xlabel('distance on x-axis');
ylabel('distance on y-axis');
 saveas(f1, [pwd '/Figures/network_topology_with_special_heads']);
%% Calculate the distances
dists = sqrt(sum((locs-[50 50]).^2,2));
heads_dists = sqrt(sum((heads_locs-[50 50]).^2,2));
%% Sort Nodes according to Distances
[dists, ids] = sort(dists);
locs = locs(ids,:);
%% remove sink to end of vectors
dists = [dists(2:end); dists(1)];
locs = [locs(2:end,:); locs(1,:)];
%% GO on cycles!
    % Initiate the energy
energy = E_initial*ones(1,N+1);
heads_energy = 2*E_initial*ones(1,5);
active_nodes_count =(N); active_heads_count = (5);
active_nodes = (1:N); active_heads = (1:5);
tx_threshold = 0;
rx_threshold = k*E_elec; %recieve from @ least 1 head
heads_energy_threshold = tx_threshold + rx_threshold;
t = 1; %initial time step
while 1
 t = t+1;
 if active_heads_count(end) > 0
   num_heads = active_heads_count(end);
   energy = cat(1, energy, zeros(1,N+1));
   heads_energy = cat(1, heads_energy, zeros(1,5));
   for node_index = 1:N
        % select the nearest head
        [dist_to_head, selected_head] = min(sqrt(sum((heads_locs-locs(node_index,:)).^2,2)));
        % (1) nodes send to heads
        if dist_to_head >= d0
            energy(t,node_index) = func_tx_energy(energy(t-1,node_index), eta_long, 4, dist_to_head);
        else
            energy(t,node_index) = func_tx_energy(energy(t-1,node_index), eta_short, 2, dist_to_head);
        end
        % (2) heads receives from nodes
        heads_energy(t,selected_head) = func_rx_energy(heads_energy(t-1,selected_head), 1);
        % (3) heads sends to sink
        head_dist = heads_dists(selected_head);
        if head_dist >= d0
            heads_energy(t, selected_head) = func_tx_energy(heads_energy(t, selected_head), eta_long, 4, head_dist);
        else
            heads_energy(t, selected_head) = func_tx_energy(heads_energy(t, selected_head), eta_short, 2, head_dist);
        end
        % (4) sink receives from heads
        energy(t,N+1) = func_rx_energy(energy(t-1,N+1), 1);
   end
   % Check energies!
    active_nodes = find(energy(t,1:N) > tx_threshold);
    active_nodes_count = cat(2, active_nodes_count, sum(energy(t,1:N) > tx_threshold));
    
    if active_nodes_count(end) == 0 || energy(t,N+1) < rx_threshold 
        break;
    end
    % Check the remaining heads
    active_heads = find(heads_energy(t,1:5) > heads_energy_threshold);
    active_heads_count = cat(2, active_heads_count, sum(heads_energy(t,1:5) > heads_energy_threshold));
T_end_clusters = t; % save this for later!
 else %else_if
    % Operate in the very normal situation!
    for i=1:N
        dist = dists(i);
        if dist <= d0
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_short, 2, dist);
        else
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_long, 4, dist);
        end
    end
    energy(t,N+1) = func_rx_energy(energy(t-1,N+1), active_nodes(end));

    % Check energies!
    active_nodes = find(energy(t,1:N) > tx_threshold);
    active_nodes_count = cat(2, active_nodes_count, sum(energy(t,1:N) > tx_threshold));
    rx_threshold = k*E_elec; %recieve from @ least 1 node
    if active_nodes_count(end) == 0 || energy(t,N+1) < rx_threshold
        break;
    end
end %end_if
end

%% Plotting the number of active nodes per cycle
% remove the first element in the active_nodes vector (intialization)
active_nodes_count = active_nodes_count(2:end);
f2 = figure('Name', 'Curve of Active Nodes');
plot(linspace(1,length(active_nodes_count),length(active_nodes_count)), active_nodes_count);
grid on
xlabel('no. of cycle');
ylabel('no. of active nodes');
title('The number of active nodes at each cycle');
saveas(f2, [pwd '/Figures/The number of active nodes at each cycle - Dedicated Heads'],'fig')
%% Part C: Identifing T1
% Find the lifetime T1
T1 = find(active_nodes_count < N, 1);
f3= figure('Name', 'Curve of Active Nodes with 1st lifetime');
plot(linspace(1,length(active_nodes_count),length(active_nodes_count)), active_nodes_count);
grid on
hold on
plot(T1, active_nodes_count(T1),'r*');
xlabel('no. of cycle');
ylabel('no. of active nodes');
legend('no. active nodes','T_1');
title('The number of active nodes at each cycle');
saveas(f3, [pwd '/Figures/curve of active nodes with 1st lifetimes - Dedicated Heads']);

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
title('Remaining energies of the N nodes after T_1 cycles')
saveas(f4, [pwd '/Figures/Remaining energies of the N nodes after T1 cycles - Dedicated Heads'])

%% Part D: Optimizing R parameter
Rs = 5:1:50*sqrt(2);
T1s= [];
%% GO on cycles!
for R = Rs
heads_locs = [R*cos(uniform_angles) + 50; R*sin(uniform_angles) + 50]';
heads_dists = sqrt(sum((heads_locs-[50 50]).^2,2));
    % Initiate the energy
energy = E_initial*ones(1,N+1);
heads_energy = 2*E_initial*ones(1,5);
active_nodes_count =(N); active_heads_count = (5);
active_nodes = (1:N); active_heads = (1:5);
tx_threshold = 0;
rx_threshold = k*E_elec; %recieve from @ least 1 head
heads_energy_threshold = tx_threshold + rx_threshold;
t = 1; %initial time step
while 1
 t = t+1;
 if active_heads_count(end) > 0
   num_heads = active_heads_count(end);
   energy = cat(1, energy, zeros(1,N+1));
   heads_energy = cat(1, heads_energy, zeros(1,5));
   heads_reputation = zeros(1,5);
   for node_index = 1:N
        % select the nearest head
        [dist_to_head, selected_head] = min(sqrt(sum((heads_locs-locs(node_index,:)).^2,2)));
        % (1) nodes send to heads
        if dist_to_head >= d0
            energy(t,node_index) = func_tx_energy(energy(t-1,node_index), eta_long, 4, dist_to_head);
        else
            energy(t,node_index) = func_tx_energy(energy(t-1,node_index), eta_short, 2, dist_to_head);
        end
        heads_reputation(1,selected_head) = heads_reputation(1,selected_head) + 1;
   end
   % (2) heads receives from nodes
   heads_energy(t,:) = func_rx_energy(heads_energy(t-1,:), heads_reputation);
   % (3) heads sends to sink
   for head_index = 1:5
        head_dist = heads_dists(head_index);
        if head_dist >= d0
            heads_energy(t, head_index) = func_tx_energy(heads_energy(t, head_index), eta_long, 4, head_dist);
        else
            heads_energy(t, head_index) = func_tx_energy(heads_energy(t, head_index), eta_short, 2, head_dist);
        end
   end
   % (4) sink receives from heads
   energy(t,N+1) = func_rx_energy(energy(t-1,N+1), active_heads_count(end));
   % Check energies!
    active_nodes = find(energy(t,1:N) > tx_threshold);
    active_nodes_count = cat(2, active_nodes_count, sum(energy(t,1:N) > tx_threshold));
    
    if active_nodes_count(end) == 0 || energy(t,N+1) < rx_threshold 
        break;
    end
    % Check the remaining heads
    active_heads = find(heads_energy(t,1:5) > heads_energy_threshold);
    active_heads_count = cat(2, active_heads_count, sum(heads_energy(t,1:5) > heads_energy_threshold));
 else %else_if
    % Operate in the very normal situation!
    for i=1:N
        dist = dists(i);
        if dist <= d0
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_short, 2, dist);
        else
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_long, 4, dist);
        end
    end
    energy(t,N+1) = func_rx_energy(energy(t-1,N+1), active_nodes(end));

    % Check energies!
    active_nodes = find(energy(t,1:N) > tx_threshold);
    active_nodes_count = cat(2, active_nodes_count, sum(energy(t,1:N) > tx_threshold));
    rx_threshold = k*E_elec; %recieve from @ least 1 node
    if active_nodes_count(end) == 0 || energy(t,N+1) < rx_threshold
        break;
    end
end %end_if
end
T1 = find(active_nodes_count < N, 1);
T1s = cat(2,T1s, T1);
end % end_Cs
% Plotting te relation
f5 = figure('Name','Relation between T1 and R');
stem(Rs, T1s);
grid on
xlabel('R');
ylabel('T_1');
title('Relation between T_1 and R')
saveas(f5, [pwd '/Figures/Relation between T1 and R'])


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
C = 5; portion = 5;
func_tx_energy = @(E,eta,exp,d) max(0, E - (k*(E_elec + E_agg) + k*eta*d.^exp));
func_tx_energy_head = @(E,eta,exp,d) max(0, E - ((500*8)*(E_elec + E_agg) + (500*8)*eta*d.^exp));
func_rx_energy = @(E, n) max(0, E - (n*(500*8)*E_elec));
func_rx_energy_normal = @(E, n) max(0, E - (n*k*E_elec));

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
tx_threshold = 0; %energy threshold for transmitters
active_nodes_count =(N);
active_nodes = (1:N);
t = 1; %initial time step
while 1
energy = cat(1, energy, zeros(1,N+1));
    % Elect the heads and decipate their enegies
num_heads = floor(portion/100 * active_nodes_count(end));
if num_heads > 0 %condition to go into the election mode
assigned_heads  = CH_election(C, active_nodes, num_heads, energy(t,1:end-1), locs(1:end-1,:), dists(1:end-1), d0);
heads = cell2mat(assigned_heads(:,1));
for j=1:C
    t = t + 1;
for head_index = 1:length(heads)
    head = heads(head_index);
    head_loc = locs(head);
    head_dist = dists(head);
    % (1) nodes send to heads
    for node_index = cell2mat(assigned_heads(head_index,3))
        node = locs(node_index);
        distance = sqrt(sum((head_loc-node).^2,2));
        if distance >= d0
            energy(t,node_index) = func_tx_energy(energy(t-1,node_index), eta_long, 4, distance);
        else
            energy(t,node_index) = func_tx_energy(energy(t-1,node_index), eta_short, 2, distance);
        end
    end
    % (2) heads receives from nodes
    energy(t, head) = func_rx_energy_normal(energy(t-1, head), assigned_heads{head_index,2});
    % (3) heads sends to sink
    if head_dist >= d0
        energy(t, head) = func_tx_energy_head(energy(t, head), eta_long, 4, head_dist);
    else
        energy(t, head) = func_tx_energy_head(energy(t, head), eta_short, 2, head_dist);
    end
end
     % (4) sink receives from heads
    energy(t,N+1) = func_rx_energy(energy(t-1,N+1), length(heads));

    % Check energies!
    active_nodes = find(energy(t,1:N) > tx_threshold);
    active_nodes_count = cat(2, active_nodes_count, sum(energy(t,1:N) > tx_threshold));
    rx_threshold = 500*8*E_elec; %recieve from @ least 1 head
    if active_nodes_count(end) == 0 || energy(t,N+1) < rx_threshold
        break;
    end
T_end_clusters = t; % save this for later!
end
else %end_if
% Operate in the very normal situation!
% C = 1;
t = t + 1;
    for i=1:N
        dist = dists(i);
        if dist <= d0
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_short, 2, dist);
        else
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_long, 4, dist);
        end
    end
    energy(t,N+1) = func_rx_energy_normal(energy(t-1,N+1), active_nodes(end));

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
ylabel('no. of active nodes -sink centered');
title('The number of active nodes at each cycle');

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
title('Remaining energies of the N nodes after T1 cycles - sink centered')
%saveas(f4, [pwd '/Figures/Remaining energies of the N nodes after T1 cycles'])

%% Part D: Optimizing C parameter
Cs = 2:1:40;
T1s= [];
for C = Cs
%% GO on cycles!
    % Initiate the energy
energy = E_initial*ones(1,N+1);
tx_threshold = 0; %energy threshold for transmitters
active_nodes_count =(N);
active_nodes = (1:N);
t = 1; %initial time step
while 1
energy = cat(1, energy, zeros(1,N+1));
    % Elect the heads and decipate their enegies
num_heads = floor(portion/100 * active_nodes_count(end));
if num_heads > 0 %condition to go into the election mode
assigned_heads  = CH_election(C, active_nodes, num_heads, energy(t,1:end-1), locs(1:end-1,:), dists(1:end-1), d0);
heads = cell2mat(assigned_heads(:,1));
for j=1:C
    t = t + 1;
for head_index = 1:length(heads)
    head = heads(head_index);
    head_loc = locs(head);
    head_dist = dists(head);
    % (1) nodes send to heads
    for node_index = cell2mat(assigned_heads(head_index,3))
        node = locs(node_index);
        distance = sqrt(sum((head_loc-node).^2,2));
        if distance >= d0
            energy(t,node_index) = func_tx_energy(energy(t-1,node_index), eta_long, 4, distance);
        else
            energy(t,node_index) = func_tx_energy(energy(t-1,node_index), eta_short, 2, distance);
        end
    end
    % (2) heads receives from nodes
    energy(t, head) = func_rx_energy_normal(energy(t-1, head), assigned_heads{head_index,2});
    % (3) heads sends to sink
    if head_dist >= d0
        energy(t, head) = func_tx_energy_head(energy(t, head), eta_long, 4, head_dist);
    else
        energy(t, head) = func_tx_energy_head(energy(t, head), eta_short, 2, head_dist);
    end
end
     % (4) sink receives from heads
    energy(t,N+1) = func_rx_energy(energy(t-1,N+1), length(heads));

    % Check energies!
    active_nodes = find(energy(t,1:N) > tx_threshold);
    active_nodes_count = cat(2, active_nodes_count, sum(energy(t,1:N) > tx_threshold));
    rx_threshold = 500*8*E_elec; %recieve from @ least 1 head
    if active_nodes_count(end) == 0 || energy(t,N+1) < rx_threshold
        break;
    end
T_end_clusters = t; % save this for later!
end
else %end_if
% Operate in the very normal situation!
t = t + 1;
    for i=1:N
        dist = dists(i);
        if dist <= d0
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_short, 2, dist);
        else
            energy(t,i) = func_tx_energy(energy(t-1,i), eta_long, 4, dist);
        end
    end
    energy(t,N+1) = func_rx_energy_normal(energy(t-1,N+1), active_nodes(end));

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
f5 = figure('Name','Relation between T_1 and C');
stem(Cs, T1s);
grid on
xlabel('C');
ylabel('T_1');
title('Relation between T_1 and C')
%saveas(f5, [pwd '/Figures/Relation between T1 and C'])


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
saveas(f1, [pwd '/Figures/network_topology_Sink_Centered']);
%% Calculate the distances
dists = sqrt(sum((locs-[50 50]).^2,2));
%% Sort Nodes according to Distances
[dists, ids] = sort(dists);
locs = locs(ids,:);
%% remove sink to end of vectors
dists = [dists(2:end); dists(1)];
locs = [locs(2:end,:); locs(1,:)];
%% Initiate the energy
energy = E_initial*ones(1,N+1);
%% GO on cycles!
active_nodes =(N);

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
saveas(f2, [pwd '/Figures/curve of active nodes Sink Centered']);
%% Find the lifetimes T1, T2, T3
T1 = find(active_nodes < 100, 1);
T2 = find(active_nodes < 50, 1);
T3 = find (active_nodes == 0, 1);

if isempty(T3) % in case of the death of the sink before one of the source nodes
   T3 = length(active_nodes); 
end

f3= figure('Name', 'Curve of Active Nodes with lifetimes');
plot(linspace(1,length(active_nodes),length(active_nodes)), active_nodes);
grid on
hold on
plot(T1, active_nodes(T1),'r*');
hold on
plot(T2, active_nodes(T2),'m*');
hold on
plot(T3, active_nodes(T3),'k*');
xlabel('no. of cycle');
ylabel('no. of active nodes');
legend('no. active nodes','T1', 'T2','T3');
title('The number of active nodes at each cycle - Sink Centered');
saveas(f3, [pwd '/Figures/curve of active nodes with lifetimes Sink Centered']);
%% Reamining Energies After T1
active_nodes = N;
energyT1 = E_initial*ones(1,N+1);
for t = 1:T1
    for i=1:N
    dist = dists(i);
        if dist <= d0
            energyT1(i) = func_tx_energy(energyT1(i), eta_short, 2, dist);
        else
            energyT1(i) = func_tx_energy(energyT1(i), eta_long, 4, dist);
        end
    end
    energyT1(N+1) = func_rx_energy(energyT1(N+1), active_nodes(end));
    % Get the no. of active nodes
    tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
    active_nodes = cat(2, active_nodes, sum(energyT1(1:N) >= tx_threshold));
end
f4 = figure('Name','Remaining energies of the N nodes after T1 cycles - Sink Centered');
stem(linspace(1,N,N),energyT1(1:end-1),"filled")
hold on
stem(N+1,energyT1(end),"filled",'LineWidth',1,'Color','r')
hold off
grid on
xlabel('Node Index');
ylabel('Energy (nJ)');
title('Remaining energies of the N nodes after T1 cycles - Sink Centered')
saveas(f4, [pwd '/Figures/Remaining energies of the N nodes after T1 cycles Sink Centered'])
%% Reamining Energies After T2
active_nodes = N;
energyT2 = E_initial*ones(1,N+1);
for t = 1:T2
    for i=1:N
    dist = dists(i);
        if dist <= d0
            energyT2(i) = func_tx_energy(energyT2(i), eta_short, 2, dist);
        else
            energyT2(i) = func_tx_energy(energyT2(i), eta_long, 4, dist);
        end
    end
    energyT2(N+1) = func_rx_energy(energyT2(N+1), active_nodes(end));
    % Get the no. of active nodes
    tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
    active_nodes = cat(2, active_nodes, sum(energyT2(1:N) >= tx_threshold));    
end
f5 = figure('Name','Remaining energies of the N nodes after T2 cycles - Sink Centered');
stem(linspace(1,N,N),energyT2(1:end-1),"filled")
hold on
stem(N+1,energyT2(end),"filled",'LineWidth',1,'Color','r')
hold off
grid on
xlabel('Node Index');
ylabel('Energy (nJ)');
title('Remaining energies of the N nodes after T2 cycles - Sink Centered')
saveas(f5, [pwd '/Figures/Remaining energies of the N nodes after T2 cycles Sink Centered'])
%% Reamining Energies After T3
active_nodes = N;
energyT3 = E_initial*ones(1,N+1);
for t = 1:T3
    for i=1:N
    dist = dists(i);
        if dist <= d0
            energyT3(i) = func_tx_energy(energyT3(i), eta_short, 2, dist);
        else
            energyT3(i) = func_tx_energy(energyT3(i), eta_long, 4, dist);
        end
    end
    energyT3(N+1) = func_rx_energy(energyT3(N+1), active_nodes(end));
    % Get the no. of active nodes
    tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
    active_nodes = cat(2, active_nodes, sum(energyT3(1:N) >= tx_threshold));
end
f6 = figure('Name','Remaining energies of the N nodes after T3 cycles - Sink Centered');
stem(linspace(1,N,N),energyT3(1:end-1),"filled")
hold on
stem(N+1,energyT3(end),"filled",'LineWidth',1,'Color','r')
hold off
grid on
xlabel('Node Index');
ylabel('Energy (nJ)');
title('Remaining energies of the N nodes after T3 cycles - Sink Centered')
saveas(f6, [pwd '/Figures/Remaining energies of the N nodes after T3 cycles Sink Centered'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For Sink in a Corner

locs = cat(1, randperm(area_width), randperm(area_height))';
% assert no node existance @ the network center
for i=1:size(locs,2)
    if locs(i,:) == [0 0]
        locs(i,:) = locs(i,:) + randi([1,20], 1 , 2);
    end
end
% Add the sink node @ the network center
f1 = figure('Name', 'Network Topology Sink Cornered');
scatter(locs(:,1),locs(:,2));
hold on
scatter(0,0,'filled');
locs = cat(1, locs, [0, 0]);
legend('source', 'sink');
title('Network Topology Sink Cornered');
xlabel('distance on x-axis');
ylabel('distance on y-axis');
saveas(f1, [pwd '/Figures/network topology - Sink Cornered']);
%% Calculate the distances
dists = sqrt(sum((locs).^2,2));
%% Sort Nodes according to Distances
[dists, ids] = sort(dists);
locs = locs(ids,:);
%% remove sink to end of vectors
dists = [dists(2:end); dists(1)];
locs = [locs(2:end,:); locs(1,:)];
%% Initiate the energy
energy = E_initial*ones(1,N+1);
%% GO on cycles!
active_nodes =(N);

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
title('The number of active nodes at each cycle - Sink Cornered');
saveas(f2, [pwd '/Figures/curve of active nodes Sink Cornered']);
%% Find the lifetimes T1, T2, T3
T1 = find(active_nodes < 100, 1);
T2 = find(active_nodes < 50, 1);
T3 = find (active_nodes == 0, 1);

if isempty(T3) % in case of the death of the sink before one of the source nodes
   T3 = length(active_nodes); 
end

f3= figure('Name', 'Curve of Active Nodes with lifetimes');
plot(linspace(1,length(active_nodes),length(active_nodes)), active_nodes);
grid on
hold on
plot(T1, active_nodes(T1),'r*');
hold on
plot(T2, active_nodes(T2),'m*');
hold on
plot(T3, active_nodes(T3),'k*');
xlabel('no. of cycle');
ylabel('no. of active nodes');
legend('no. active nodes','T1', 'T2','T3');
title('The number of active nodes at each cycle - Sink Cornered');
saveas(f3, [pwd '/Figures/curve of active nodes with lifetimes Sink Cornered']);
%% Reamining Energies After T1
active_nodes = N;
energyT1 = E_initial*ones(1,N+1);
for t = 1:T1
    for i=1:N
    dist = dists(i);
        if dist <= d0
            energyT1(i) = func_tx_energy(energyT1(i), eta_short, 2, dist);
        else
            energyT1(i) = func_tx_energy(energyT1(i), eta_long, 4, dist);
        end
    end
    energyT1(N+1) = func_rx_energy(energyT1(N+1), active_nodes(end));
    % Get the no. of active nodes
    tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
    active_nodes = cat(2, active_nodes, sum(energyT1(1:N) >= tx_threshold));
end
f4 = figure('Name','Remaining energies of the N nodes after T1 cycles - Sink Cornered');
stem(linspace(1,N,N),energyT1(1:end-1),"filled")
hold on
stem(N+1,energyT1(end),"filled",'LineWidth',1,'Color','r')
hold off
grid on
xlabel('Node Index');
ylabel('Energy (nJ)');
title('Remaining energies of the N nodes after T1 cycles - Sink Cornered')
saveas(f4, [pwd '/Figures/Remaining energies of the N nodes after T1 cycles Sink Cornered'])
%% Reamining Energies After T2
active_nodes = N;
energyT2 = E_initial*ones(1,N+1);
for t = 1:T2
    for i=1:N
    dist = dists(i);
        if dist <= d0
            energyT2(i) = func_tx_energy(energyT2(i), eta_short, 2, dist);
        else
            energyT2(i) = func_tx_energy(energyT2(i), eta_long, 4, dist);
        end
    end
    energyT2(N+1) = func_rx_energy(energyT2(N+1), active_nodes(end));
    % Get the no. of active nodes
    tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
    active_nodes = cat(2, active_nodes, sum(energyT2(1:N) >= tx_threshold));
end
f5 = figure('Name','Remaining energies of the N nodes after T2 cycles - Sink Cornered');
stem(linspace(1,N,N),energyT2(1:end-1),"filled")
hold on
stem(N+1,energyT2(end),"filled",'LineWidth',1,'Color','r')
hold off
grid on
xlabel('Node Index');
ylabel('Energy (nJ)');
title('Remaining energies of the N nodes after T2 cycles - Sink Cornered')
saveas(f5, [pwd '/Figures/Remaining energies of the N nodes after T2 cycles Sink Cornered'])
%% Reamining Energies After T3
active_nodes = N;
energyT3 = E_initial*ones(1,N+1);
for t = 1:T3
    for i=1:N
    dist = dists(i);
        if dist <= d0
            energyT3(i) = func_tx_energy(energyT3(i), eta_short, 2, dist);
        else
            energyT3(i) = func_tx_energy(energyT3(i), eta_long, 4, dist);
        end
    end
    energyT3(N+1) = func_rx_energy(energyT3(N+1), active_nodes(end));
    % Get the no. of active nodes
    tx_threshold = k*(E_elec+E_agg) + k*eta_short*1^2;
    active_nodes = cat(2, active_nodes, sum(energyT3(1:N) >= tx_threshold));
end
f6 = figure('Name','Remaining energies of the N nodes after T3 cycles - Sink Cornered');
stem(linspace(1,N,N),energyT3(1:end-1),"filled")
hold on
stem(N+1,energyT3(end),"filled",'LineWidth',1,'Color','r')
hold off
grid on
xlabel('Node Index');
ylabel('Energy (nJ)');
title('Remaining energies of the N nodes after T3 cycles - Sink Cornered')
saveas(f6, [pwd '/Figures/Remaining energies of the N nodes after T3 cycles Sink Cornered'])
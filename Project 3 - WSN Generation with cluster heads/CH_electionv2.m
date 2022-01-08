function assigned_heads = CH_electionv2(C, active_nodes, num_heads, energies, locations, distances, do)
% This function perform the election of cluster heads.
% @ Inputs:
    % C: The number of cycles for the heads
    % active_nodes: A vector of the active nodes
    % num_heads: The number of heads to be elected    
    % energies: The remaining energies for all nodes
    % locations: The location of each node
    % distances: The distances between each node & the sink
    % do: The threshold for determining the path loss exponent @
    %     transmission
% @ Outputs:
    % assigned_heads: A cell contains 3 columns (the elected heads, number
    %                 of nodes assigned to each head, a vector of the assigned 
    %                 nodes to each head)

% ---------------------------------------------------------------------%
% Constants & Parameters
E_elec = 50; E_agg = 50; 
eta_short = 10; eta_long = 0.0013;    %nJ/bit/m2;
k = 625*8; %number_of_bits per cycle
N = length(active_nodes); %number of current active nodes

% Defining the rx and tx functions
func_rx_energy = @(E, n) max(0, E - C*(n*k*E_elec));
func_tx_energy = @(E,eta,exp,d) max(0, E - C*(k*(E_elec + E_agg) + k*eta*d.^exp));
% ---------------------------------------------------------------------%
% Create Ranges'
% ranges  = cell(Sector Number, sector angle range, Assigned Nodes)
assigned_heads = cell(num_heads, 3);
for i = 1:  num_heads
    assigned_heads{i,2} = [(i-1)*(2*pi/num_heads) i*2*pi/num_heads];
end

% Detect nodes in each range
for i = 1:N
    temp = locations(i,:) - [50 50];
    angle = atan2(temp(2),temp(1));
    if(angle < 0)
        angle = angle + 2*pi;
    end
    for j  = 1:num_heads
        range = assigned_heads{j,2};
        if (angle >= range(1) && angle < range(2))
            assigned_heads{j,3} = [assigned_heads{j,3}, i];
            break;
        end
    end
end

% Elect Head for Each cluster
for i = 1:num_heads
    if ~isempty(assigned_heads{i,3})
        nodes = assigned_heads{i,3}; % nodes in the cluster
        clusterDists = distances(assigned_heads{i,3});
        clusterEnergies = energies(assigned_heads{i,3});
         % energies after recpetion
        result_energies = func_rx_energy(clusterEnergies,length(nodes));
        % energies after transmission
        for j = 1:length(nodes)
            head = nodes(j);
            if distances(head) >= do && ~isempty(nodes)
                result_energies(i) = func_tx_energy(result_energies(j),eta_long,4,distances(head));
            elseif assigned_heads{i, 2} ~= 0
                result_energies(i) = func_tx_energy(result_energies(j),eta_short,2,distances(head));
            end
        end
        nodes = nodes(result_energies > 0); %nodes that satisfy the energy constraint
        if ~isempty(nodes)
%             headsDists = clusterDists(result_energies > 0); %dists of potential heads
%             m = mean(headsDists); %mean of the distances of those nodes
%             [~, headIndex] = min(abs(headsDists-m));
%             if length(headIndex)>1
%                 headIndex = headIndex(1);
%             end

           headIndex = randi(length(nodes));
            assigned_heads{i,1} = nodes(headIndex);
            assigned_heads{i,3} = assigned_heads{i,3}(assigned_heads{i,3} ~= nodes(headIndex));
            assigned_heads{i,2} = length(assigned_heads{i,3});
        end
    end
end
end
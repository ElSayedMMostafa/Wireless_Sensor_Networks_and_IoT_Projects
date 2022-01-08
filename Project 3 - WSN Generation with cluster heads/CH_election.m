function assigned_heads = CH_election(C, active_nodes, num_heads, energies, locations, distances, do)
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

 % loop over all the possibilities (it should be N^C instead but we changed it for easier simulation)
while 1
    % Get a random sample
    heads = randsample(active_nodes, num_heads);
    heads_locations = locations(heads, :);
    heads_distances = distances(heads);
    %if any(heads_distances >  mean(distances(active_nodes))) 
    if any(heads_distances >  35) 
        continue
    end
    % assigned_heads = cell[heads, no. assigned nodes,vector of assigned nodes]
    assigned_heads = cell(num_heads, 3);
    assigned_heads(:,1) = num2cell(heads');
    assigned_heads(:,2) = num2cell(0);
    % loop over the remaining nodes to assign nodes to cluster heads
    for index = 1:N
        if any(index == heads)
            continue
        end
        node_location = locations(index,:);
        [~, selected_head] = min(sqrt(sum((heads_locations-node_location).^2,2)));
        % head_row = find(assigned_heads(:,1) == selected_head);
        assigned_heads{selected_head, 2} = assigned_heads{selected_head, 2} + 1;
        assigned_heads{selected_head, 3} = cat(2,assigned_heads{selected_head, 3}, index);
    end
    % ------- after assignment processing ---------- %
    % Refuse the case where a head is not assigned to any nodes
    if any(cell2mat(assigned_heads(:,2)) == 0) %&& r ~= N*C
        continue
    end
    % energies after recpetion
    result_energies = func_rx_energy(energies(heads),cell2mat(assigned_heads(:,2))');
    % energies after transmission
    for i = 1:length(heads)
        head = heads(i);
        if distances(head) >= do && assigned_heads{i, 2} ~= 0
            result_energies(i) = func_tx_energy(result_energies(i),eta_long,4,distances(head));
        elseif assigned_heads{i, 2} ~= 0
            result_energies(i) = func_tx_energy(result_energies(i),eta_short,2,distances(head));
        end
    end
    % Check the remaining energies
    if all(result_energies > 0)
        return
    end
end

end
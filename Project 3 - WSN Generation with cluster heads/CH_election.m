function heads = CH_election(C, portion, energies, locations, distances, do)
% This function perform the election of cluster heads.
% @ Inputs:
    % C: The number of cycles for the heads
    % portion: The portion of nodes to be elected as heads (in percentage)
    % energies: The remaining energies for all nodes
    % locations: The location of each node
    % distances: The distances between each node & the sink
    % do: The threshold for determining the path loss exponent @
    %     transmission
% @ Outputs:
    % selected_heads: Indecies of the elected heads
% ---------------------------------------------------------------------%
% Constants & Parameters
E_elec = 50; E_agg = 50; 
eta_short = 10; eta_long = 0.0013;    %nJ/bit/m2;
N = length(energies);
num_heads = floor(portion/100 * N);
% Defining the rx and tx functions
func_rx_energy = @(E, n) max(0, E - C*(n*500*8*E_elec));
func_tx_energy = @(E,eta,exp,d) max(0, E - C*((500*8)*(E_elec + E_agg) + (500*8)*eta*d.^exp));

% Get all the possible combinations & shuffle them
all_possibilities = nchoosek(1:N, num_heads);
all_possibilities = all_possibilities(randperm(size(all_possibilities, 1)), :);
for r = 1:size(all_possibilities, 1) % loop over all the possibilities
    heads = all_possibilities(r,:);
    heads_locations = locations(heads, :);
    assigned_heads = zeros(num_heads,2); assigned_heads(:,1) = heads;
    
    % loop over the remaining nodes
    for index = 1:N
        if any(index == heads)
            continue
        end
        node_location = locations(index);
        [~, selected_head] = min(sqrt(sum((heads_locations-node_location).^2,2)));
        % head_row = find(assigned_heads(:,1) == selected_head);
        assigned_heads(selected_head, 2) = assigned_heads(selected_head, 2) + 1;
    end    
    % energies after recpetion
    temp_energies = func_rx_energy(energies(heads),assigned_heads(:,2));
    % energies after transmission
    for i = length(heads)
        head = heads(i);
        if distances(head) >= do
            temp_energies(i) = func_tx_energy(temp_energies(i),eta_long,4,distances(head));
        else
            temp_energies(i) = func_tx_energy(temp_energies(i),eta_short,2,distances(head));
        end
    end

    % Check the remaining energies
    if all(temp_energies > 0)
        return
    end
end

end
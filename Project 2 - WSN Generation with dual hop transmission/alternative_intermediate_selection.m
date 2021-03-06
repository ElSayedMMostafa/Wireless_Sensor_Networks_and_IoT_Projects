%function [min_Index, dist_to_node] = alternative_intermediate_selection(node_index, locations, dists_to_sink)
%function [max_Index, dist_to_node] = alternative_intermediate_selection(node_index, locations, remaining_enrgies)
function [max_Index, dist_to_node] = alternative_intermediate_selection(node_index, locations, dists_to_sink, remaining_enrgies)

% This function selects the best intermediate node for dual hop
% transmission given a specific transmitter node assuming that the sink in
% located in [50, 225]

% The choice of the intermediate node here is dependant on the
% remaining enrgies for the possible intermediate nodes and the minimum
% sum of the distance between the tx node and the intermediate node and the
% intermediate node and the sink

% @ Inputs:
    % node_index: The index of tx node location
    % locations: locations of the all nodes (except for the sink)
    % dists_to_sink: 1D vector contains distances from nodes to sink
% @ Outputs:
    % max_Index: The index of the chosen node
    % dist_to_node: The distance between the tx node and its intermediate
node = locations(node_index, :);
locations(node_index, :) = [inf, inf];
dists_to_sink(node_index) = inf;
remaining_enrgies(node_index) = -inf;

dists_to_node = sqrt(sum((locations-node).^2,2));
dists_sum = dists_to_sink + dists_to_node;

min_Index = zeros(1,5);
for i = 1:5
    [~, min_Index(i)] = min(dists_sum);
    dists_sum(min_Index(i)) = inf;
end
[~, max_Index] = max(remaining_enrgies(min_Index)); 
max_Index = min_Index(max_Index);
dist_to_node = dists_to_node(max_Index);
end
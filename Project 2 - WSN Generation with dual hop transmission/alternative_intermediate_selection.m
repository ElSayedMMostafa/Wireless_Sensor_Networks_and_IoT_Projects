%function [min_Index, dist_to_node] = alternative_intermediate_selection(node_index, locations, dists_to_sink)
function [max_Index, dist_to_node] = alternative_intermediate_selection(node_index, locations, remaining_enrgies)

% This function selects the best intermediate node for dual hop
% transmission given a specific transmitter node assuming that the sink in
% located in [50, 225]

% The choice of the intermediate node here is only dependant on the
% remaining enrgies for the possible intermediate nodes.

% @ Inputs:
    % node_index: The index of tx node location
    % locations: locations of the all nodes (except for the sink)
    % dists_to_sink: 1D vector contains distances from nodes to sink
% @ Outputs:
    % max_Index: The index of the chosen node
    % dist_to_node: The distance between the tx node and its intermediate
node = locations(node_index, :);
locations(node_index, :) = [inf, inf];
remaining_enrgies(node_index) = -inf;


dists_to_node = sqrt(sum((locations-node).^2,2));


[~, max_Index] = max(remaining_enrgies); 
dist_to_node = dists_to_node(max_Index);
end
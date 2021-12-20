function [min_Index, dist_to_node] = alternative_intermediate_selection(d0, node_index, locations, dists_to_sink)
% This function selects the best intermediate node for dual hop
% transmission given a specific transmitter node assuming that the sink in
% located in [50, 50]
% @ Inputs:
    % d0: The threshold between sending using path loss of 2 and 4
    % node_index: The index of tx node location
    % locations: locations of the all nodes (except for the sink)
    % dists_to_sink: 1D vector contains distances from nodes to sink
% @ Outputs:
    % min_Index: The index of the chosen node
    % dist_to_node: The distance between the tx node and its intermediate
node = locations(node_index, :);
locations(node_index, :) = [inf, inf];
dists_to_sink(node_index) = inf;


dists_to_node = sqrt(sum((locations-node).^2,2));
dists_sum = dists_to_sink + dists_to_node;
dists_sum(dists_to_sink > d0) = inf;

[~, min_Index] = min(dists_sum); 
dist_to_node = dists_to_node(min_Index);

end
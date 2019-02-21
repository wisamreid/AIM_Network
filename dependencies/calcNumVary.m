function nVary = calcNumVary(vary)

nVary = cellfun(@length, vary(:,3), 'uni',0);
nVary = prod([nVary{:}]);

fprintf('nVary = %i \n', nVary);
end
import networkx as nx

def bfs_endpoint(graph: nx.Graph, compound: str, element: str) -> dict[str, set]:
    # remove every node that is not the element (if data is not clean, we might have wrong symmetry or direction edges here)
    cleaned_graph = graph.copy()
    for key, data in graph.nodes(data=True):
        if data["element"] != element:
            cleaned_graph.remove(key)
        
    starting_nodes = []
    # holds subgraphs for every compound, only SYMMETRY edges are present
    compound_dict: dict[str, nx.Graph] = {}
    # holds list of nodes for every compound
    compound_node_dict: dict[str, set[any]] = {}
    
    # for every node in the given graph
    for [key, data] in cleaned_graph.nodes(data=True):
        # if the current compound is already present, add the node to the compound list
        if data['compound_name'] in compound_node_dict:
            compound_node_dict[data['compound_name']].add(key)
        # otherwise create the compound with the given node in the compound list
        else:
            compound_node_dict[data['compound_name']] = set([key])
        # if the current compound is the starting compound and current element is the starting element
        # add them to the list of starting nodes
        if data['compound_name'] == compound:
            starting_nodes.append(key)
    
    # for every set of nodes of each compound
    for [key, nodes] in compound_node_dict.items():
        # create subgraph only containing nodes of this compound
        compound_graph: nx.Graph = cleaned_graph.subgraph(nodes)
        # create symmetry graph, that only contains edges labed as "SYMMETRY"
        symmetry_graph: nx.Graph = cleaned_graph.copy()
        # remove other edges
        for u,v,edge_type in compound_graph.edges(data="transition"):
            if edge_type != "TransitionType.SYMMETRY":
                symmetry_graph.remove_edge(u,v)
        # save the symmetry graph
        compound_dict[key] = symmetry_graph
    
    endpoints: dict[str, set[any]] = {}
    
    # for every starting node
    for start_node in starting_nodes:
        # save a dictionary of potential endpoint compounds and their potential endpoint atoms
        current_endpoints: dict[str, list[any]] = {}
        # normal bfs from here on
        visited = set()
        queue = [start_node]

        while queue:
            next_node = queue.pop(0)
            if next_node in visited:
                continue
            visited.add(next_node)
            for neighbor in cleaned_graph.neighbors(next_node):
                if neighbor not in visited and neighbor not in queue:
                    queue.append(neighbor)
                    
                    # the source atom is always already in the current_endpoints list
                    # except for the starting compound. We assume the starting compound is not an endpoint
                    compound_target = cleaned_graph.nodes[neighbor]['compound_name']
                    compound_source = cleaned_graph.nodes[next_node]['compound_name']
                    # if we traversed to a new compound
                    if compound_source != compound_target:
                        # if the new compound is in our current endpoints
                        if compound_target in current_endpoints:
                            # retrieve all recorded potential endpoint atoms for this compound
                            nodes = current_endpoints[compound_target]
                            # retrieve the current subgraph for this compound
                            sym_graph = compound_dict[compound_target]
                            to_be_removed = None
                            # for every node we labeled as "potential endpoint atom"
                            for node in nodes:
                                try:
                                    # find a shortest path between the target and the nodes
                                    nx.shortest_path(sym_graph, neighbor, node)
                                    # if a path is found, this symmetry graph has two connections, the current potential endpoint is not an endpoint anymore
                                    # there can never be more then one shortest path, since all potential endpoint atoms are disconnected between themselfs
                                    to_be_removed = node
                                    break
                                except nx.NetworkXNoPath:
                                    pass
                            # if we haven't found any potential endpoint, that is connected, this current target atom is a new potential endpoint atom
                            if to_be_removed is None:
                                nodes.append(neighbor)
                            # if we found a potential endpoint, this symmetry subgraph has two connections/no endpoint, the potential endpoint needs to be removed
                            else:
                                nodes.remove(to_be_removed)
                        else:
                            current_endpoints[compound_target] = [neighbor]
            for [key, nodes] in current_endpoints.items():
                if len(nodes) == 0:
                    continue
                if key not in endpoints:
                    endpoints[key] = set(nodes)
                else:
                    for node in nodes:
                        endpoints[key].add(node)
        return endpoints
                            
                            

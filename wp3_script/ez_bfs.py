import networkx as nx

def bfs_endpoint(graph: nx.Graph, molecule, element):
    start_points = []
    molecule_dict: dict[str, nx.Graph] = dict()
    molecule_dict_node: dict[str, set[any]] = dict()
    for node in graph.nodes(data=True):
        if node[1]['compound_name'] in  molecule_dict_node:
            molecule_dict_node[node[1]['compound_name']].add(node[0])
        else:
            molecule_dict_node[node[1]['compound_name']] = set([node[0]])
        if node[1]['compound_name'] == molecule and node[1]['element'] == element:
            start_points.append(node[0])
    
    for node_set in molecule_dict_node.items():
        compound_graph: nx.Graph = graph.subgraph(node_set[1])
        symmetry_graph: nx.Graph = graph.copy()
        for u,v,edge_type in compound_graph.edges(data="transition"):
            if edge_type != "TransitionType.SYMMETRY":
                symmetry_graph.remove_edge(u,v)
        molecule_dict[node[0]] = symmetry_graph
    
    for start_point in start_points:
        end_points: dict[str, list[any]] = dict[]
        visited = set()
        queue = [start_point]

        while queue:
            next = queue.pop(0)
            if next in visited:
                continue
            visited.add(next)
            for neighbor in graph.neighbors(next):
                if neighbor not in visited and neighbor not in queue:
                    compound_b = graph.nodes[neighbor]['compound_name']
                    compound_a = graph.nodes[next]['compound_name']
                    if compound_a != compound_b:
                        if compound_a in end_points:
                            nodes = end_points[compound_a]
                            graph = molecule_dict[compound_a]
                            is_deleted = False
                            for node in nodes:
                                try:
                                    nx.shortest_path(graph, compound_a, node)
                                    del end_points[compound_a]
                                    is_deleted = True
                                    break
                                except nx.NetworkXNoPath:
                                    pass
                            if is_deleted is not False:
                                
                            
                            

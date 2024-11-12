def topological_sort(graph):
    # Initialize states for each node
    state = {v: "unvisited" for v in graph}
    has_cycle = [False]  # Use a list to keep track of the cycle status (mutable inside DFS)
    result = []

    # DFS function
    def dfs(v):
        if state[v] == "visiting":
            # Cycle detected
            has_cycle[0] = True
            return
        if state[v] == "visited":
            return
        
        # Mark as visiting
        state[v] = "visiting"
        
        # Visit all successors
        for neighbor in graph[v]:
            dfs(neighbor)
        
        # Mark as visited
        state[v] = "visited"
        
        # Add to result list
        result.append(v)

    # Run DFS for each unvisited node
    for vertex in graph:
        if state[vertex] == "unvisited":
            dfs(vertex)
        
        # If cycle is detected, break early
        if has_cycle[0]:
            break
    
    # Check if a cycle was detected
    if has_cycle[0]:
        return "No topological order, cycle detected"
    else:
        return result[::-1]  # Reverse the result list for topological ordering


# Example usage:
# Define a graph as an adjacency list
graph = {
    0: [1, 2],
    1: [3],
    2: [3],
    3: [4],
    4: []
}

result = topological_sort(graph)
print(result)

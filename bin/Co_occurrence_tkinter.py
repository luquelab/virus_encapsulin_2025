import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import sys

# Global variables to hold data and domain assignments.
loaded_df = None         # The full CSV read as DataFrame.
extracted_header = []    # List of domain names extracted from row 1.
cleaned_data = None      # The numeric data extracted starting from the data start row.
domain_assignment_vars = {}  # Dictionary for group assignments.
domain_discard_vars = {}     # Dictionary for discard checkboxes.
group_assignments = {}       # Final domain group assignments.
full_graph = None            # The full co-occurrence graph (G).

# Default file paths.
default_input_path = "/Users/abelardoaguilar/projects/github_repos/mini-devel/mini-devel/results/verse_filesystem/Effect_of_DTRs_verse/data/Tree_Annotation/minimal_plus_no_DTR.csv"
default_output_path = "domain_cooccurrence_graph.svg"

###############################################################################
# Phase 1: Load File and Set Parameters
###############################################################################
def load_file():
    global loaded_df
    try:
        input_path = input_file_var.get()
        loaded_df = pd.read_csv(input_path, header=None)
        print("File preview (first 10 rows):")
        print(loaded_df.head(10))
        messagebox.showinfo("Info", "File loaded successfully. Now click 'Assign Domain Groups' to proceed.")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load file: {e}")

###############################################################################
# Phase 2: Assign Domain Groups (with Discard, Uncategorized Options)
###############################################################################
def assign_groups():
    global extracted_header, cleaned_data, domain_assignment_vars, domain_discard_vars
    if loaded_df is None:
        messagebox.showerror("Error", "Load a file first.")
        return

    try:
        # Get parameters.
        skip_rows = [int(x.strip()) for x in skip_rows_var.get().split(',') if x.strip() != '']
        data_start_row = int(data_start_row_var.get())
        
        # Extract the header from row 1 (Field labels) using columns 2 up to second-to-last.
        extracted_header = loaded_df.iloc[1, 2:-1].tolist()
        print("\nExtracted domain names from 'Field labels':")
        print(extracted_header)
        
        # Extract the data rows (starting at data_start_row) from columns 2:-1.
        cleaned_data = loaded_df.iloc[data_start_row:, 2:-1].copy()
        cleaned_data.columns = extracted_header  # Set domain names as header.
        
        # Convert all data to numeric.
        cleaned_data = cleaned_data.apply(pd.to_numeric, errors='coerce')
        # Replace -1 with 0.
        cleaned_data = cleaned_data.replace(-1, 0)
        # Keep only numeric columns.
        cleaned_data = cleaned_data.select_dtypes(include=['number'])
        
        print("\nCleaned Data Shape:", cleaned_data.shape)
        print(cleaned_data.head())
        
        # Open a new window for domain assignment.
        domain_window = tk.Toplevel(root)
        domain_window.title("Assign Domain Groups")
        tk.Label(domain_window, text="For each domain, choose its group (Encapsulin, MCP, or Uncategorized) \nand whether to discard it:").grid(row=0, column=0, columnspan=4, pady=5)
        
        # Initialize dictionaries for assignments and discard options.
        domain_assignment_vars = {}
        domain_discard_vars = {}
        for idx, domain in enumerate(extracted_header):
            tk.Label(domain_window, text=domain).grid(row=idx+1, column=0, sticky="w")
            group_var = tk.StringVar(value="MCP")  # Default group is MCP.
            domain_assignment_vars[domain] = group_var
            # Three radio buttons.
            tk.Radiobutton(domain_window, text="Encapsulin", variable=group_var, value="Encapsulin").grid(row=idx+1, column=1)
            tk.Radiobutton(domain_window, text="MCP", variable=group_var, value="MCP").grid(row=idx+1, column=2)
            tk.Radiobutton(domain_window, text="Uncategorized", variable=group_var, value="Uncategorized").grid(row=idx+1, column=3)
            # Checkbox to discard the domain.
            discard_var = tk.BooleanVar(value=False)
            domain_discard_vars[domain] = discard_var
            tk.Checkbutton(domain_window, text="Discard", variable=discard_var).grid(row=idx+1, column=4, padx=10)
        
        tk.Button(domain_window, text="Run Analysis", command=run_analysis).grid(row=len(extracted_header)+1, column=1, columnspan=2, pady=10)
    except Exception as e:
        messagebox.showerror("Error", f"Error during group assignment: {e}")

###############################################################################
# Phase 3: Run Analysis & Build Graph; then Compute Evolutionary Paths
###############################################################################
def run_analysis():
    global group_assignments, full_graph
    try:
        # Drop discarded domains.
        domains_to_keep = [d for d in extracted_header if not domain_discard_vars[d].get()]
        if not domains_to_keep:
            messagebox.showerror("Error", "No domains selected (all discarded).")
            return
        
        # Update cleaned_data to keep only non-discarded domains.
        data_subset = cleaned_data[domains_to_keep]
        
        # Compute the co-occurrence matrix.
        co_occurrence_matrix = data_subset.T.dot(data_subset)
        print("\nCo-occurrence Matrix Shape:", co_occurrence_matrix.shape)
        print("Co-occurrence Matrix Columns:", list(co_occurrence_matrix.columns))
        print(co_occurrence_matrix.head())
        
        # Use updated domains.
        domains = list(co_occurrence_matrix.columns)
        
        # Get user-assigned groups.
        group_assignments = {d: domain_assignment_vars[d].get() for d in extracted_header if d in domains_to_keep}
        print("\nUser-assigned domain groups (for kept domains):")
        print(group_assignments)
        
        # Build the full graph G.
        G = nx.Graph()
        for domain in domains:
            G.add_node(domain, type="domain")
        for domain1 in domains:
            for domain2 in domains:
                if domain1 != domain2:
                    value = co_occurrence_matrix.loc[domain1, domain2]
                    if value > 0:
                        for i in range(int(value)):
                            protein_node = f"Protein_{domain1}_{domain2}_{i}"
                            cat1 = group_assignments.get(domain1)
                            cat2 = group_assignments.get(domain2)
                            if cat1 == "Uncategorized" or cat2 == "Uncategorized":
                                node_color = "gray"
                            elif cat1 != cat2:
                                node_color = "red"
                            else:
                                node_color = "green"
                            G.add_node(protein_node, type="protein", color=node_color)
                            G.add_edge(protein_node, domain1)
                            G.add_edge(protein_node, domain2)
        
        # Remove isolated domain nodes.
        for domain in domains:
            if G.degree(domain) == 0:
                G.remove_node(domain)
        
        # Save the full graph version (with Uncategorized).
        pos = nx.spring_layout(G, seed=42)
        domain_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "domain"]
        red_protein_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "red"]
        green_protein_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "green"]
        gray_protein_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "gray"]
        node_sizes = [max(100, min(500, G.degree(n)*100)) for n in domain_nodes]
        
        plt.figure()
        nx.draw_networkx_nodes(G, pos, nodelist=domain_nodes, node_size=node_sizes, node_color="lightblue", label="Domains")
        nx.draw_networkx_nodes(G, pos, nodelist=red_protein_nodes, node_size=2, node_color="red", label="Non-canonical")
        nx.draw_networkx_nodes(G, pos, nodelist=green_protein_nodes, node_size=1, node_color="lightgreen", label="Canonical")
        nx.draw_networkx_nodes(G, pos, nodelist=gray_protein_nodes, node_size=2, node_color="gray", label="Uncategorized")
        nx.draw_networkx_edges(G, pos, style="solid", width=0.2)
        nx.draw_networkx_labels(G, pos, labels={n: n for n in domain_nodes}, font_size=1, font_color="black", font_weight="bold")
        output_with = output_file_var.get()
        plt.savefig(output_with, format="svg", bbox_inches="tight")
        plt.close()
        
        # Store full graph.
        full_graph = G.copy()
        
        # Compute evolutionary paths.
        compute_evolutionary_paths()
        
    except Exception as e:
        messagebox.showerror("Error", str(e))

###############################################################################
# Phase 4: Compute Evolutionary Paths (Hamiltonian/Spanning Tree Approach)
###############################################################################
def compute_evolutionary_paths():
    """
    Build a weighted, undirected domain graph H from the co-occurrence matrix.
    Each edge (u,v) exists if u and v co-occur.
    Weight is defined as:
      - 3 if either domain is Uncategorized,
      - 5 if domains belong to different groups,
      - 1 if both domains are in the same group.
    Then, if the number of domains is small (≤12), enumerate all spanning trees.
    Export the best spanning tree plus up to five additional alternative evolutionary paths.
    If there are more than 12 domains, use MST as an approximation.
    """
    global full_graph, group_assignments
    domain_nodes = [n for n, d in full_graph.nodes(data=True) if d.get("type") == "domain"]
    H = nx.Graph()
    H.add_nodes_from(domain_nodes)
    
    def edge_cost(u, v):
        cat_u = group_assignments.get(u)
        cat_v = group_assignments.get(v)
        if cat_u == "Uncategorized" or cat_v == "Uncategorized":
            return 3
        elif cat_u != cat_v:
            return 5
        else:
            return 1
    
    co_occurrence_matrix = cleaned_data.T.dot(cleaned_data)
    for i, u in enumerate(domain_nodes):
        for v in domain_nodes[i+1:]:
            if co_occurrence_matrix.loc[u, v] > 0:
                cost = edge_cost(u, v)
                H.add_edge(u, v, weight=cost)
    
    if not nx.is_connected(H):
        messagebox.showerror("Error", "The domain graph is not connected; cannot compute spanning trees.")
        return
    
    candidate_trees = []
    if len(domain_nodes) > 12:
        messagebox.showinfo("Info", "Too many domains to enumerate all spanning trees. Using MST as approximation.")
        best_tree = nx.minimum_spanning_tree(H, weight='weight')
        best_cost = sum(data["weight"] for u, v, data in best_tree.edges(data=True))
        candidate_trees = [(best_tree, best_cost)]
    else:
        candidate_trees = list(enumerate_spanning_trees(H))
        if not candidate_trees:
            messagebox.showerror("Error", "No spanning trees found.")
            return
        candidate_trees = [(T, sum(H[u][v]["weight"] for u, v in T)) for T in candidate_trees]
        candidate_trees.sort(key=lambda x: x[1])
    
    # Export up to 6 alternative evolutionary paths.
    num_paths = min(6, len(candidate_trees))
    for i in range(num_paths):
        tree_edges, cost = candidate_trees[i]
        T = nx.Graph()
        T.add_nodes_from(domain_nodes)
        T.add_edges_from(tree_edges)
        pos_tree = nx.spring_layout(T, seed=42)
        plt.figure()
        nx.draw_networkx_nodes(T, pos_tree, node_size=300, node_color="lightblue")
        nx.draw_networkx_edges(T, pos_tree, edge_color="black", width=2)
        nx.draw_networkx_labels(T, pos_tree, font_size=8)
        output_path_evol = output_file_var.get().replace(".svg", f"_evol_{i+1}.svg")
        plt.title(f"Evolutionary Path {i+1} (Cost: {cost})")
        plt.savefig(output_path_evol, format="svg", bbox_inches="tight")
        plt.close()
        print(f"Evolutionary path {i+1} saved as: {output_path_evol}")

###############################################################################
# Helper: Enumerate All Spanning Trees of H
###############################################################################
def enumerate_spanning_trees(H):
    """
    Enumerate all spanning trees of an undirected connected graph H using backtracking.
    Returns a list of spanning trees, each represented as a set of edges (sorted tuples).
    Suitable only for small graphs.
    """
    trees = []
    nodes = list(H.nodes())
    if not nodes:
        return trees
    root = nodes[0]
    def backtrack(current_tree, used):
        if len(used) == len(nodes):
            trees.append(set(current_tree))
            return
        for u in list(used):
            for v in H.neighbors(u):
                if v not in used:
                    edge = tuple(sorted((u, v)))
                    current_tree.append(edge)
                    used.add(v)
                    backtrack(current_tree, used)
                    used.remove(v)
                    current_tree.pop()
    backtrack([], set([root]))
    return trees

###############################################################################
# GUI Setup
###############################################################################
root = tk.Tk()
root.title("Co-occurrence Graph Generator")

# Default parameters.
input_file_var = tk.StringVar(value=default_input_path)
output_file_var = tk.StringVar(value=default_output_path)
skip_rows_var = tk.StringVar(value="0,2,3,4")
data_start_row_var = tk.StringVar(value="5")

tk.Label(root, text="Input CSV File:").grid(row=0, column=0, sticky="e")
tk.Entry(root, textvariable=input_file_var, width=80).grid(row=0, column=1)
tk.Button(root, text="Browse", command=lambda: input_file_var.set(filedialog.askopenfilename(title="Select Input CSV File", filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]))).grid(row=0, column=2)

tk.Label(root, text="Output SVG File:").grid(row=1, column=0, sticky="e")
tk.Entry(root, textvariable=output_file_var, width=80).grid(row=1, column=1)
tk.Button(root, text="Browse", command=lambda: output_file_var.set(filedialog.asksaveasfilename(title="Save Output SVG File", defaultextension=".svg", filetypes=[("SVG Files", "*.svg"), ("All Files", "*.*")]))).grid(row=1, column=2)

tk.Label(root, text="Skip Rows (comma-separated):").grid(row=2, column=0, sticky="e")
tk.Entry(root, textvariable=skip_rows_var, width=40).grid(row=2, column=1, sticky="w")

tk.Label(root, text="Data Start Row:").grid(row=3, column=0, sticky="e")
tk.Entry(root, textvariable=data_start_row_var, width=10).grid(row=3, column=1, sticky="w")

tk.Button(root, text="Load File", command=load_file).grid(row=4, column=1, pady=5)
tk.Button(root, text="Assign Domain Groups", command=assign_groups).grid(row=5, column=1, pady=5)

root.mainloop()

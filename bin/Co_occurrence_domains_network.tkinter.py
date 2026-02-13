import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import networkx as nx
import itertools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ============================================================
# Helpers: Scrollable frame
# ============================================================
class ScrollableFrame(ttk.Frame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.canvas = tk.Canvas(self, highlightthickness=0)
        self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)

        self.inner.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        self.canvas_window = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")

        # Make inner frame resize with canvas width
        self.canvas.bind("<Configure>", self._on_canvas_configure)

    def _on_canvas_configure(self, event):
        self.canvas.itemconfigure(self.canvas_window, width=event.width)


# ============================================================
# Main App
# ============================================================
class CooccurrenceApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Domain Co-occurrence Explorer")

        # ---- Style
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass

        style.configure("TLabel", padding=(4, 2))
        style.configure("TButton", padding=(6, 4))
        style.configure("Header.TLabel", font=("Helvetica", 12, "bold"))

        # ---- State
        self.loaded_df = None
        self.extracted_header = []
        self.cleaned_data = None
        self.data_subset_binary = None
        self.protein_ids = None

        self.domain_group_vars = {}     # domain -> tk.StringVar
        self.domain_discard_vars = {}   # domain -> tk.BooleanVar
        self.group_assignments = {}     # domain -> group (Encapsulin/MCP/Uncategorized)

        self.pair_df = pd.DataFrame(columns=["domain1", "domain2", "number_proteins"])
        self.full_graph = None

        # Pair override store (domain1, domain2) -> int
        self.pair_overrides = {}

        # ---- Defaults
        self.default_input_path = "/Users/abelardoaguilar/projects/github_repos/mini-devel/mini-devel/results/verse_filesystem/Effect_of_DTRs_verse/data/Tree_Annotation/minimal_plus_no_DTR.csv"

        # ---- Top-level layout: Notebook tabs + status bar
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill="both", expand=True)

        self.status_var = tk.StringVar(value="Ready.")
        status = ttk.Label(root, textvariable=self.status_var, anchor="w")
        status.pack(fill="x", side="bottom")

        # ---- Tabs
        self.tab_load = ttk.Frame(self.notebook)
        self.tab_groups = ttk.Frame(self.notebook)
        self.tab_network = ttk.Frame(self.notebook)
        self.tab_evol = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_load, text="1) Load")
        self.notebook.add(self.tab_groups, text="2) Assign groups")
        self.notebook.add(self.tab_network, text="3) Co-occurrence")
        self.notebook.add(self.tab_evol, text="4) Evolutionary paths")

        self._build_load_tab()
        self._build_groups_tab()
        self._build_network_tab()
        self._build_evol_tab()

    # ============================================================
    # TAB 1: Load
    # ============================================================
    def _build_load_tab(self):
        frm = ttk.Frame(self.tab_load, padding=10)
        frm.pack(fill="both", expand=True)

        ttk.Label(frm, text="Load input CSV", style="Header.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 8), columnspan=3)

        self.input_path_var = tk.StringVar(value=self.default_input_path)
        ttk.Label(frm, text="Input CSV:").grid(row=1, column=0, sticky="e")
        ttk.Entry(frm, textvariable=self.input_path_var, width=90).grid(row=1, column=1, sticky="we", padx=(6, 6))
        ttk.Button(frm, text="Browse", command=self._browse_input).grid(row=1, column=2, sticky="w")

        self.data_start_row_var = tk.StringVar(value="5")
        ttk.Label(frm, text="Data start row (0-index):").grid(row=2, column=0, sticky="e", pady=(8, 0))
        ttk.Entry(frm, textvariable=self.data_start_row_var, width=10).grid(row=2, column=1, sticky="w", pady=(8, 0), padx=(6, 0))

        btns = ttk.Frame(frm)
        btns.grid(row=3, column=1, sticky="w", pady=(10, 10))
        ttk.Button(btns, text="Load file", command=self.load_file).pack(side="left")
        ttk.Button(btns, text="Go to group assignment →", command=lambda: self.notebook.select(self.tab_groups)).pack(side="left", padx=(10, 0))

        # Preview box
        ttk.Label(frm, text="Preview (top rows):").grid(row=4, column=0, sticky="ne", pady=(8, 0))
        self.preview_text = tk.Text(frm, height=16, width=110, wrap="none")
        self.preview_text.grid(row=4, column=1, columnspan=2, sticky="nsew", pady=(8, 0))

        frm.columnconfigure(1, weight=1)
        frm.rowconfigure(4, weight=1)

    def _browse_input(self):
        path = filedialog.askopenfilename(
            title="Select Input CSV File",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if path:
            self.input_path_var.set(path)

    def load_file(self):
        try:
            path = self.input_path_var.get().strip()
            if not path:
                messagebox.showerror("Error", "Please choose an input file.")
                return

            self.loaded_df = pd.read_csv(path, header=None)

            self.preview_text.delete("1.0", "end")
            self.preview_text.insert("end", self.loaded_df.head(12).to_string(index=False))
            self.status_var.set(f"Loaded file: {path}")

            # Try to auto-extract header + data immediately (so tab 2 can populate)
            self._prepare_clean_data()

            messagebox.showinfo("Loaded", "File loaded successfully.\nGo to 'Assign groups' to proceed.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file:\n{e}")
            self.status_var.set("Error loading file.")

    def _prepare_clean_data(self):
        """Extract domain header and numeric matrix from the loaded_df based on your CSV structure."""
        if self.loaded_df is None:
            return

        data_start_row = int(self.data_start_row_var.get())

        # Extract header from row 1, columns 2:-1
        header = self.loaded_df.iloc[1, 2:-1].tolist()
        self.extracted_header = header

        # Extract numeric data from data_start_row onward, columns 2:-1
        cleaned = self.loaded_df.iloc[data_start_row:, 2:-1].copy()
        cleaned.columns = header
        cleaned = cleaned.apply(pd.to_numeric, errors="coerce").fillna(0)
        cleaned = cleaned.replace(-1, 0)
        cleaned = cleaned.select_dtypes(include=["number"])

        self.cleaned_data = cleaned

        # Protein ids from column 0
        protein_ids = self.loaded_df.iloc[data_start_row:, 0].astype(str).tolist()
        if len(protein_ids) != len(cleaned):
            protein_ids = [f"ProteinRow_{i}" for i in range(len(cleaned))]
        self.protein_ids = protein_ids

    # ============================================================
    # TAB 2: Assign groups
    # ============================================================
    def _build_groups_tab(self):
        outer = ttk.Frame(self.tab_groups, padding=10)
        outer.pack(fill="both", expand=True)

        ttk.Label(outer, text="Assign domain groups", style="Header.TLabel").pack(anchor="w", pady=(0, 8))

        topbar = ttk.Frame(outer)
        topbar.pack(fill="x", pady=(0, 8))

        ttk.Button(topbar, text="Refresh domains from loaded file", command=self.populate_domain_list).pack(side="left")
        ttk.Button(topbar, text="Compute co-occurrence →", command=self.compute_pairs_and_go_network).pack(side="left", padx=(10, 0))

        note = ttk.Label(
            outer,
            text="Tip: Double-check the Data start row in tab 1 if domains look wrong.\n"
                 "Discarded domains will not be used for co-occurrence."
        )
        note.pack(anchor="w", pady=(0, 8))

        # Scrollable list of domains
        self.dom_scroll = ScrollableFrame(outer)
        self.dom_scroll.pack(fill="both", expand=True)

        # Header row inside scroll
        header_row = ttk.Frame(self.dom_scroll.inner)
        header_row.grid(row=0, column=0, sticky="we", pady=(0, 4))
        ttk.Label(header_row, text="Domain", width=40).grid(row=0, column=0, sticky="w")
        ttk.Label(header_row, text="Group", width=18).grid(row=0, column=1, sticky="w")
        ttk.Label(header_row, text="Discard", width=10).grid(row=0, column=2, sticky="w")

        self.populate_domain_list()

    def populate_domain_list(self):
        if self.loaded_df is None:
            messagebox.showwarning("No file", "Load a file first in tab 1.")
            return

        self._prepare_clean_data()
        domains = list(self.extracted_header)

        # Clear previous widgets except header row
        for w in self.dom_scroll.inner.winfo_children():
            # keep header_row at row 0 (we detect by children count, safer to just rebuild all)
            pass

        # Rebuild scroll content
        for w in self.dom_scroll.inner.winfo_children():
            w.destroy()

        header_row = ttk.Frame(self.dom_scroll.inner)
        header_row.grid(row=0, column=0, sticky="we", pady=(0, 4))
        ttk.Label(header_row, text="Domain", width=40).grid(row=0, column=0, sticky="w")
        ttk.Label(header_row, text="Group", width=18).grid(row=0, column=1, sticky="w")
        ttk.Label(header_row, text="Discard", width=10).grid(row=0, column=2, sticky="w")

        self.domain_group_vars = {}
        self.domain_discard_vars = {}

        for i, d in enumerate(domains, start=1):
            row = ttk.Frame(self.dom_scroll.inner)
            row.grid(row=i, column=0, sticky="we", pady=1)

            ttk.Label(row, text=d, width=40).grid(row=0, column=0, sticky="w")

            gv = tk.StringVar(value="MCP")
            self.domain_group_vars[d] = gv
            cb = ttk.Combobox(row, textvariable=gv, values=["Encapsulin", "MCP", "Uncategorized"], width=16, state="readonly")
            cb.grid(row=0, column=1, sticky="w")

            dv = tk.BooleanVar(value=False)
            self.domain_discard_vars[d] = dv
            ttk.Checkbutton(row, variable=dv).grid(row=0, column=2, sticky="w", padx=(6, 0))

        self.status_var.set(f"Loaded {len(domains)} domains for assignment.")

    def compute_pairs_and_go_network(self):
        try:
            self.compute_pair_table()
            self.notebook.select(self.tab_network)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # ============================================================
    # Pair counting & table creation
    # ============================================================
    def compute_pair_table(self):
        if self.cleaned_data is None:
            raise RuntimeError("No cleaned data available. Load the file first.")

        # Determine kept domains
        domains_to_keep = [d for d in self.cleaned_data.columns if d in self.domain_discard_vars and not self.domain_discard_vars[d].get()]
        if not domains_to_keep:
            raise RuntimeError("No domains selected (all discarded).")

        # Save group assignments
        self.group_assignments = {d: self.domain_group_vars[d].get() for d in domains_to_keep}

        # Binarize domain presence
        subset = self.cleaned_data[domains_to_keep].copy()
        self.data_subset_binary = (subset > 0).astype(int)

        # Compute pair counts: for each protein, count all combinations of present domains
        pair_counts = {}
        for i in range(len(self.data_subset_binary)):
            present = self.data_subset_binary.columns[self.data_subset_binary.iloc[i].values > 0].tolist()
            if len(present) < 2:
                continue
            present_sorted = sorted(present)
            for u, v in itertools.combinations(present_sorted, 2):
                pair_counts[(u, v)] = pair_counts.get((u, v), 0) + 1

        # Build pair_df
        rows = []
        for (u, v), c in sorted(pair_counts.items(), key=lambda x: (-x[1], x[0][0], x[0][1])):
            rows.append({"domain1": u, "domain2": v, "number_proteins": int(c)})

        self.pair_df = pd.DataFrame(rows, columns=["domain1", "domain2", "number_proteins"])

        # Reset overrides (fresh computation)
        self.pair_overrides = {(r["domain1"], r["domain2"]): int(r["number_proteins"]) for _, r in self.pair_df.iterrows()}

        # Push to UI
        self._populate_pair_treeview()
        self.draw_network_from_overrides()

        self.status_var.set(f"Computed {len(self.pair_df)} domain pairs from {len(self.data_subset_binary)} proteins.")

    # ============================================================
    # TAB 3: Co-occurrence (table + plot + export)
    # ============================================================
    def _build_network_tab(self):
        outer = ttk.Frame(self.tab_network, padding=10)
        outer.pack(fill="both", expand=True)

        ttk.Label(outer, text="Co-occurrence network (editable pair counts)", style="Header.TLabel").pack(anchor="w", pady=(0, 8))

        # Controls
        controls = ttk.Frame(outer)
        controls.pack(fill="x", pady=(0, 8))

        ttk.Button(controls, text="Recompute pairs from data", command=self.compute_pair_table).pack(side="left")
        ttk.Button(controls, text="Redraw network (use table overrides)", command=self.apply_table_overrides_and_redraw).pack(side="left", padx=(8, 0))
        ttk.Button(controls, text="Export SVG…", command=self.export_current_svg).pack(side="left", padx=(8, 0))

        # Split area: left table, right plot
        paned = ttk.Panedwindow(outer, orient="horizontal")
        paned.pack(fill="both", expand=True)

        left = ttk.Frame(paned)
        right = ttk.Frame(paned)
        paned.add(left, weight=1)
        paned.add(right, weight=2)

        # ---- Pair table
        ttk.Label(left, text="Pair summary (double-click count to edit):").pack(anchor="w")

        self.pair_tree = ttk.Treeview(left, columns=("domain1", "domain2", "number_proteins"), show="headings", height=18)
        self.pair_tree.heading("domain1", text="domain1")
        self.pair_tree.heading("domain2", text="domain2")
        self.pair_tree.heading("number_proteins", text="number_proteins")

        self.pair_tree.column("domain1", width=180, stretch=True)
        self.pair_tree.column("domain2", width=180, stretch=True)
        self.pair_tree.column("number_proteins", width=120, stretch=False)

        yscroll = ttk.Scrollbar(left, orient="vertical", command=self.pair_tree.yview)
        self.pair_tree.configure(yscrollcommand=yscroll.set)

        self.pair_tree.pack(side="left", fill="both", expand=True)
        yscroll.pack(side="right", fill="y")

        self.pair_tree.bind("<Double-1>", self._on_double_click_pair_tree)

        # ---- Plot area
        ttk.Label(right, text="Network preview:").pack(anchor="w")

        self.fig = plt.Figure(figsize=(7, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.axis("off")

        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill="both", expand=True)

        # Keep last drawn graph (for exporting)
        self.full_graph = None
        self.last_pos = None

    def _populate_pair_treeview(self):
        # clear
        for item in self.pair_tree.get_children():
            self.pair_tree.delete(item)

        if self.pair_df is None or self.pair_df.empty:
            return

        # insert
        for _, r in self.pair_df.iterrows():
            u, v = r["domain1"], r["domain2"]
            c = int(r["number_proteins"])
            self.pair_tree.insert("", "end", values=(u, v, c))

    def _on_double_click_pair_tree(self, event):
        """Allow editing number_proteins by overlaying an Entry widget."""
        region = self.pair_tree.identify("region", event.x, event.y)
        if region != "cell":
            return

        row_id = self.pair_tree.identify_row(event.y)
        col = self.pair_tree.identify_column(event.x)
        if not row_id or col != "#3":  # only edit number_proteins (3rd col)
            return

        x, y, width, height = self.pair_tree.bbox(row_id, col)
        old_val = self.pair_tree.set(row_id, "number_proteins")

        entry = ttk.Entry(self.pair_tree)
        entry.place(x=x, y=y, width=width, height=height)
        entry.insert(0, old_val)
        entry.focus_set()

        def finish_edit(_event=None):
            new_val = entry.get().strip()
            try:
                n = int(new_val)
                if n < 0:
                    raise ValueError
            except Exception:
                messagebox.showerror("Invalid value", "number_proteins must be a non-negative integer.")
                entry.destroy()
                return

            self.pair_tree.set(row_id, "number_proteins", str(n))
            entry.destroy()

        entry.bind("<Return>", finish_edit)
        entry.bind("<FocusOut>", finish_edit)

    def apply_table_overrides_and_redraw(self):
        """Read the Treeview table and update overrides, then redraw."""
        overrides = {}
        for item in self.pair_tree.get_children():
            u, v, c = self.pair_tree.item(item, "values")
            try:
                c_int = int(c)
            except Exception:
                c_int = 0
            # ensure consistent ordering
            u2, v2 = sorted((u, v))
            overrides[(u2, v2)] = max(0, c_int)

        self.pair_overrides = overrides
        self.draw_network_from_overrides()
        self.status_var.set("Redrew network using overridden pair counts.")

    def draw_network_from_overrides(self):
        """Build and draw the bipartite graph using the current pair_overrides counts."""
        if not self.pair_overrides:
            self.ax.clear()
            self.ax.axis("off")
            self.canvas.draw()
            return

        # Collect all domains present in pairs with count > 0
        domains = set()
        for (u, v), c in self.pair_overrides.items():
            if c > 0:
                domains.add(u)
                domains.add(v)

        # Build bipartite graph: domain nodes + synthetic protein nodes per pair
        G = nx.Graph()
        for d in sorted(domains):
            G.add_node(d, type="domain")

        # Protein node color determined by groups of the two domains (plus Uncategorized)
        def protein_color(u, v):
            gu = self.group_assignments.get(u, "Uncategorized")
            gv = self.group_assignments.get(v, "Uncategorized")
            if gu == "Uncategorized" or gv == "Uncategorized":
                return "gray"
            if gu != gv:
                return "red"
            return "green"

        # add protein nodes
        for (u, v), c in self.pair_overrides.items():
            if c <= 0:
                continue
            u, v = sorted((u, v))
            col = protein_color(u, v)
            for i in range(int(c)):
                pn = f"ProteinPair:{u}|{v}|{i+1}"
                G.add_node(pn, type="protein", color=col)
                G.add_edge(pn, u)
                G.add_edge(pn, v)

        # draw
        self.full_graph = G

        self.ax.clear()
        self.ax.axis("off")

        # layout (stable seed)
        pos = nx.spring_layout(G, seed=42)
        self.last_pos = pos

        domain_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "domain"]
        red_p = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "red"]
        green_p = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "green"]
        gray_p = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "gray"]

        node_sizes = [max(120, min(1200, G.degree(n) * 90)) for n in domain_nodes]

        nx.draw_networkx_nodes(G, pos, nodelist=domain_nodes, node_size=node_sizes, node_color="lightblue", ax=self.ax)
        nx.draw_networkx_nodes(G, pos, nodelist=red_p, node_size=6, node_color="red", ax=self.ax)
        nx.draw_networkx_nodes(G, pos, nodelist=green_p, node_size=5, node_color="lightgreen", ax=self.ax)
        nx.draw_networkx_nodes(G, pos, nodelist=gray_p, node_size=6, node_color="gray", ax=self.ax)

        nx.draw_networkx_edges(G, pos, width=0.25, alpha=0.8, ax=self.ax)
        nx.draw_networkx_labels(G, pos, labels={n: n for n in domain_nodes}, font_size=7, ax=self.ax)

        self.canvas.draw()

    def export_current_svg(self):
        if self.full_graph is None or self.last_pos is None:
            messagebox.showwarning("Nothing to export", "Draw the network first.")
            return

        out_path = filedialog.asksaveasfilename(
            title="Save network as SVG",
            defaultextension=".svg",
            filetypes=[("SVG Files", "*.svg"), ("All Files", "*.*")]
        )
        if not out_path:
            return

        try:
            # Create a clean export figure (don’t depend on the embedded canvas)
            fig = plt.Figure(figsize=(9, 7), dpi=120)
            ax = fig.add_subplot(111)
            ax.axis("off")

            G = self.full_graph
            pos = self.last_pos

            domain_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "domain"]
            red_p = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "red"]
            green_p = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "green"]
            gray_p = [n for n, d in G.nodes(data=True) if d.get("type") == "protein" and d.get("color") == "gray"]

            node_sizes = [max(120, min(1400, G.degree(n) * 95)) for n in domain_nodes]

            nx.draw_networkx_nodes(G, pos, nodelist=domain_nodes, node_size=node_sizes, node_color="lightblue", ax=ax)
            nx.draw_networkx_nodes(G, pos, nodelist=red_p, node_size=6, node_color="red", ax=ax)
            nx.draw_networkx_nodes(G, pos, nodelist=green_p, node_size=5, node_color="lightgreen", ax=ax)
            nx.draw_networkx_nodes(G, pos, nodelist=gray_p, node_size=6, node_color="gray", ax=ax)

            nx.draw_networkx_edges(G, pos, width=0.25, alpha=0.85, ax=ax)
            nx.draw_networkx_labels(G, pos, labels={n: n for n in domain_nodes}, font_size=8, ax=ax)

            fig.savefig(out_path, format="svg", bbox_inches="tight")
            self.status_var.set(f"Exported SVG: {out_path}")
            messagebox.showinfo("Exported", f"Saved:\n{out_path}")

        except Exception as e:
            messagebox.showerror("Export failed", str(e))

    # ============================================================
    # TAB 4: Evolutionary paths (kept but made more “button-driven”)
    # ============================================================
    def _build_evol_tab(self):
        outer = ttk.Frame(self.tab_evol, padding=10)
        outer.pack(fill="both", expand=True)

        ttk.Label(outer, text="Evolutionary paths (domain graph spanning trees)", style="Header.TLabel").pack(anchor="w", pady=(0, 8))

        ttk.Label(
            outer,
            text="This uses the current domain grouping + kept domains.\n"
                 "Edges exist when two domains co-occur in at least one protein.\n"
                 "Weights: 1 (same group), 5 (different groups), 3 (Uncategorized involved)."
        ).pack(anchor="w", pady=(0, 10))

        btns = ttk.Frame(outer)
        btns.pack(anchor="w", pady=(0, 10))
        ttk.Button(btns, text="Compute evolutionary paths", command=self.compute_evolutionary_paths).pack(side="left")
        ttk.Button(btns, text="Export selected path SVG…", command=self.export_selected_evol_svg).pack(side="left", padx=(8, 0))

        # Listbox for paths
        self.evol_list = tk.Listbox(outer, height=10)
        self.evol_list.pack(fill="x", pady=(0, 10))

        self.evol_paths = []  # list of (edge_set, cost, H_nodes)

        # Preview for evol (simple)
        self.evol_fig = plt.Figure(figsize=(7, 5), dpi=100)
        self.evol_ax = self.evol_fig.add_subplot(111)
        self.evol_ax.axis("off")
        self.evol_canvas = FigureCanvasTkAgg(self.evol_fig, master=outer)
        self.evol_canvas.get_tk_widget().pack(fill="both", expand=True)

        self.evol_list.bind("<<ListboxSelect>>", self._on_select_evol_path)

    def compute_evolutionary_paths(self):
        if self.data_subset_binary is None or not self.group_assignments:
            messagebox.showwarning("Missing data", "Compute co-occurrence first (tab 3).")
            return

        kept_domains = list(self.data_subset_binary.columns)
        if len(kept_domains) < 2:
            messagebox.showerror("Error", "Not enough domains to build a domain graph.")
            return

        H = nx.Graph()
        H.add_nodes_from(kept_domains)

        def edge_cost(u, v):
            gu = self.group_assignments.get(u, "Uncategorized")
            gv = self.group_assignments.get(v, "Uncategorized")
            if gu == "Uncategorized" or gv == "Uncategorized":
                return 3
            return 1 if gu == gv else 5

        # Build H from real protein co-occurrence (not overridden table)
        for i in range(len(self.data_subset_binary)):
            present = self.data_subset_binary.columns[self.data_subset_binary.iloc[i].values > 0].tolist()
            if len(present) < 2:
                continue
            for u, v in itertools.combinations(sorted(present), 2):
                w = edge_cost(u, v)
                if H.has_edge(u, v):
                    if w < H[u][v]["weight"]:
                        H[u][v]["weight"] = w
                else:
                    H.add_edge(u, v, weight=w)

        if not nx.is_connected(H):
            messagebox.showerror("Error", "Domain graph is not connected; cannot compute spanning trees.")
            return

        # Enumerate spanning trees for small graphs, else MST
        trees = []
        if len(kept_domains) > 12:
            mst = nx.minimum_spanning_tree(H, weight="weight")
            edge_set = set(tuple(sorted(e)) for e in mst.edges())
            cost = sum(H[u][v]["weight"] for u, v in edge_set)
            trees = [(edge_set, cost)]
        else:
            all_trees = self._enumerate_spanning_trees(H)
            for edge_set in all_trees:
                cost = sum(H[u][v]["weight"] for u, v in edge_set)
                trees.append((edge_set, cost))
            trees.sort(key=lambda x: x[1])

        # keep top 6
        self.evol_paths = []
        self.evol_list.delete(0, "end")

        for i, (edge_set, cost) in enumerate(trees[:6], start=1):
            self.evol_paths.append((edge_set, cost, kept_domains))
            self.evol_list.insert("end", f"Path {i}  |  cost = {cost}")

        if self.evol_paths:
            self.evol_list.selection_set(0)
            self._draw_evol_path(0)

        self.status_var.set(f"Computed {len(self.evol_paths)} evolutionary path(s).")

    def _enumerate_spanning_trees(self, H):
        """Backtracking enumeration (small graphs only). Returns list of edge sets."""
        trees = []
        nodes = list(H.nodes())
        if not nodes:
            return trees
        root = nodes[0]

        def backtrack(current_edges, used):
            if len(used) == len(nodes):
                trees.append(set(current_edges))
                return
            for u in list(used):
                for v in H.neighbors(u):
                    if v not in used:
                        edge = tuple(sorted((u, v)))
                        current_edges.append(edge)
                        used.add(v)
                        backtrack(current_edges, used)
                        used.remove(v)
                        current_edges.pop()

        backtrack([], set([root]))
        return trees

    def _on_select_evol_path(self, _event):
        sel = self.evol_list.curselection()
        if not sel:
            return
        self._draw_evol_path(sel[0])

    def _draw_evol_path(self, idx):
        if idx < 0 or idx >= len(self.evol_paths):
            return
        edge_set, cost, nodes = self.evol_paths[idx]

        T = nx.Graph()
        T.add_nodes_from(nodes)
        T.add_edges_from(edge_set)

        self.evol_ax.clear()
        self.evol_ax.axis("off")

        pos = nx.spring_layout(T, seed=42)
        nx.draw_networkx_nodes(T, pos, node_size=300, node_color="lightblue", ax=self.evol_ax)
        nx.draw_networkx_edges(T, pos, width=2, ax=self.evol_ax)
        nx.draw_networkx_labels(T, pos, font_size=8, ax=self.evol_ax)

        self.evol_ax.set_title(f"Evolutionary path (cost={cost})")
        self.evol_canvas.draw()

    def export_selected_evol_svg(self):
        sel = self.evol_list.curselection()
        if not sel:
            messagebox.showwarning("No selection", "Select a path first.")
            return

        idx = sel[0]
        edge_set, cost, nodes = self.evol_paths[idx]

        out_path = filedialog.asksaveasfilename(
            title="Save evolutionary path as SVG",
            defaultextension=".svg",
            filetypes=[("SVG Files", "*.svg"), ("All Files", "*.*")]
        )
        if not out_path:
            return

        try:
            T = nx.Graph()
            T.add_nodes_from(nodes)
            T.add_edges_from(edge_set)

            fig = plt.Figure(figsize=(8, 6), dpi=120)
            ax = fig.add_subplot(111)
            ax.axis("off")

            pos = nx.spring_layout(T, seed=42)
            nx.draw_networkx_nodes(T, pos, node_size=320, node_color="lightblue", ax=ax)
            nx.draw_networkx_edges(T, pos, width=2, ax=ax)
            nx.draw_networkx_labels(T, pos, font_size=8, ax=ax)

            ax.set_title(f"Evolutionary path (cost={cost})")
            fig.savefig(out_path, format="svg", bbox_inches="tight")

            self.status_var.set(f"Exported evolutionary path SVG: {out_path}")
            messagebox.showinfo("Exported", f"Saved:\n{out_path}")
        except Exception as e:
            messagebox.showerror("Export failed", str(e))


# ============================================================
# Run
# ============================================================
if __name__ == "__main__":
    root = tk.Tk()
    app = CooccurrenceApp(root)
    root.mainloop()


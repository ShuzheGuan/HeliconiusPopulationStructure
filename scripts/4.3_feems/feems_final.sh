#!/usr/bin/env bash
set -euo pipefail

# Purpose: fit the final FEEMS model with selected lambda and lambda_q.
# Input: LD-pruned PLINK files, sample coordinate file, and FEEMS grid shapefile.
# Output: fitted FEEMS object arrays and migration-surface plot.
# Software: FEEMS, pandas-plink, scikit-learn, matplotlib

PLINK_PREFIX="<path/to/ld_pruned_plink_prefix_without_extension>"
COORD_PATH="<path/to/sample_coordinates_lon_lat.txt>"
GRID_PATH="<path/to/grid.shp>"
OUT_DIR="<path/to/final_feems_output_dir>"
OUT_PREFIX="<output_prefix>"

BUFFER=0
THREADS=16
LAMBDA="<selected_lambda>"
LAMBDA_Q="<selected_lambda_q>"

mkdir -p "$OUT_DIR"

export MPLBACKEND=Agg
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS"
export PLINK_PREFIX COORD_PATH GRID_PATH OUT_DIR OUT_PREFIX BUFFER LAMBDA LAMBDA_Q

python3 - <<'PY'
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pandas_plink import read_plink
from sklearn.impute import SimpleImputer
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz

plink = os.environ["PLINK_PREFIX"]
coord = np.loadtxt(os.environ["COORD_PATH"])
out_dir = os.environ["OUT_DIR"]
prefix = os.environ["OUT_PREFIX"]

_, _, G = read_plink(plink)
G = SimpleImputer(missing_values=np.nan, strategy="mean").fit_transform(np.asarray(G).T).astype(np.float32)
outer, edges, grid_nodes, _ = prepare_graph_inputs(coord=coord, ggrid=os.environ["GRID_PATH"], translated=False, buffer=int(os.environ["BUFFER"]), outer=None)

graph = SpatialGraph(G, coord, grid_nodes, edges, scale_snps=True)
graph.fit(lamb=float(os.environ["LAMBDA"]), lamb_q=float(os.environ["LAMBDA_Q"]), optimize_q="n-dim")

np.savetxt(os.path.join(out_dir, f"{prefix}_grid_nodes.tsv"), grid_nodes, delimiter="\t")
np.savetxt(os.path.join(out_dir, f"{prefix}_edges.tsv"), edges, fmt="%d", delimiter="\t")
np.savetxt(os.path.join(out_dir, f"{prefix}_edge_weights.tsv"), graph.w, delimiter="\t")

fig, ax = plt.subplots(figsize=(7, 7), dpi=300)
viz = Viz(ax, graph, edge_width=0.6, edge_alpha=1.0, sample_pt_size=12, obs_node_size=5)
viz.draw_edges(use_weights=True)
viz.draw_obs_nodes(use_ids=False)
viz.draw_edge_colorbar()
fig.tight_layout()
fig.savefig(os.path.join(out_dir, f"{prefix}_feems_migration_surface.png"))
PY

#!/usr/bin/env bash
set -euo pipefail

# Purpose: cross-validate FEEMS lambda and lambda_q values.
# Input: LD-pruned PLINK files, sample coordinate file, and FEEMS grid shapefile.
# Output: CV matrix and best FEEMS parameter table.
# Software: FEEMS, pandas-plink, scikit-learn

PLINK_PREFIX="<path/to/ld_pruned_plink_prefix_without_extension>"
COORD_PATH="<path/to/sample_coordinates_lon_lat.txt>"
GRID_PATH="<path/to/grid.shp>"
OUT_DIR="<path/to/feems_parameter_testing_output_dir>"
OUT_PREFIX="<output_prefix>"

BUFFER=0
THREADS=16
N_FOLDS=5
FACTR="1e10"
LAMBDA_VALUES="1e-3,3.6e-3,1.3e-2,4.6e-2,1.7e-1,6e-1,2.2,7.7,28,100"
LAMBDA_Q_VALUES="1e-2,1e-1,1,10,100"

mkdir -p "$OUT_DIR"

export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS"
export PLINK_PREFIX COORD_PATH GRID_PATH OUT_DIR OUT_PREFIX BUFFER N_FOLDS FACTR LAMBDA_VALUES LAMBDA_Q_VALUES

python3 - <<'PY'
import os
import numpy as np
import pandas as pd
from pandas_plink import read_plink
from sklearn.impute import SimpleImputer
from feems.utils import prepare_graph_inputs
from feems.cross_validation import run_cv_joint
from feems import SpatialGraph

plink = os.environ["PLINK_PREFIX"]
coord = np.loadtxt(os.environ["COORD_PATH"])
out_dir = os.environ["OUT_DIR"]
prefix = os.environ["OUT_PREFIX"]

_, _, G = read_plink(plink)
G = SimpleImputer(missing_values=np.nan, strategy="mean").fit_transform(np.asarray(G).T).astype(np.float32)
outer, edges, grid_nodes, _ = prepare_graph_inputs(coord=coord, ggrid=os.environ["GRID_PATH"], translated=False, buffer=int(os.environ["BUFFER"]), outer=None)
graph = SpatialGraph(G, coord, grid_nodes, edges, scale_snps=True)

lambdas = np.array([float(x) for x in os.environ["LAMBDA_VALUES"].split(",")])
lambda_qs = np.array([float(x) for x in os.environ["LAMBDA_Q_VALUES"].split(",")])
cv = run_cv_joint(graph, lambdas, lambda_qs, n_folds=int(os.environ["N_FOLDS"]), factr=float(os.environ["FACTR"]))
mean_cv = cv.mean(axis=0)
iq, il = np.unravel_index(np.argmin(mean_cv), mean_cv.shape)

pd.DataFrame(mean_cv, index=lambda_qs, columns=lambdas).to_csv(
    os.path.join(out_dir, f"{prefix}_feems_cv.tsv"), sep="\t"
)
pd.DataFrame(
    [{"lambda": lambdas[il], "lambda_q": lambda_qs[iq], "mean_cv_error": mean_cv[iq, il]}]
).to_csv(os.path.join(out_dir, f"{prefix}_best_feems_parameters.tsv"), sep="\t", index=False)
PY

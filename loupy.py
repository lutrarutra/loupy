import os
import argparse
import uuid
from typing import Optional

import rpy2.robjects
import rpy2.robjects.packages
import scanpy as sc
import scipy
import rpy2

loupy_script_path = os.path.join(os.path.dirname(__file__), "loupy.R")


def export_to_cloupe(adata: sc.AnnData, output: str, projection_names: Optional[list[str]] = None, categorical_names: Optional[list[str]] = None, layer: Optional[str] = None):
    if projection_names is None:
        projection_names = [projection.removeprefix("X_") for projection in adata.obsm_keys() if projection != "X_pca"]

    if categorical_names is None:
        categorical_names = [col for col in adata.obs.columns if adata.obs[col].dtype.name == "category"]

    temp_name = str(uuid.uuid4())

    if layer is not None:
        adata.X = adata.layers[layer]

    if isinstance(adata.X, scipy.sparse.csc_matrix):
        adata.X = adata.X.tocsr()
        
    del adata.layers

    adata.write_h5ad(f"{temp_name}.h5ad")

    with open(loupy_script_path, "r") as f:
        r_code = f.read()
    
    r_loupy = rpy2.robjects.packages.STAP(r_code, "r_loupy")

    output_dir = os.path.dirname(output)
    output_name = os.path.basename(output).removesuffix(".cloupe")

    print("Exporting to Loupe Cell Browser's .cloupe file format...")
    try:
        r_loupy.create_loupe_file(temp_name, output_dir, output_name, rpy2.robjects.StrVector(projection_names), rpy2.robjects.StrVector(categorical_names))
        print("Done!")
    except Exception as e:
        print(e)

    print("Removing temporary files...")
    os.remove(f"{temp_name}.h5ad")
    if os.path.exists(f"{temp_name}.h5seurat"):
        os.remove(f"{temp_name}.h5seurat")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export clusters, annotations and projections from ScanPy's AnnData object to Loupe Cell Browser's .cloupe file format.")
    parser.add_argument("input", type=str, help="Input file")
    parser.add_argument("output", type=str, help="Output file")
    parser.add_argument("--projections", type=str, help="Comma-separated names (Without 'X_'-prefix) of the projections to be exported. Should be found in 'adata.obsm_keys()' and be with two dimensions. Exports all projections by default.", default=None)
    parser.add_argument("--categoricals", type=str, help="Comma-separated names of the categorical annotations to be exported. Should be found in 'adata.obs_keys()'. Exports all categorical features from adata.obs by default.", default=None)
    parser.add_argument("--layer", type=str, help="Name of the layer to be exported. If not provided, adata.X will be exported.", default=None)
    args = parser.parse_args()

    adata = sc.read(args.input)

    if args.projections is not None:
        projection_names = [projection.removeprefix("X_") for projection in args.projections.split(",")]
        for proj in projection_names:
            if f"X_{proj}" not in adata.obsm_keys():
                raise ValueError(f"{proj} is not found in 'adata.obsm'")
            if adata.obsm[f"X_{proj}"].shape[1] != 2:
                raise ValueError(f"{proj} is not a 2D projection.")
    else:
        projection_names = None

    if args.categoricals is not None:
        categorical_names = args.categoricals.split(",")
        for cat in categorical_names:
            if cat not in adata.obs_keys():
                raise ValueError(f"{cat} is not found in 'adata.obs'")
            if adata.obs[cat].dtype.name != "category":
                raise ValueError(f"{cat} is not a categorical feature. If it should be, you can use 'adata.obs['{cat}'] = adata.obs['{cat}'].astype('category')' to convert it to categorical dtype.")
    else:
        categorical_names = None

    if (layer := args.layer) is not None:
        if layer not in adata.layers.keys():
            raise ValueError(f"{layer} is not found in 'adata.layers'")

    export_to_cloupe(
        adata, args.output,
        projection_names=projection_names, categorical_names=categorical_names,
        layer=layer
    )

import os
import uuid
from importlib import resources

import anndata as ad
import rpy2
from rpy2.robjects import conversion, default_converter
import rpy2.robjects
import rpy2.robjects.packages


def export_to_cloupe(
    adata: ad.AnnData, output: str, projection_names: list[str] | None = None,
    categorical_names: list[str] | None = None, layer: str | None = None
):
    loupy_script_path = os.path.join(str(resources.files("loupy")), "loupy.R")
    import scipy
    import pandas as pd
    if projection_names is None:
        projection_names = [projection.removeprefix("X_") for projection in adata.obsm_keys() if projection != "X_pca"]

    if categorical_names is None:
        categorical_names = [col for col in adata.obs.columns if adata.obs[col].dtype.name == "category"]

    temp_name = str(uuid.uuid4())

    if layer is not None:
        adata.X = adata.layers[layer]

    if isinstance(adata.X, scipy.sparse.csc_matrix):
        adata.X = adata.X.tocsr()  # type: ignore
        
    del adata.layers

    for obsm in adata.obsm_keys():
        if isinstance(adata.obsm[obsm], pd.DataFrame):
            del adata.obsm[obsm]

    adata.write_h5ad(f"{temp_name}.h5ad")

    with open(loupy_script_path, "r") as f:
        r_code = f.read()
        
    with conversion.localconverter(default_converter):
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

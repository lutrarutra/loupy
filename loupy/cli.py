import anndata as ad
import argparse

from .loupy import export_to_cloupe


def main():
    parser = argparse.ArgumentParser(description="Export clusters, annotations and projections from ScanPy's AnnData object to Loupe Cell Browser's .cloupe file format.")
    parser.add_argument("input", type=str, help="Input file (Anndata .h5ad-file)")
    parser.add_argument("output", type=str, help="Output file (.cloupe)")
    parser.add_argument("--projections", type=str, help="Comma-separated names (Without 'X_'-prefix) of the projections to be exported. Should be found in 'adata.obsm_keys()' and be with two dimensions. Exports all projections by default.", default=None)
    parser.add_argument("--categoricals", type=str, help="Comma-separated names of the categorical annotations to be exported. Should be found in 'adata.obs_keys()'. Exports all categorical features from adata.obs by default.", default=None)
    parser.add_argument("--layer", type=str, help="Name of the layer to be exported. If not provided, adata.X will be exported.", default=None)
    args = parser.parse_args()

    adata = ad.read_h5ad(args.input)

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


if __name__ == "__main__":
    main()
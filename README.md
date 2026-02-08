# CoExpress

**Gene Expression Correlation Explorer**

Explore co-expression patterns across ~1,700 cancer cell lines using RNA-seq data (log2 TPM+1) from [DepMap](https://depmap.org).

**[Launch CoExpress](https://fredrikwermeling.github.io/CoExpress/)**

## What it does

CoExpress finds genes with correlated expression patterns across cancer cell lines. Enter a set of genes and discover which ones are co-expressed, visualize correlation networks, and explore expression distributions by cancer type or mutation status.

## Features

- **Correlation analysis** within a gene list or genome-wide discovery mode
- **Interactive network** visualization with customizable layout and coloring
- **Scatter plots** for any gene pair with tissue and mutation overlays
- **Mutation analysis** comparing expression between wild-type and mutant cell lines
- **Expression distributions** by cancer type with statistical testing
- **Synonym/ortholog lookup** for gene name resolution
- **Export** networks (PNG/SVG), tables (CSV), and full results (ZIP)

## Data

All data from [DepMap 25Q3](https://depmap.org/portal/data_page/?tab=currentRelease):

- **Expression:** OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv (19,215 genes x 1,699 cell lines)
- **Metadata:** Model.csv (cell line annotations, lineage, disease)
- **Mutations:** OmicsSomaticMutationsMatrixHotspot.csv (49 hotspot genes)

Please [acknowledge DepMap](https://depmap.org/portal/data_page/?tab=overview#how-to-cite) if you use this data.

## Related tools

- [Correlate](https://fredrikwermeling.github.io/correlation-app/) - Same analysis engine for CRISPR gene effect data
- [Green Listed](https://greenlisted.cmm.se) - Essential gene explorer

## License

MIT - [Wermeling Lab](https://wermelinglab.com), Karolinska Institutet

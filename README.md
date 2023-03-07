# Adjustments to the reference dataset design improve cell type label transfer
This is the GitHub repository holding the scripts and result files for the paper "Adjustments to the reference dataset design improve cell type label transfer" by Moelbert et al.  The preprint can be found on BioRxiv (https://www.biorxiv.org/content/10.1101/2023.01.25.525533v1).

The transfer of cell type labels from prior annotated (reference) to newly collected data is an important task in single-cell data analysis. As the number of publicly available annotated datasets which can be used as a reference, as well as the number of computational methods for cell type label transfer are constantly growing, rationals to understand and decide which reference design and which method to use for a particular query dataset is needed. Here, we benchmark a set of five popular cell type annotation methods, study the performance on different cell types and highlight the importance of the design of the reference data (number of cell samples for each cell type, inclusion of multiple datasets in one reference, gene set selection, etc.) for more reliable predictions. 


We introduce a weighted bootstrapping-based approach that allows to make use of the entrie reference dataset while still keeping the benefits from not under-representing any of the existing cell types. 

![bootstrapping_algorithm](https://user-images.githubusercontent.com/35259165/223424034-ddfd186d-db55-47bc-b9a7-4249bfc130b5.jpg)

The approach can be applied in combination with any lable transfer method and is especially helpful when working with model-based approaches. 

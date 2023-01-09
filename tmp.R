import(CellChat)
import(ComplexHeatmap)
import(NMF)
import(Seurat)
import(SeuratObject)
import(circlize)
import(data.table)
import(ggalluvial)
import(ggplot2)

usethis::use_package("CellChat")
usethis::use_package("ComplexHeatmap")
usethis::use_package("NMF")
usethis::use_package("Seurat")
usethis::use_package("SeuratObject")
usethis::use_package("circlize")
usethis::use_package("data.table")
usethis::use_package("ggalluvial")
usethis::use_package("ggplot2")



usethis::use_vignette("Update_NeuronChat_database")
usethis::use_vignette("Spatial_analysis")
usethis::use_readme_md(open = rlang::is_interactive())
rmarkdown::html_document()

roxygen2::roxygenise()

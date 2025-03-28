# Small pancreatic Neuroendocrine Tumors
R package to predict probability of metastatic pancreatic NET from CT scan features



# Install
```R
install.packages("BiocManager")
BiocManager::install(c("devtools"))
devtools::install_github("GeNeHetX/smallpanetradio")

```

# Web app
To run the interactive Shiny webapp :
```R
library(smallpanetradio)
runShinyApp()
```


# Warning

**This tool is intended for research purposes only and should not be used for clinical or diagnostic decision-making.**
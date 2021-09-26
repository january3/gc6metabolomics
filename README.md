To compile the manuscript, please use the following code.


```
library(remotes)
install_github("january3/myfuncs")
library(myfuncs)
source("functions/manuscript_helpers.R")
manu_process("manuscript.rmd", "manuscript_out.rmd")
rmarkdown::render("manuscript_out.rmd")
```



active.analysis <- "all_domains"
markdown.dir <- file.path("markdown",active.analysis)
if(!dir.exists(markdown.dir)){
  dir.create(markdown.dir)
}

knitr::spin("code/r-scripts/001-phyloseq-qiime2.R", knit = FALSE)
file.rename("code/r-scripts/001-phyloseq-qiime2.Rmd",
            file.path(markdown.dir,"001-phyloseq-qiime2.Rmd"))

rmarkdown::render(file.path(markdown.dir,"001-phyloseq-qiime2.Rmd"),
                  'html_document',
                  knit_root_dir=here::here(),
                  params=list(active.analysis=active.analysis))

knitr::spin("code/r-scripts/002-summary-stats-qiime2.R", knit = FALSE)
file.rename("code/r-scripts/002-summary-stats-qiime2.Rmd",
            file.path(markdown.dir,"002-summary-stats-qiime2.Rmd"))
rmarkdown::render(file.path(markdown.dir,"002-summary-stats-qiime2.Rmd"), 
                  'html_document',
                  knit_root_dir=here::here(),
                  params=list(active.analysis=active.analysis))

knitr::spin("code/r-scripts/003-compare-data.R", knit = FALSE)
file.rename("code/r-scripts/003-compare-data.Rmd",
            "markdown/003-compare-data.Rmd")
rmarkdown::render('./markdown/003-compare-data.Rmd', 'html_document',
                  knit_root_dir=here::here())

knitr::spin("code/r-scripts/004-barplots-qiime2.R", knit = FALSE)
file.rename("code/r-scripts/004-barplots-qiime2.Rmd",
            file.path(markdown.dir, "004-barplots-qiime2.Rmd"))
rmarkdown::render(file.path(markdown.dir, "004-barplots-qiime2.Rmd"), 
                  'html_document',
                  knit_root_dir=here::here(),
                  params=list(active.analysis=active.analysis))

knitr::spin("code/r-scripts/005-pca.R", knit = FALSE)
file.rename("code/r-scripts/005-pca.Rmd", 
            file.path(markdown.dir,"005-pca.Rmd"))
rmarkdown::render(file.path(markdown.dir,"005-pca.Rmd"),
                  'html_document',
                  knit_root_dir=here::here(),
                  params=list(active.analysis=active.analysis))

knitr::spin("code/r-scripts/006-alpha-diversity.R", knit = FALSE)
file.rename("code/r-scripts/006-alpha-diversity.Rmd",
            file.path(markdown.dir,"006-alpha-diversity.Rmd"))
rmarkdown::render(file.path(markdown.dir,"006-alpha-diversity.Rmd"), 
                  'html_document',
                  knit_root_dir=here::here(),
                  params=list(active.analysis=active.analysis))

knitr::spin("code/r-scripts/007-diffabund-tests.R", knit = FALSE)
file.rename("code/r-scripts/007-diffabund-tests.Rmd",
            file.path(markdown.dir,"007-diffabund-tests.Rmd"))

rmarkdown::render(file.path(markdown.dir,"007-diffabund-tests.Rmd"), 
                  'html_document',
                  knit_root_dir=here::here(),
                  params=list(active.analysis=active.analysis))

knitr::spin("code/r-scripts/008-diversity-inside-custom-host.R", knit = FALSE)
file.rename("code/r-scripts/008-diversity-inside-custom-host.Rmd",
            file.path(markdown.dir,"008-diversity-inside-custom-host.Rmd"))

rmarkdown::render(file.path(markdown.dir,"008-diversity-inside-custom-host.Rmd"), 
                  'html_document',
                  knit_root_dir=here::here(),
                  params=list(active.analysis=active.analysis))


render_book(input = file.path(markdown.dir), 
            output_format = "bookdown::html_document2", clean = TRUE,
            envir = parent.frame(),
            output_dir =file.path("markdown", paste(active.analysis,"16S_bind", sep = "_")) , 
            new_session = TRUE, preview = FALSE,
            config_file = file.path(markdown.dir,"_bookdown.yml"))

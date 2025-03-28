FROM rocker/shiny:latest
COPY . /srv/shiny-app
RUN R -e "install.packages(c('ggpubr','caret','FactoMineR','glmnet','openxlsx'), repos='http://cran.rstudio.com/')"
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/srv/shiny-app')"]


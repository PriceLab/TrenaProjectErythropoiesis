FROM rocker/shiny
WORKDIR /home/shiny
COPY installTrenaCodePrequisites.R .
RUN R -f ./installTrenaCodePrequisites.R
COPY linePlot.R .
COPY linePlot.js .
COPY linePlot.css .
COPY srm.rna.averaged.clean.RData .
CMD ["R", "-e", "shiny::runApp('/home/shiny/linePlot.R', port=3838, host='0.0.0.0')"]






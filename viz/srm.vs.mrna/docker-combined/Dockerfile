FROM rocker/shiny
WORKDIR /home/shiny
COPY installTrenaCodePrequisites.R .
RUN R -f ./installTrenaCodePrequisites.R
COPY combined.R .
COPY multiPlot.js .
COPY multiPlot.css .
COPY linePlot.js .
COPY linePlot.css .
COPY srm.rna.averaged.clean.RData .
CMD ["R", "-e", "shiny::runApp('/home/shiny/combined.R', port=3838, host='0.0.0.0')"]






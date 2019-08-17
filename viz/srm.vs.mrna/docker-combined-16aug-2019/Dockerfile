FROM rocker/shiny
WORKDIR /home/shiny
COPY installTrenaCodePrequisites.R .
RUN R -f ./installTrenaCodePrequisites.R
RUN mkdir www
COPY www/* www/
COPY app.R .
COPY multiPlot.js .
COPY multiPlot.css .
COPY linePlot.js .
COPY linePlot.css .
COPY srm.rna.averaged.clean.RData .
CMD ["R", "-e", "shiny::runApp('/home/shiny/app.R', port=3838, host='0.0.0.0')"]






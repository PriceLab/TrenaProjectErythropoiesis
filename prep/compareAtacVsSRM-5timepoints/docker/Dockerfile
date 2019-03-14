FROM pshannon/trenashinybase
USER root
WORKDIR /home/trena/app
COPY simpleShinyApp.R .
COPY KLF1.FLI1.atac.srm.RData .
COPY README.html .
USER trena
CMD ["R", "-e", "shiny::runApp('/home/trena/app/simpleShinyApp.R', port=3838, host='0.0.0.0')"]






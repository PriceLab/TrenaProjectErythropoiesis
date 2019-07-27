//------------------------------------------------------------------------------------------------------------------------
r2d3.onRender(function(data, svg, width, height, options){

   r2d3.svg.selectAll("g").remove()

    var text = r2d3.svg.selectAll("text")
                       .data("placeholder")
                       .enter()
                       .append("text");

   rnaData = data.rna
   srmData = data.srm
   xMax = data.xMax
   yMax = data.yMax
  
   var d3Div = document.getElementById("d3");
   var actual_width = d3Div.clientWidth;
   var actual_height = d3Div.clientHeight;
  
   width = actual_width * 0.95;
   height = actual_height * 0.95;
   margin = 40;


   // console.log("cmd: " + data.cmd)
   var vectorNames = Object.keys(data.vectors)

   var xScalingFunction = d3.scaleLinear()
       .domain([0, xMax * 1.1])  // the range of the values to plot
       .range([0, width]);             // the pixel range of the x-axis
    
   var yScalingFunction = d3.scaleLinear()
       .domain([0, yMax * 1.1])
       .range([height, 0]);
    
   var lineGenerator = d3.line()
       .x(function(d, i) {return xScalingFunction(d.x);})
       .y(function(d)    {return yScalingFunction(d.y);})
       .curve(d3.curveCardinal)
       //.curve(d3.curveLinear) // apply smoothing to the line

   var xAxis = d3.axisBottom()
       .scale(xScalingFunction);
   
   var yAxis = d3.axisLeft()
        .scale(yScalingFunction);
    
    //------------------------------
    // axes
    //------------------------------
    
    xShift = margin;
    yShift = height;
    translationString = `translate(${xShift}, ${yShift})`
    
    r2d3.svg.append('g')
        .attr('transform', translationString)
        .call(xAxis);
    
    xShift = margin;
    yShift = 0
    translationString = `translate(${xShift}, ${yShift})`
    
    r2d3.svg.append('g')
        .attr('transform', translationString)
        .call(yAxis);
    
    //------------------------------
    // the plotting surface
    //------------------------------
    
    xShift = margin;
    yShift = 0;
    translationString = `translate(${xShift}, ${yShift})`
    var plottingSurface = r2d3.svg.append('g')
        .attr('transform', translationString)
        .attr('width', width)
        .attr('height', height)
        .attr('class', 'plottingSurface')   
    
    var margin = {top: 50, right: 50, bottom: 50, left: 50}
    var width = window.innerWidth - margin.left - margin.right // Use the window's width 
    var height = window.innerHeight - margin.top - margin.bottom; // Use the window's height
       
    var colorNumber = 0;
    var colors = d3.schemeCategory10;
    var colorCount = colors.length;
    
    for(vectorName of vectorNames){
       colorNumber = colorNumber + 1
       console.log("adding " + vectorName);
       var oneDataset = data.vectors[vectorName];

       plottingSurface.append("path")
         .datum(oneDataset)
         .attr("d", lineGenerator)
         .attr("fill", "none")
         .attr("stroke-width", 2)
         .attr("vName", vectorName)
         .attr("stroke", colors[(colorNumber % colorCount) - 1])
         .on("mouseover", function(d, i){
             console.log("mouse over: " + d3.select(this).attr("vName"));
             d3.select(this).attr('stroke-width', 10);
             Shiny.setInputValue("currentlySelectedVector", d3.select(this).attr("vName"));
             })
           .on("mouseout", function(d, i) {
              d3.select(this).attr('stroke-width', 2)
              Shiny.setInputValue("currentlySelectedVector", " ");
              })


       plottingSurface.selectAll("dot")
         .data(oneDataset)
         .enter().append("circle") // Uses the enter().append() method
           .attr("class", "dot") // Assign a class for styling
           .attr("cx", function(d) {return xScalingFunction(d.x)})
           .attr("cy", function(d) {return yScalingFunction(d.y)})
           .attr("r", 10)
        
           //.on("mouseover", function(d,i) {
           //   console.log("mouse!");
           //   d3.select(this).append("text")
           //  .text(vectorName)
           //  .attr("x", x(xScalingFunction(d.x)))
           //  .attr("y", y(yScalingFunction(d.y))); 
           //  });


       } // for vectorName

}) // onRender
//------------------------------------------------------------------------------------------------------------------------


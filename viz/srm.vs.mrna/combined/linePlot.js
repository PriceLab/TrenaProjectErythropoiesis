//------------------------------------------------------------------------------------------------------------------------
r2d3.onRender(function(data, svg, width, height, options){

   //console.log("=== linePlot.js, r2d3.onRender")
    
   r2d3.svg.selectAll("g").remove()

   rnaData = data.rna
   srmData = data.srm
   xMax = data.xMax
   yMax = data.yMax
   y2Max = data.y2Max
  
   var d3Div = document.getElementById("srm.rna.d3");
   var actual_width = d3Div.clientWidth;
   var actual_height = d3Div.clientHeight;
  
   sideMargin = 150; // 40;
   bottomMargin = 50;
   topMargin = 20;

   //width = actual_width * 0.80; // 0.95
   width = actual_width - (2 * sideMargin);
   height = actual_height - (1 * (bottomMargin + topMargin)); //* 0.90;
   //height = actual_height * 0.90;

   // console.log("cmd: " + data.cmd)
   var dataset = data.srm

   var xScalingFunction = d3.scaleLinear()
       .domain([0, xMax * 1.0])  // the range of the values to plot
       .range([0, width]);             // the pixel range of the x-axis
    
   var yScalingFunction = d3.scaleLinear()
       .domain([0, yMax * 1.05])
       .range([height, 0]);

   var y2ScalingFunction = d3.scaleLinear()
       .domain([0, y2Max * 1.05])
       .range([height, 0]);
    
   var lineFunction1 = d3.line()
       .x(function(d, i) { return xScalingFunction(d.x); }) // set the x values for the line generator
       .y(function(d) { return yScalingFunction(d.y); }) // set the y values for the line generator 
       //.curve(d3.curveMonotoneX)
       .curve(d3.curveLinear) // apply smoothing to the line

   var lineFunction2 = d3.line()
       .x(function(d, i) { return xScalingFunction(d.x); }) // set the x values for the line generator
       .y(function(d) { return y2ScalingFunction(d.y); }) // set the y values for the line generator 
       //.curve(d3.curveMonotoneX)
       .curve(d3.curveLinear) // apply smoothing to the line

   var xAxis = d3.axisBottom()
       .scale(xScalingFunction);
   
   var yAxis = d3.axisLeft()
        .scale(yScalingFunction);

   var y2Axis = d3.axisRight()
        .scale(y2ScalingFunction);
    
    //------------------------------
    // axes
    //------------------------------
    
    xShift = sideMargin;
    yShift = height;
    translationString = `translate(${xShift}, ${yShift})`
    
    r2d3.svg.append('g')
        .attr("class", "axis")
        .attr('transform', translationString)
        .call(xAxis);


    xShift = sideMargin;
    yShift = 0
    translationString = `translate(${xShift}, ${yShift})`
    
    r2d3.svg.append('g')
        .attr("class", "axis")
        //.style("stroke", "blue")
        .attr('transform', translationString)
        .call(yAxis);

    xShift = sideMargin + width;
    yShift = 0;
    translationString = `translate(${xShift}, ${yShift})`
    
    r2d3.svg.append('g')
        .attr("class", "axis")
        //.style("stroke", "red")
        .attr('transform', translationString)
        .call(y2Axis);
    

    //------------------------------
    // axis labels
    //------------------------------
    let xPos = 30;
    let yPos = height/2;
    translationString = `translate(${xPos}, ${yPos})`
    r2d3.svg.append("text")             
      .style("font-size", 14)
      .attr("transform", translationString)
      .style("text-anchor", "middle")
      .style("stroke", "blue")
      .text("Protein");

    xPos = width + (1.5 * sideMargin);
    yPos = height/2;
    translationString = `translate(${xPos}, ${yPos})`
    r2d3.svg.append("text")             
      .style("font-size", 14)
      .attr("transform", translationString)
      .style("text-anchor", "middle")
      .style("stroke", "red")
      .text("RNA");

    
    xPos = (width/2) + sideMargin;
    yPos = height + (bottomMargin);
    translationString = `translate(${xPos}, ${yPos})`
    r2d3.svg.append("text")             
      .style("font-size", 18)
      .attr("transform", translationString)
      .style("text-anchor", "middle")
      .text("Day");

      // Day0: MPP
      // Day2 to Day4: MEP
      // Day4 to Day11.5: CFU-E
      // Day11.5 : ProEB
      // Day14: BasoEB
    
    let developmentalStage_fontSize = 24;
    let developmentalStage_fontStyle = "normal";
    let developmentalStage_fontColor = "red";
    let developmentalStage_textAnchor = "middle"

   xPos = sideMargin;
    yPos = height + (bottomMargin) + 20;
    translationString = `translate(${xPos}, ${yPos})`
    r2d3.svg.append("text")             
      .style("font-size", developmentalStage_fontSize)
      .style("font-style", developmentalStage_fontStyle)
      .style("font-color", developmentalStage_fontColor)
      .style("text-anchor", developmentalStage_textAnchor)
      .attr("transform", translationString)
      .text("MPP");

    xPos = sideMargin + (3 * (width/14));
    yPos = height + (bottomMargin) + 20;
    translationString = `translate(${xPos}, ${yPos})`
    r2d3.svg.append("text")             
      .style("font-size", developmentalStage_fontSize)
      .style("font-style", developmentalStage_fontStyle)
      .style("font-color", developmentalStage_fontColor)
      .style("text-anchor", developmentalStage_textAnchor)
      .attr("transform", translationString)
      .text("MEP");

    xPos = sideMargin + (7.75 * (width/14));
    yPos = height + (bottomMargin) + 20;
    translationString = `translate(${xPos}, ${yPos})`
    r2d3.svg.append("text")             
      .style("font-size", developmentalStage_fontSize)
      .style("font-style", developmentalStage_fontStyle)
      .style("font-color", developmentalStage_fontColor)
      .style("text-anchor", developmentalStage_textAnchor)
      .attr("transform", translationString)
      .text("CFU-E");

    xPos = sideMargin + (11.5 * (width/14));
    yPos = height + (bottomMargin) + 20;
    translationString = `translate(${xPos}, ${yPos})`
    r2d3.svg.append("text")             
      .style("font-size", developmentalStage_fontSize)
      .style("font-style", developmentalStage_fontStyle)
      .style("font-color", developmentalStage_fontColor)
      .style("text-anchor", developmentalStage_textAnchor)
      .attr("transform", translationString)
      .text("ProEB");

    xPos = sideMargin + (14 * (width/14));
    yPos = height + (bottomMargin) + 20;
    translationString = `translate(${xPos}, ${yPos})`
    r2d3.svg.append("text")             
      .style("font-size", developmentalStage_fontSize)
      .style("font-style", developmentalStage_fontStyle)
      .style("font-color", developmentalStage_fontColor)
      .style("text-anchor", developmentalStage_textAnchor)
      .attr("transform", translationString)
      .text("BasoEB");



    //------------------------------
    // the plotting surface
    //------------------------------
    
    xShift = sideMargin;
    yShift = 0;
    translationString = `translate(${xShift}, ${yShift})`
    
    var plottingSurface = r2d3.svg.append('g')
        .attr('transform', translationString)
        .attr('width', width)
        .attr('height', height)
        .attr('class', 'plottingSurface')   

    //console.log("=== about to remove all from plottingSurface")
    //plottingSurface.selectAll("*").remove();

    //var margin = {top: 50, right: 50, bottom: 50, left: 50}
    //var width = window.innerWidth - margin.left - margin.right // Use the window's width 
    //var height = window.innerHeight - margin.top - margin.bottom; // Use the window's height
    
    
    plottingSurface.append("path")
        .datum(dataset) // 10. Binds data to the line 
        .attr("class", "line_srm") // Assign a class for styling 
        .attr("d", lineFunction1); // 11. Calls the line generator 
    
    plottingSurface.selectAll("dot")
        .data(dataset)
        .enter().append("circle") // Uses the enter().append() method
        .attr("class", "dot") // Assign a class for styling
        .attr("cx", function(d) { return xScalingFunction(d.x) })
        .attr("cy", function(d) { return yScalingFunction(d.y) })
        .attr("r", 5)
    
    var dataset = data.rna
    
    plottingSurface.append("path")
        .datum(dataset) // 10. Binds data to the line 
        .attr("class", "line_rna") // Assign a class for styling 
        .attr("d", lineFunction2); // 11. Calls the line generator 
    
    plottingSurface.selectAll("dot")
        .data(dataset)
        .enter().append("circle") // Uses the enter().append() method
        .attr("class", "dotrna") // Assign a class for styling
        .attr("cx", function(d) { return xScalingFunction(d.x) })
        .attr("cy", function(d) { return y2ScalingFunction(d.y) })
        .attr("r", 5)

}) // onRender
//------------------------------------------------------------------------------------------------------------------------


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

var xScalingFunction = d3.scaleLinear()
  .domain([0, xMax * 1.1])  // the range of the values to plot
  .range([0, width]);             // the pixel range of the x-axis
 
var yScalingFunction = d3.scaleLinear()
     .domain([0, yMax * 1.1])
     .range([height, 0]);

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

var n = 21;

var line = d3.line()
   .x(function(d, i) { return xScalingFunction(d.x); }) // set the x values for the line generator
   .y(function(d) { return yScalingFunction(d.y); }) // set the y values for the line generator 
   .curve(d3.curveLinear) // apply smoothing to the line

//------------------------------------------------------------------------------------------------------------------------
r2d3.onRender(function(data, svg, width, height, options){

  console.log("cmd: " + data.cmd)
  var dataset = data.srm

  plottingSurface.selectAll("*").remove()
    
  plottingSurface.append("path")
    .datum(dataset) // 10. Binds data to the line 
    .attr("class", "line_srm") // Assign a class for styling 
    .attr("d", line); // 11. Calls the line generator 

  plottingSurface.selectAll("dot")
    .data(dataset)
    .enter().append("circle") // Uses the enter().append() method
      .attr("class", "dot") // Assign a class for styling
      .attr("cx", function(d) { return xScalingFunction(d.x) })
      .attr("cy", function(d) { return yScalingFunction(d.y) })
      .attr("r", 10)

  var dataset = data.rna

  plottingSurface.append("path")
    .datum(dataset) // 10. Binds data to the line 
    .attr("class", "line_rna") // Assign a class for styling 
    .attr("d", line); // 11. Calls the line generator 

  plottingSurface.selectAll("dot")
   .data(dataset)
   .enter().append("circle") // Uses the enter().append() method
    .attr("class", "dotrna") // Assign a class for styling
    .attr("cx", function(d) { return xScalingFunction(d.x) })
    .attr("cy", function(d) { return yScalingFunction(d.y) })
    .attr("r", 10)

}) // onRender
//------------------------------------------------------------------------------------------------------------------------


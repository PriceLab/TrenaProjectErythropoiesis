
plotData = data.values;
xMax = data.xMax
yMax = data.yMax

// size and margins for the chart

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
    // the x axis
    //------------------------------
    
    xShift = margin;
    yShift = height;
    translationString = `translate(${xShift}, ${yShift})`
    
    r2d3.svg.append('g')
        .attr('transform', translationString)
        .call(xAxis);
    
    //------------------------------
    // the y axis
    //------------------------------
    
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
    

//------------------------------------------------------------------------------------------------------------------------
r2d3.onRender(function(data, svg, width, height, options){

   var actual_width = d3Div.clientWidth;
   var actual_height = d3Div.clientHeight;
   var width = actual_width * 0.95;
   var height = actual_height * 0.95;

   var plotData = data.values;
   var xMax = data.xMax
   var yMax = data.yMax
       
   var circle = plottingSurface.selectAll("circle").data(plotData);
   circle.exit().remove();
    
   var focalEntry = plotData.length;  // highlight this entry

   circle.enter().append("circle")
       .attr("r", 10)
     .merge(circle)
       .attr("cx", function (d) {return xScalingFunction(d.x);})
       .attr("cy", function (d) {return yScalingFunction(d.y);})
        .style("opacity", function(d, i){
            //console.log("i: " + i + "  length of data: " + plotData.length);
            opacity = 0.3;
            if((focalEntry-1) == i){
               opacity = 1.0
               }
            //console.log("opacity: " + opacity)
            return(opacity)
         }) // style
    
   }) // onRender
//------------------------------------------------------------------------------------------------------------------------
              

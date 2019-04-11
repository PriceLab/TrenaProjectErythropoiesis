vizmap = [

   {selector: "node", css: {
      "shape": "ellipse",
      "text-valign":"center",
      "text-halign":"center",
      "content": "data(label)",
      "background-color": "#FFFFFF",
      "border-color":"black","border-width":"1px",
      "width": "mapData(degree, 0.0, 100.0, 20.0, 200.0)",
      "height":"mapData(degree, 0.0, 100.0, 20.0, 200.0)",
      "font-size":"32px"}},


   {selector: "node[type='info']", css: {
       "shape": "roundrectangle",
       "font-size": "72px",
       "width": "360px",
       "height": "120px",
       "border-width": "3px",
       "background-color": "beige"
       }},

   {selector: "node[type='regulatoryRegion']", css: {
       "shape": "roundrectangle",
       "width": "5px",
       "height": "30px",
       "font-size": "0px",
       "background-color": "lightblue"
       }},

   {selector: "node[type='regulatoryRegion']:selected", css: {
       "shape": "roundrectangle",
       "width": "10px",
       "height": "30px",
       "font-size": "40px",
       "background-color": "beige"
       }},

   {selector: "node[type='TF'][pearson>0]", css: {
      "shape": "ellipse",
      "background-color": "mapData(pearson, 0, 1.0, white, red)",
      "width": "mapData(expression, 0, 1.0, 20.0, 200.0)",
      "height":"mapData(expression, 0, 1.0, 20.0, 200.0)"
       }},

   {selector: "node[type='targetGene']", css: {
      "shape": "roundrectangle",
      "background-color": "white",
      "width": "mapData(expression, 0.0, 1.0, 20.0, 200.0)",
      "height":"mapData(expression, 0.0, 1.0, 20.0, 200.0)",
      "border-width": 5,
      "font-size": 48
       }},

   {selector: "node[type='TF'][pearson<=0]", css: {
      "shape": "ellipse",
      "background-color": "mapData(pearson, -1.0, 0, green, white)",
      "width": "mapData(expression, 0.0, 1.0, 20.0, 200.0)",
      "height":"mapData(expression, 0.0, 1.0, 20.0, 200.0)"
       }},

   {selector: "node[type='tf']", css: {
       "shape": "rectangle"
       }},


   {selector: "node[score=1]", css: {
       "background-color": "#AAFFAA"
       }},

   {selector: "node[score=2]", css: {
       "background-color": "#DDDDFF"
       }},

   {selector: "node[score=3]", css: {
       "background-color": "#FF2222"
       }},

   {selector: 'edge', css: {
      "line-color": "gray",
      "target-arrow-shape": "triangle",
      "target-arrow-color": "rgb(0, 0, 0)",
      "curve-style": "bezier",
      "width": "3px"
      }},


   {selector: 'edge:selected', css: {
      "line-color": "red",
      "target-arrow-shape": "triangle",
      "target-arrow-color": "rgb(0, 0, 0)",
       width: "5px",
      'curve-style': 'bezier'
      }},

   {selector:"node:selected", css: {
       "text-valign":"center",
       "text-halign":"center",
       "border-color": "black",
       "content": "data(label)",
       "border-width": "3px",
       "overlay-opacity": 0.2,
       "overlay-color": "gray"
        }}

   ];

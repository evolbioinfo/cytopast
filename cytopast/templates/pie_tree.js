var cy = cytoscape({
  container: document.getElementById('cy'),

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'width': 'data(node_size)',
        'height': 'data(node_size)',
        'content': 'data(node_name)',
        'shape': 'data(shape)',
        'pie-size': '94%',
        'background-color': '#909090',
        'color': 'black',
        'text-opacity': 1,
        'text-valign': 'center',
        'text-halign': 'center',
        'font-size': 'data(node_fontsize)'
      })
    {% for (clazz, css) in clazz2css %}
    .selector(".{{clazz}}")
        .css({
        {{css}}
        })
    {% endfor %}
    .selector('edge')
      .css({
        'width': 'data(edge_size)',
        'font-size': 'data(edge_size)',
        'color': '#cccccc',
        'content': 'data(edge_name)',
        'curve-style': 'bezier',
        'target-arrow-shape': 'data(interaction)',
        'target-arrow-color': 'data(color)',
        'opacity': 0.8,
        'text-opacity': 1,
        'line-color': 'data(color)'
      })
    .selector(':selected')
      .css({
        'background-color': 'black',
        'line-color': 'black',
        'target-arrow-color': 'black',
        'source-arrow-color': 'black',
        'pie-size': '60%',
        'opacity': 1
      })
    .selector('.faded')
      .css({
        'opacity': 0.25,
        'text-opacity': 0
      }),

  elements: {{elements}},

  layout: {
    name: '{{layout}}',
    nodesep: 10,
    minLen: function( edge ){ return edge.data('minLen') | 1; }
  },

  ready: function(){
    window.cy = this;
  }
});

cy.on('mouseover', 'node', function(event) {
    var node = event.target;
    if (node.data('tooltip') !== undefined) {
        node.qtip({
            content: node.data('tooltip'),
            show: {
               event: event.type,
               ready: true
            },
            hide: {
               event: 'mouseout unfocus'
            },
            style: {
                classes: 'qtip-bootstrap',
            }
        }, event);
    }
});


function to_image(){
    document.getElementById("downloader").href = cy.jpg({ full: false, quality: 1.0, scale: 2}).replace(/^data:image\/[^;]/, 'data:application/octet-stream');
}

function fit() {
    cy.fit();
}



var slider = document.getElementById("myRange");
var output = document.getElementById("demo");
output.innerHTML = slider.value; // Display the default slider value

// Update the current slider value (each time you drag the slider handle)
var removed = cy.collection("[date>{{max_date}}]");

slider.oninput = function() {
    output.innerHTML = this.value;
    removed.restore();
    removed = cy.remove("[date>" + this.value + "]");
    var list = cy.$("");
    for (var i=0, node; node = list[i]; i++) {
        node.data('node_name', node.data('node_name_' + this.value));
        node.data('node_size', node.data('node_size_' + this.value));
        node.data('edge_name', node.data('edge_name_' + this.value));
        node.data('edge_size', node.data('edge_size_' + this.value));
        node.data('node_fontsize', node.data('node_fontsize_' + this.value));
    }
}
var cy = cytoscape({
  container: document.getElementById('cy'),

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'width': 'data(size)',
        'height': 'data(size)',
        'content': 'data(name)',
        'shape': 'data(shape)',
        'pie-size': '94%',
        'background-color': '#909090',
        'color': 'black',
        'text-opacity': 1,
        'text-valign': 'center',
        'text-halign': 'center',
        'font-size': 'data(fontsize)'
      })
    {% for (clazz, css) in clazz2css %}
    .selector(".{{clazz}}")
        .css({
        {{css}}
        })
    {% endfor %}
    .selector('edge')
      .css({
        'width': 'data(size)',
        'font-size': 'data(size)',
        'color': '#cccccc',
        'content': 'data(name)',
        'curve-style': 'bezier',
        'target-arrow-shape': 'data(interaction)',
        'target-arrow-color': 'data(color)',
        'opacity': 0.8,
        'text-opacity': 1,
        'content': 'data(name)',
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
    document.getElementById("downloader").download = "{{title}}.png";
    document.getElementById("downloader").href = cy.png({ full: true }).replace(/^data:image\/[^;]/, 'data:application/octet-stream');
}

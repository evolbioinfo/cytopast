var cy = cytoscape({
  container: document.getElementById('cy'),

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'width': 'data(size)',
        'height': 'data(size)',
        'content': 'data(name)',
        'shape': 'data(shape)',
        'pie-size': '100%',
        {% for (name, colour) in name2colour %}
            'pie-{{loop.index}}-background-color': "{{colour}}",
            'pie-{{loop.index}}-background-size': 'mapData({{name}}, 0, 100, 0, 100)',
        {% endfor %}
        'text-opacity': 0.5,
        'text-valign': 'center',
        'text-halign': 'center',
        'font-size': 'data(fontsize)'
      })
    .selector('edge')
      .css({
        'width': 'data(size)',
        'font-size': 'data(size)',
        'content': 'data(name)',
        'curve-style': 'bezier',
        'target-arrow-shape': 'data(interaction)',
        'opacity': 0.8,
        'content': 'data(name)',
      })
    .selector(':selected')
      .css({
        'background-color': 'black',
        'line-color': 'black',
        'target-arrow-color': 'black',
        'source-arrow-color': 'black',
        'pie-size': '50%',
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
             }
        }, event);
    }
});

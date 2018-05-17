var cy = cytoscape({
  container: document.getElementById('cy'),

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'width': 'data(node_size)',
        'height': 'data(node_size)',
        'content': 'data(node_name)',
        'shape': 'data(shape)',
        'pie-size': '95%',
        'background-color': '#909090',
        'color': '#383838',
        'text-opacity': 1,
        'text-valign': 'center',
        'text-halign': 'center',
        'font-size': 'data(node_fontsize)',
        'text-halign' : 'center',
        'text-valign' : 'center',
        'min-zoomed-font-size': 12,
        'text-background-color' : '#ffffff',
        'text-background-shape' : 'roundrectangle',
        'text-background-opacity': .3,
        'text-background-padding' : 1,
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
        'color': 'black',
        'content': 'data(edge_name)',
        'curve-style': 'bezier',
        'target-arrow-shape': 'data(interaction)',
        'target-arrow-color': 'data(color)',
        'opacity': 0.8,
        'text-opacity': 1,
        'line-color': 'data(color)',
        'text-background-color' : '#ffffff',
        'text-background-shape' : 'roundrectangle',
        'text-background-opacity': 1,
        'text-background-padding' : 4,
        'min-zoomed-font-size': 10,
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
    nodesep: 0, // the separation between adjacent nodes in the same rank
    edgeSep: 1, // the separation between adjacent edges in the same rank
    rankSep: 40, // the separation between adjacent ranks
    rankDir: 'TB', // 'TB' for top to bottom flow, 'LR' for left to right,
    ranker: 'longest-path', // Type of algorithm to assign a rank to each node in the input graph. Possible values: 'network-simplex', 'tight-tree' or 'longest-path'
    minLen: function( edge ){ return edge.data('minLen') | 0; }, // number of ranks to keep between the source and target of the edge

    // general layout options
    fit: true, // whether to fit to viewport
    padding: 1, // fit padding
    spacingFactor: undefined, // Applies a multiplicative factor (>0) to expand or compress the overall area that the nodes take up
    nodeDimensionsIncludeLabels: true, // whether labels should be included in determining the space used by a node
    animate: false, // whether to transition the node positions
    animateFilter: function( node, i ){ return true; }, // whether to animate specific nodes when animation is on; non-animated nodes immediately go to their final positions
    animationDuration: 500, // duration of animation in ms if enabled
    animationEasing: undefined, // easing of animation if enabled
    boundingBox: undefined, // constrain layout bounds; { x1, y1, x2, y2 } or { x1, y1, w, h }
    transform: function( node, pos ){ return pos; }, // a function that applies a transform to the final node position
    ready: function(){}, // on layoutready
    stop: function(){} // on layoutstop
  },

  ready: function(){
    window.cy = this;
  }
});

cy.filter(function(ele, i, eles) {
    return ele.isNode() && ele.data('tooltip') !== undefined;
} ).qtip({
    content: function(){
            var tooltip = this.data('tooltip');
            tooltip += '<br>id' + (this.data('node_meta') !== undefined ? 's: ': ': ') + this.data('node_root_id')
            + (this.data('node_meta') !== undefined ? ', ...': '');
            if (this.data('node_in_tips') !== undefined) {
                tooltip += '<br>tips inside: ' + this.data('node_in_tips');
            }
            if (this.data('node_all_tips') !== undefined) {
                tooltip += '<br>total tips in the subtree: ' + this.data('node_all_tips');
            }
            return tooltip;
        },
    show: {event: 'mouseover'},
    hide: {event: 'mouseout'},
    style: {
            classes: 'qtip-bootstrap',
    },
    position: {
        at: 'center center',
    }
});

function to_image(){
    document.getElementById("downloader").href = cy.jpg({ full: false, quality: 1.0, scale: 2}).replace(/^data:image\/[^;]/, 'data:application/octet-stream');
}

function fit() {
    cy.fit();
}

var slider = document.getElementById("myRange");
if (slider !== null) {
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
            if (node.data('node_in_tips_' + this.value) !== undefined) {
                node.data('node_in_tips', node.data('node_in_tips_' + this.value));
            }
            if (node.data('node_all_tips_' + this.value) !== undefined) {
                node.data('node_all_tips', node.data('node_all_tips_' + this.value));
            }
            node.data('node_size', node.data('node_size_' + this.value));
            node.data('edge_name', node.data('edge_name_' + this.value));
            node.data('edge_size', node.data('edge_size_' + this.value));
            if (node.data('edge_meta_' + this.value) !== undefined) {
                node.css('line-color', '#383838');
                node.css('target-arrow-color', '#383838');
                node.data('edge_meta', node.data('edge_meta_' + this.value))
            } else {
                node.css('line-color', '#909090');
                node.css('target-arrow-color', '#909090');
                node.removeData('edge_meta');
            }
            if (node.data('node_meta_' + this.value) !== undefined) {
                node.data('node_meta', node.data('node_meta_' + this.value));
            } else {
                node.removeData('node_meta');
            }
            node.data('node_fontsize', node.data('node_fontsize_' + this.value));
        }
    }
}
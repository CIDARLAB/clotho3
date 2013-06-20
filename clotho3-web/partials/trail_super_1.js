'use strict';

var tempModule = {
    "moduleName" : "widgetApp3",
    "moduleUrl" : "widgets/widget3-module.js"
};

//add it, get div id
var selectors = $clotho.api.bootstrap(tempModule)
    .then(function(selectors) {
        //jquery to move it to the right spot
        $(selectors[1]).appendTo("[myWidget]");
    });




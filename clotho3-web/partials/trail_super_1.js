'use strict';

var tempModule = {
    "moduleName" : "widgetApp3",
    "moduleUrl" : "widget/widgets/widget3-module.js"
};

//add it, get div id
var divId = $clotho.api.bootstrap(tempModule);

//jquery to move it to the right spot
$("[clotho-widget-uuid=" + divId + "]").appendTo("[myWidget]");

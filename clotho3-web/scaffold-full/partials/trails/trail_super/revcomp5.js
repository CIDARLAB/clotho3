'use strict';

var tempModule = {
    "moduleName" : "revcompApp",
    "moduleUrl" : "partials/trails/trail_super/widget/revcompApp.js"
};

//add it, append where desired
Application.bootstrap(tempModule)
    .then(function(selectors) {
        console.log(selectors);
        //jquery to move it to the right spot
        $(selectors[1]).appendTo("[insertWidgetHere]");
    });


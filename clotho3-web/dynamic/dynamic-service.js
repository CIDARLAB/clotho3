'use strict';

Application.Dynamic.service('Dynamic', ['PubSub', 'Clotho', function(PubSub, Clotho) {
    /*
     Purpose
     - This is a sort of test-ground for Clotho.show()... ability to dyanmically display content given some data asynchronously.

     Notes
     -

     Resources
     - http://www.bennadel.com/blog/2441-Nested-Views-Routing-And-Deep-Linking-With-AngularJS.htm
     - http://blog.freeside.co/post/42872615955/dynamic-templates-in-angular-routes
     */


    //adapted: http://stackoverflow.com/questions/14974271/can-you-change-a-path-without-reloading-the-controller-in-angularjs
    //future - $routeProvider - set reloadOnSearch = false
    var views = {};

    // these can be added asynchronously
    views.dyn = {"name" : "Dyn", "template" : "partials/dyn-partial.html"};
    views.dyn2 = {"name" : "Dyn2", "template" : "partials/dyn2-partial.html"};

    //e.g. usage
    //<div class="content" ng-include src="Dynamic.views.{{template_name}}"></div>

    return {
        views : views
    }

}]);


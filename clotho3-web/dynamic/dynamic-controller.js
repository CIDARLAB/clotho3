'use strict';

//note works in 1.1.4, can't pass function in 1.0.6

var dynamicCtrl = Application.Dynamic.controller('DynamicCtrl', ['$scope', 'Clotho', '$route', '$rootScope',  function($scope, Clotho, $route, $rootScope) {

    /*
     Purpose
     - This is a sort of test-ground for Clotho.show()... ability to dynamically display content given some data asynchronously.

     Notes
     -

     Resources
     - http://www.bennadel.com/blog/2441-Nested-Views-Routing-And-Deep-Linking-With-AngularJS.htm
     - http://blog.freeside.co/post/42872615955/dynamic-templates-in-angular-routes
     - http://stackoverflow.com/questions/12356185/angular-js-delaying-controller-initialization
     - http://stackoverflow.com/questions/13037051/routeprovider-injecting-controller-dependencies-depending-on-url
     */

    var resolved = $route.current.locals.resolve;
    $scope.model = resolved.model;
    $scope.template_url = resolved.template_url;


    $scope.editModel = function(uuid) {
        Clotho.edit(uuid);
    };

    $scope.showTemplate = function(uuid) {
        $scope.customTemplate = 'partials/' + uuid + '.html';
    };

    $scope.showDynamic = function(uuid) {
        $scope.dynTemplate = 'partials/' + uuid + '.html';
    };

    //bootstrapping new apps

    $scope.bootstrapWidgetOne = function() {
        var tempModule = {
            "moduleName" : "widgetApp",
            "moduleUrl" : "widget/widgets/widget-module.js"
        };
        console.log(Clotho.bootstrap(tempModule));
    };
    $scope.bootstrapWidgetTwo = function() {
        var tempModule = {
            "moduleName" : "widgetApp2",
            "moduleUrl" : "widget/widgets/widget2-module.js"
        };
        console.log(Clotho.bootstrap(tempModule));
    };

    //testing
    $scope.testData = "this is a new string";
    Clotho.listen('testMessage', function(data) {
        $scope.testData = data;
    }, '$scope.$id');


    //inserting new dependencies

    //testing
    //note - working
    $scope.addDirective = function() {
        $('newDirective').html('<simple-text></simple-text>');
        Application.Widgets.mixin('widget/dependencies/simpleText-directive.js', 'newDirective');
    };
    //note-working
    $scope.addController = function() {
        $('newController').html('<div ng-controller="SimpleCtrl">{{test}}<br />{{serviceText}}</div>');

        //note this controller requires a service to show we can do that too (dependency injection)

        //add service first, don't recompile (don't pass in element), just add dependency
        //Application.Widgets.mixin('widget/dependencies/simple-service.js');
        // then pass in the controller
        //Application.Widgets.mixin('widget/dependencies/simple-controller.js', 'newController');

        //or just pass it in all together
        Application.Widgets.mixin(['widget/dependencies/simple-controller.js', 'widget/dependencies/simple-service.js'], 'newController');

    };
    //note - not working
    $scope.addFilter = function() {
        $('newFilter').html("{{'tester' | capitalize}}");
        console.log("1");
        Application.Widgets.mixin('widget/dependencies/capitalize-filter.js', 'newFilter');
        console.log("2");
    };

}]);
/*dynamicCtrl.template = function(Clotho) {
    //return Clotho.get_url('show_template.html');
    return 'dynamic/dynamic-partial.html';
};*/

dynamicCtrl.template = 'dynamic/dynamic-partial.html';

dynamicCtrl.resolve = function(Clotho, $q, $timeout) {
    var resolved = {};
    var deferred = $q.defer();

    $q.all([
        resolved.model = Clotho.get('inst_first'),
        resolved.template_url = Clotho.get_url('show_template.html')

        //testing
        //$timeout(function () {console.log('timeout')}, 1500)

    ]).then(function(values) {
        deferred.resolve(resolved);
    });

    return deferred.promise;
};
'use strict';

//note works in 1.1.4, can't pass function in 1.0.6

var dynamicCtrl = Application.Dynamic.controller('DynamicCtrl', ['$scope', 'Clotho', '$route', '$rootScope', '$position', '$dialog', function($scope, Clotho, $route, $rootScope, $position, $dialog) {

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

    // testing show / display
    //files must use Application.Widgets.__ syntax
    $scope.showSimple = function(event) {

        var pos = $position.position(angular.element(event.target));

        Clotho.show({
            "template" : 'extensions/simple-template.html',
            "controller" : 'extensions/simple-controller.js',
            "dependencies" : [
                'extensions/simple-service.js'
            ],
            "styles" : {
                "background" : "#FF0000",
                "position" : "absolute",
                "top" : pos.top,
                "left" : pos.left
            }
        });
    };
    $scope.showEditor = function() {
        Clotho.show({
            "template" : 'extensions/editor-template.html',
            "args" : {
                "id" : "inst_first"
            }
            //"target" : ".editorCatcher"
        });
    };
    $scope.showEditorModal = function() {

        var editor_template = '<form sharable-editor name="sharableEditor" id="inst_first" class=" form-horizontal well" novalidate></form>';

        var dialog_opts = {
            backdrop: true,
            keyboard: true,
            backdropClick: true,
            template:  editor_template
        };
        var d = $dialog.dialog(dialog_opts);
        d.open();
    };

    //bootstrapping new apps

    $scope.bootstrapWidgetOne = function() {
        var tempModule = {
            "moduleName" : "widgetApp",
            "moduleUrl" : "widgets/widget-module.js"
        };
        console.log(Clotho.bootstrap(tempModule));
    };
    $scope.bootstrapWidgetTwo = function() {
        var tempModule = {
            "moduleName" : "widgetApp2",
            "moduleUrl" : "widgets/widget2-module.js"
        };
        console.log(Clotho.bootstrap(tempModule));
    };

    //testing
    $scope.testData = "this is a new string";
    Clotho.listen('testMessage', function(data) {
        $scope.testData = data;
    }, $scope);


    //inserting new dependencies

    //testing
    //note - working
    $scope.addDirective = function() {
        $('newDirective').html('<simple-text></simple-text>');
        Application.mixin('extensions/simpleText-directive.js', 'newDirective');
    };
    //note-working
    $scope.addController = function() {

        $('newController').html('<div ng-controller="SimpleCtrl">{{test}}<br />{{serviceText}}</div>');

        //note this controller requires a service to show we can do that too (dependency injection)

        //add service first, don't recompile (don't pass in element), just add dependency
        //Application.Widgets.mixin('extensions/simple-service.js');
        // then pass in the controller
        //Application.Widgets.mixin('extensions/simple-controller.js', 'newController');

        //or just pass it in all together
        Application.mixin(['extensions/simple-controller.js', 'extensions/simple-service.js'], 'newController');

    };
    //note - not working
    $scope.addFilter = function() {
        $('newFilter').html("{{'tester' | capitalize}}");
        Application.mixin('extensions/capitalize-filter.js', 'newFilter');
    };

}]);

dynamicCtrl.template = function(Clotho) {
    //return Clotho.get_url('show_template.html');
    return 'dynamic/dynamic-partial.html';
};

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
'use strict';

/**
 * The application file bootstraps the angular app by  initializing the main module and
 * creating namespaces and modules for controllers, filters, services, and directives.
 *
 * controllers should be defined in partials, not here
 */

// FUTURE -- do we want to actively cache templates in localStorage? pull them out on app loading?

var Application = Application || {};

Application.Chat = angular.module('clotho.chat', []);
Application.Dynamic = angular.module('clotho.dynamic', []);
Application.Editor = angular.module('clotho.editor', []);
Application.Search = angular.module('clotho.search', []);
Application.Trails = angular.module('clotho.trails', []);
Application.Primary = angular.module('clotho.primary', []);

Application.Widgets = angular.module('clotho.widgets', []); // lazy-loading module dependencies

Application.Foundation = angular.module('clotho.core', [])
    .run(['$rootScope', 'Clotho', function ($rootScope, Clotho) {
        //on first run, add API to $clotho object
        window.$clotho.api = Clotho;

        //extend scope with Clotho API
        $rootScope.Clotho = Clotho;

        /**
         @name $rootScope.$safeApply
         @description Particularly for 3rd party apps, when need to force digest or apply safely. Each app needs to insert this into its own run() clause

         You can run it like so:
         $scope.$safeApply(function() {
//this function is run once the apply process is running or has just finished
});

         An alternative is to use $timeout(function() {}), which will run after the previous $digest is complete. May be better to do it using timeout?
         */
        $rootScope.$safeApply = function(fn) {
            fn = fn || function() {};
            if($rootScope.$$phase) {
                //don't worry, the value gets set and AngularJS picks up on it...
                fn();
            }
            else {
                //this will fire to tell angularjs to notice that a change has happened
                //if it is outside of it's own behaviour...
                $rootScope.$apply(fn);
            }
        }
    }]);

angular.module('clothoRoot', ['clotho.core', 'clotho.chat', 'clotho.dynamic', 'clotho.editor', 'clotho.search', 'clotho.trails', 'clotho.widgets', 'clotho.primary']).
    config(['$routeProvider', function ($routeProvider) {
        $routeProvider.
            when('/', {
                templateUrl:'home/home-partial.html'
            }).
            when('/trails', {
                templateUrl:'trails/trail_browser-partial.html'
            }).
            when('/trails/:uuid', {
                templateUrl:'trails/trail-partial.html'
            }).
            when('/editor', {
                redirectTo:'/editor/inst_first'
            }).
            when('/editor/:uuid', {
                templateUrl:'editor/editor-partial.html'
                //todo - get this working, instead of doing it in the link of directive
                //resolve: []
            }).
            when('/chat', {
                templateUrl:'chat/chat-partial.html'
            }).
            when('/dynamic', {
                // note: can be function in 1.1.x, not 1.0.x (currently) - see: https://github.com/angular/angular.js/pull/1849/files
                //templateUrl:'dynamic/dynamic-partial.html',
                templateUrl: dynamicCtrl.template,
                resolve: {
                    resolve: dynamicCtrl.resolve
                },
                clotho : dynamicCtrl.clotho,
                custom : {
                    model : "inst_second"
                }
            }).
            otherwise({
                redirectTo:'/'
            });
    }])
    .run(['$rootScope', function($rootScope) {

    /**************
       CONFIG
    **************/

    //testing
    //$rootScope.$on('$destroy', console.log("\n\ndestroyed"));
    //todo - extend native $destroy() to unhook listeners (or emit event?)

    /**************
     TESTING
     **************/

    $rootScope.$on('$routeChangeError', function(event, current, previous, rejection) {
        console.log("Route Change Error: " + rejection);
    });

    /**************
     Functions
     **************/

}]);

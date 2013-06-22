'use strict';

/**
 * The application file bootstraps the angular app by  initializing the main module and
 * creating namespaces and modules for controllers, filters, services, and directives.
 *
 * controllers should be defined in partials, not here
 */

// FUTURE -- do we want to actively cache templates in localStorage? pull them out on app loading?

var Application = Application || {};

Application.Primary = angular.module('clotho.primary', ['clotho.core']);
Application.Interface = angular.module('clotho.interface', ['clotho.core']);
Application.Extensions = angular.module('clotho.extensions', ['clotho.core']);
Application.Widgets = angular.module('clotho.widgets', ['clotho.core']);

Application.Browser = angular.module('clotho.browser', ['clotho.core']);
Application.Chat = angular.module('clotho.chat', ['clotho.core']);
Application.Dynamic = angular.module('clotho.dynamic', ['clotho.core']);
Application.Editor = angular.module('clotho.editor', ['clotho.core']);
Application.Search = angular.module('clotho.search', ['clotho.core']);
Application.Trails = angular.module('clotho.trails', ['clotho.core']);

Application.Foundation = angular.module('clotho.core', [])
    .run(['$rootScope', 'Clotho', function ($rootScope, Clotho) {
        //on first run, add API to $clotho object
        window.$clotho.api = Clotho;

        //extend scope with Clotho API
        $rootScope.Clotho = Clotho;

        /**
         @name $rootScope.$safeApply
         @note Each app needs to insert this into its own run() clause
         @description Particularly for 3rd party apps, when need to force digest or apply safely.

         You can run it like so:
         $scope.$safeApply(function() {
//this function is run once the apply process is running or has just finished
});

         An alternative is to use $timeout(function() {}), which will run after the previous $digest is complete. However, this may cause UI flicker, as it will not run until the previous digest cycle is complete.
         */
        $rootScope.$safeApply = function(fn) {
            fn = fn || function() {};
            if($rootScope.$$phase) { fn(); }
            else { $rootScope.$apply(fn); }
        };
    }]);

angular.module('clothoPackage', ['clotho.browser', 'clotho.core', 'clotho.extensions', 'clotho.interface', 'clotho.primary', 'clotho.widgets', 'clotho.chat', 'clotho.dynamic', 'clotho.editor', 'clotho.search', 'clotho.trails']);

angular.module('clothoRoot', ['clothoPackage']).
    config(['$routeProvider', function ($routeProvider) {
        $routeProvider.
            when('/', {
                templateUrl:'home/home-partial.html'
            }).
            when('/trails', {
                templateUrl:'trails/trail_browser-partial.html'
            }).
            when('/trails/:uuid', {
                templateUrl:'trails/trail-partial.html',
                resolve : {
                    trail : function (Clotho, $q, $route, Trails) {
                        var deferred = $q.defer();
                        //todo - add timeout
                        Clotho.get($route.current.params.uuid).then(function(result) {
                            Trails.compile(result).then(function (compiled) {
                                deferred.resolve(compiled);
                            });
                        });
                        return deferred.promise;
                    }
                }
            }).
            when('/editor', {
                redirectTo:'/editor/inst_first'
            }).
            when('/editor/:id', {
                templateUrl:'editor/editor-partial.html'
                //todo - get this working, instead of doing it in the link of directive
                //resolve: []
            }).
            when('/browser', {
                templateUrl:'browser/browser-partial.html'
            }).
            when('/chat', {
                templateUrl:'chat/chat-partial.html'
            }).
            when('/dynamic', {
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
            when('/terminal', {
                templateUrl:'search/terminal-partial.html',
                resolve : {
                    //ctrl_dl: Application.mixin('search/terminal-controller.js')
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
        console.log("Route Change Error:");
        console.log(rejection);
    });

    /**************
     Functions
     **************/

}]);

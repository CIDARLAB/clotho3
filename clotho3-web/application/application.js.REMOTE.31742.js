'use strict';

/**
 * The application file bootstraps the angular app by  initializing the main module and
 * creating namespaces and modules for controllers, filters, services, and directives.
 *
 * controllers should be defined in partials, not here
 */

// FUTURE -- do we want to actively cache templates in localStorage? pull them out on app loading?

var Application = Application || {};

Application.Primary = angular.module('clotho.primary', []);
Application.Interface = angular.module('clotho.interface', []);
Application.Extensions = angular.module('clotho.extensions', []);
Application.Widgets = angular.module('clotho.widgets', []);

Application.Dna = angular.module('clotho.dna', []);

Application.Browser = angular.module('clotho.browser', []);
Application.Chat = angular.module('clotho.chat', []);
Application.Dynamic = angular.module('clotho.dynamic', []);
Application.Editor = angular.module('clotho.editor', []);
Application.Plasmid = angular.module('clotho.plasmid', []);
Application.Search = angular.module('clotho.search', []);
Application.Trails = angular.module('clotho.trails', []);
Application.Schemas = angular.module('clotho.schemas',[]);
Application.Functions = angular.module('clotho.functions',[]);

Application.Foundation = angular.module('clotho.setup', [])
    .run(['$rootScope', 'Clotho', function ($rootScope, Clotho) {
        //on first run, add API to $clotho object
        window.$clotho.api = Clotho;

        //extend scope with Clotho API
        $rootScope.Clotho = Clotho;
    }]);

angular.module('clotho.ng-additions', [])
    .run(['$rootScope', function($rootScope) {
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

        //angular function extensions
        var ext = {};

        ext.isEmpty = function(value) {
            return angular.isUndefined(value) || value === '' || value === null || value !== value;
        };

        ext.isScope = function(obj) {
            return obj && obj.$evalAsync && obj.$watch;
        };

        angular.extend(angular, ext);
    }]);

angular.module('clothoPackage', ['clotho.browser', 'clotho.setup', 'clotho.ng-additions', 'clotho.dna', 'clotho.extensions', 'clotho.interface', 'clotho.primary', 'clotho.widgets', 'clotho.chat', 'clotho.dynamic', 'clotho.editor', 'clotho.plasmid', 'clotho.search', 'clotho.trails', 'clotho.schemas', 'clotho.functions']);

angular.module('clothoRoot', ['clothoPackage']).
    config(['$routeProvider', function ($routeProvider) {
        $routeProvider.
            when('/', {
                templateUrl:'home/home-partial.html'
            }).
            when('/trails', {
                templateUrl:'trails/trail_browser-partial.html'
            }).
            when('/trails/:id', {
                templateUrl:'trails/trail-partial.html',
                resolve : {
                    trail : function (Clotho, $q, $route, Trails) {
                        var deferred = $q.defer();
                        //todo - add timeout
                        Clotho.get($route.current.params.id).then(function(result) {
                            Trails.compile(result).then(function (compiled) {
                                deferred.resolve(compiled);
                            });
                        });
                        return deferred.promise;
                    },
                    deps : function() {
                        return Application.mixin('https://www.youtube.com/iframe_api');
                    }
                }
            }).
            when('/editor', {
                templateUrl:'editor/editor-partial.html'
            }).
            when('/editor/:id', {
                templateUrl:'editor/editor-partial.html'
                //todo - get this working, instead of doing it in the link of directive
                //resolve: []
            }).
            when('/browser', {
                templateUrl:'browser/browser-partial.html'
            }).
            when('/plasmid', {
                templateUrl:'plasmid/plasmid-partial.html'
            }).
            when('/plasmid/:id', {
                templateUrl:'plasmid/plasmid-partial.html'
            }).
            when('/construction', {
                templateUrl:'dna/construction-partial.html'
            }).
            when('/chat', {
                templateUrl:'chat/chat-partial.html'
            }).
            when('/schemas', {
                templateUrl:'schemas/schemas-partial.html'
            }).
            when('/scripts', {
                templateUrl:'functions/functions-partial.html'
            }).
            when('/dynamic', {
                //templateUrl:'dynamic/dynamic-partial.html',
                templateUrl: dynamicCtrl.template,
                resolve: {
                    resolve : dynamicCtrl.resolve
                }
            }).
            when('/lazyload', {
                //templateUrl:'dynamic/dynamic-partial.html',
                templateUrl: '/testing/lazyLoad.html',
                resolve: {
                    deps : function() {
                        return Application.mixin('/partials/trails/trail_super/revcomp-filter.js')
                    }
                }
            }).
            when('/terminal', {
                templateUrl:'search/terminal-partial.html',
                resolve : {
                    deps : function() {
                        return Application.mixin('search/terminal-controller.js')
                    }
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
        console.log(event);
        console.log(current);
        console.log(previous);
        console.log(rejection);
    });

    /**************
     Functions
     **************/

}]);

'use strict';

/**
 * The application file bootstraps the angular app by  initializing the main module and
 * creating namespaces and modules for controllers, filters, services, and directives.
 *
 * controllers should be defined in partials, not here
 */

// FUTURE -- do we want to actively cache templates in localStorage? pull them out on app loading?

var Application = Application || {};

Application.Primary     = angular.module('clotho.primary', []);
Application.Interface   = angular.module('clotho.interface', ['ui.codemirror']);
Application.Extensions  = angular.module('clotho.extensions', []);
Application.Widgets     = angular.module('clotho.widgets', []);

Application.Dna         = angular.module('clotho.dna', []);

Application.Browser     = angular.module('clotho.browser', []);
Application.Chat        = angular.module('clotho.chat', []);
Application.Dynamic     = angular.module('clotho.dynamic', []);
Application.Editor      = angular.module('clotho.editor', []);
Application.Plasmid     = angular.module('clotho.plasmid', []);
Application.Search      = angular.module('clotho.search', []);
Application.Trails      = angular.module('clotho.trails', []);
Application.Schemas     = angular.module('clotho.schemas',[]);
Application.Functions   = angular.module('clotho.functions',[]);

Application.Foundation = angular.module('clotho.setup', [])
    .run(['$rootScope', 'Clotho', function ($rootScope, Clotho) {

		/**
		 @name $rootScope.$safeApply
		 @attribution https://github.com/yearofmoo/AngularJS-Scope.SafeApply
		 @note Just load this module and these functions will all be run
		 @description Particularly for 3rd party apps, when need to force digest or apply safely. Alternative is $timeout, but can cause flicker.

		 Can run as:
		 $scope.$safeApply();
		 $rootScope.$safeApply($scope);
		 $scope.$safeApply(<function>);
		 $rootScope.$safeApply($scope, <function>, <force>);
		 */
		$rootScope.$safeApply = function() {
			var $scope, fn, force = false;
			if(arguments.length == 1) {
				var arg = arguments[0];
				if(typeof arg == 'function') {
					fn = arg;
				}
				else {
					$scope = arg;
				}
			}
			else {
				$scope = arguments[0];
				fn = arguments[1];
				if(arguments.length == 3) {
					force = !!arguments[2];
				}
			}
			$scope = $scope || this;
			fn = fn || function() { };
			if(force || !$scope.$$phase) {
				$scope.$apply ? $scope.$apply(fn) : $scope.apply(fn);
			}
			else {
				fn();
			}
		};

		//extend scope with Clotho API so don't need to do this in each controller.
		$rootScope.Clotho = Clotho;

		//sort of init() function with server (assumes WebSocket is set up)
		Clotho.submit("clotho.run('clientSetup', [])");

    }]);

angular.module('clotho.ng-additions', [])
    .config([function() {


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
                templateUrl:'home/home-partial.html',
                title : 'Home'
            }).
            when('/about', {
                templateUrl:'about/about-partial.html',
                title : 'About'
            }).
            when('/team', {
                templateUrl:'about/team-partial.html',
                title : 'The Clotho Team'
            }).
            when('/trails', {
                templateUrl:'trails/trail_browser-partial.html',
                title : 'Trails'
            }).
            when('/trails/:id', {
                templateUrl:'trails/trail-partial.html',
                resolve : {
                    trail : function (Clotho, $q, $route, Trails) {
                        var deferred = $q.defer();
                        //todo - add timeout
                        Clotho.get($route.current.params.id).then(function(result) {
                            Trails.compile(result).then(function (compiled) {
                                $route.current.$$route.title = result.name;
                                deferred.resolve(compiled);
                            });
                        });
                        return deferred.promise;
                    },
                    deps : function() {
                        //return Application.mixin('https://www.youtube.com/player_api');
                        return Application.mixin('/lib/youtubeAPI.js');
                    }
                }
            }).
            when('/editor', {
                templateUrl:'editor/editor-partial.html',
                title : 'Editor',
				        resolve : {
					        deps : function($q) {
						        //todo - lazyload in directive, not all at once
						        //todo - test if can load library lazily
						        return Application.mixin('/lib/codemirror-3.19/lib/codemirror.js').then(function() {
							        $q.all([
								        Application.css('/lib/codemirror-3.19/lib/codemirror.css'),
								        Application.mixin('/lib/codemirror-3.19/mode/javascript/javascript.js'),
								        Application.mixin('/lib/codemirror-3.19/mode/python/python.js'),
								        Application.mixin('/lib/codemirror-3.19/mode/groovy/groovy.js'),
								        Application.mixin('/lib/codemirror-3.19/mode/clike/clike.js'),
								        Application.mixin('/lib/codemirror-3.19/mode/markdown/markdown.js')
							        ]);
						        })
					        }
				        }
            }).
            when('/editor/:id', {
                templateUrl:'editor/editor-partial.html'
                //todo - get this working, instead of doing it in the link of directive
                //use search param and set reloadOnSearch to false
                //resolve: []
            }).
            when('/browser', {
                templateUrl:'browser/browser-partial.html',
                title : 'Browser'
            }).
            when('/plasmid', {
                templateUrl:'plasmid/plasmid-partial.html'
            }).
            when('/plasmid/:id', {
                templateUrl:'plasmid/plasmid-partial.html'
            }).
            when('/dnaPlayground', {
                templateUrl:'dna/playground-partial.html',
                title : 'DNA Playground'
            }).
            when('/chat', {
                templateUrl:'chat/chat-partial.html'
            }).
            when('/schemas', {
                templateUrl:'schemas/schemas-partial.html'
            }).
            when('/dynamic', {
                //templateUrl:'dynamic/dynamic-partial.html',
                templateUrl: dynamicCtrl.template,
                resolve: {
                    resolve : dynamicCtrl.resolve
                }
            }).
            when('/terminal', {
                templateUrl:'search/terminal-partial.html',
                title : 'Terminal',
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
    .run(['$rootScope', '$route', '$window', function($rootScope, $route, $window) {

    /**************
       CONFIG
    **************/

    $rootScope.$on('$routeChangeSuccess', function() {
        var title = $route.current.$$route.title;
        //can't use interpolation in document title because ng-app is within body
        $window.document.title = 'Clotho' + (angular.isDefined(title) ? ' | ' + title : '');
    });


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

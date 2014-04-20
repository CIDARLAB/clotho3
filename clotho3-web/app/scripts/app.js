'use strict';

//full default angular stack
angular.module('clotho.fullPackage', [
	'clotho.foundation',    //core components as described above
	'clotho.commandbar',        //command bar
	'clotho.webapp',        //general web app components
	//additional webapp modules
	'clotho.editor',
	'clotho.interface',
	'clotho.trails'
]);

//web application set up
angular.module('clothoRoot', ['clotho.fullPackage'])
.config(function ($routeProvider) {
	$routeProvider
	.when('/', {
		templateUrl: 'views/home.html',
		controller: 'HomeCtrl',
		title : 'Home',
		hotkeys : [
			['h', 'Show Intro Modal', 'showHelp = true']
		]
	})
	.when('/about', {
	  templateUrl: 'views/about.html',
	  controller: 'AboutCtrl'
	})
	.when('/team', {
	  templateUrl: 'views/team.html',
	  controller: 'TeamCtrl'
	})
	.when('/browser', {
	  templateUrl: 'views/browser.html',
	  controller: 'BrowserCtrl'
	})
	.when('/edit', {
		redirectTo: '/editor'
	})
	.when('/editor/query/:queryTerm', {
		templateUrl: 'views/editor.html',
		controller: 'EditorCtrl',
		resolve : {
			queryResult : ['Clotho', '$route', '$q', function (Clotho, $route, $q) {
				var query = '';
				try {
					query = angular.fromJson($route.current.params.queryTerm);
				} catch (err) {
					console.warn('editor query malformed: ' + err)
				}

				console.log('querying for editor: ', query);

				if (!angular.isEmpty(query)) {
					return Clotho.query(query).then(function (data) {
						return data;
					});
				} else {
					return $q.when();
				}
			}],
			deps : ['codemirrorLoader', function(loader) {
				return loader.loadMain();
			}]
		}
	})
	.when('/editor/query', {
		redirectTo : '/editor'
	})
	.when('/editor/:id?', {
	  templateUrl: 'views/editor.html',
	  controller: 'EditorCtrl',
		resolve : {
			deps : ['codemirrorLoader', function(loader) {
				return loader.loadMain();
			}]
		}
	})
	.when('/trails', {
	  templateUrl: 'views/trails.html',
	  controller: 'TrailsCtrl'
	})
	.when('/trails/:id', {
		templateUrl: 'views/trail.html',
		controller: 'TrailCtrl',
		reloadOnSearch: false,
		resolve : {
			trail : ['Clotho', '$q', '$http', '$route', 'Trails', function (Clotho, $q, $http, $route, Trails) {
				var deferred = $q.defer();
				Clotho.get($route.current.params.id).then(function(result) {
					Trails.compile(result).then(function (compiled) {
						$route.current.$$route.title = result.name;
						deferred.resolve(compiled);
					});
				});
				return deferred.promise;
			}]
		},
		hotkeys : [
			['alt+left', 'Previous page of Trail', 'prev()'],
			['alt+right', 'Next page of Trail', 'next()']
		]
	})
	.when('/terminal', {
		templateUrl:'views/_command/terminal.html',
		title : 'Terminal',
		resolve : {
			deps : function() {
				return $clotho.extensions.mixin('_command/terminal.js')
			}
		}
	})


	//testing
	.when('/widgets', {
	  templateUrl: 'views/widgets.html',
	  controller: 'WidgetsCtrl'
	})


.when('/test/tokenizer', {
  templateUrl: 'views/test/tokenizer.html',
  controller: 'TestTokenizerCtrl'
})
	.otherwise({
		redirectTo:'/'
	});

})
.run(function($rootScope, $route, $window, $location, Clotho, CommandBar, hotkeys) {

	/****** Init ******/

	//sort of init() function with server, for easier scripting
	Clotho.submit("clotho.run('clientSetup', [])");

	/****** Config *****/

	$rootScope.$on('$routeChangeSuccess', function(event, current, previous) {
		var title = current.$$route.title;
		//can't use interpolation in document title because ng-app is within body
		$window.document.title = 'Clotho' + (angular.isDefined(title) ? ' | ' + title : '');
	});

	/******* Global Hotkeys ******/

	hotkeys.add('f', 'Focus Command Bar', function (event) {event.preventDefault(); CommandBar.focusInput(); } );
	hotkeys.add('a', 'Show Activity Log', function (event) {event.preventDefault(); CommandBar.showActivityLog(); } );
	hotkeys.add('g h', 'Go to Homepage', function () {$location.path('/') });
	hotkeys.add('g e', 'Go to Editor', function () { $location.path('/editor') });
	hotkeys.add('g t', 'Go to Trails', function () { $location.path('/trails') });

});

/*
 angular.element(document).ready(function() {
 angular.bootstrap(document, ['clothoRoot']);
 });
 */

'use strict';

//full default angular stack
angular.module('clotho.fullPackage', [
	'clotho.foundation',    //core components as described above
	'clotho.commandbar',        //command bar
	'clotho.webapp',        //general web app components
	//additional webapp modules
	'clotho.editor', 'clotho.interface', 'clotho.trails'
]);

//web application set up
angular.module('clothoRoot', ['clotho.fullPackage'])
.config(function ($routeProvider) {
	$routeProvider
	.when('/', {
		templateUrl: 'views/home.html',
		controller: 'HomeCtrl',
		title : 'Home'
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
	.when('/editor', {
	  templateUrl: 'views/editor.html',
	  controller: 'EditorCtrl',
		resolve : {
			deps : ['$q', function($q) {
				return $clotho.extensions.mixin('bower_components/codemirror/lib/codemirror.js').then(function() {
					$q.all([
						$clotho.extensions.css('bower_components/codemirror/lib/codemirror.css'),
						$clotho.extensions.mixin('bower_components/codemirror/mode/javascript/javascript.js'),
						$clotho.extensions.mixin('bower_components/codemirror/mode/css/css.js'),
						$clotho.extensions.mixin('bower_components/codemirror/mode/htmlmixed/htmlmixed.js'),
						$clotho.extensions.mixin('bower_components/codemirror/mode/python/python.js'),
						$clotho.extensions.mixin('bower_components/codemirror/mode/groovy/groovy.js'),
						$clotho.extensions.mixin('bower_components/codemirror/mode/clike/clike.js'),
						$clotho.extensions.mixin('bower_components/codemirror/mode/markdown/markdown.js')
					]);
				})
			}]
		}
	})
	.when('/trails', {
	  templateUrl: 'views/trails.html',
	  controller: 'TrailsCtrl'
	})
	.when('/trails/:id/:position?', {
		templateUrl: 'views/trail.html',
		controller: 'TrailCtrl',
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
		}
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


	.otherwise({
		redirectTo:'/'
	});
	//routes go in here

})
.run(function($rootScope, $route, $window, Clotho) {

	/****** Init ******/

	//sort of init() function with server (assumes WebSocket is set up)
	Clotho.submit("clotho.run('clientSetup', [])");

	/****** Config *****/

	$rootScope.$on('$routeChangeSuccess', function(event, current, previous) {
		var title = current.$$route.title;
		//can't use interpolation in document title because ng-app is within body
		$window.document.title = 'Clotho' + (angular.isDefined(title) ? ' | ' + title : '');
	});
});

angular.element(document).ready(function() {
	angular.bootstrap(document, ['clothoRoot']);
});
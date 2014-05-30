'use strict';

//full default angular stack
angular.module('clotho.fullPackage', [
	'clotho.foundation',    //core components as described above
	'clotho.commandbar',        //command bar
	'clotho.dna',         //in here temporarily until server can handle all this
	'clotho.interface',
	'clotho.trails',
	'clotho.construction',
	'ngSanitize',
	'ngRoute'
]);

//web application set up
angular.module('clothoRoot', ['clotho.fullPackage'])
	.config(function ($routeProvider, $locationProvider) {

		/*
		 simulate legacy browser not supporting pushstate (add $provide to DI clause above)
		 $provide.decorator('$sniffer', function($delegate) {
		 $delegate.history = false;
		 return $delegate;
		 });
		 */

		$locationProvider
			.html5Mode(false)
			.hashPrefix('!');

		$routeProvider
			.when('/', {
				templateUrl: 'views/trail-splash.html',
				controller: 'TrailSplashCtrl',
				title: 'Home',
				hotkeys: [
					['h', 'Show Intro Modal', 'showHelp = !showHelp']
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
			.when('/trail', {
				templateUrl: 'views/trail.html',
				controller: 'TrailCtrl',
				reloadOnSearch: false,
				resolve: {
					trail: ['Clotho', '$q', '$http', '$route', 'Trails', function (Clotho, $q, $http, $route, Trails) {
						if (angular.isUndefined($route.current.params.id)) {
							return $q.when(null);
						}
						var deferred = $q.defer();
						Clotho.get($route.current.params.id).then(function (result) {
							Trails.compile(result).then(function (compiled) {
								$route.current.$$route.title = result.name;
								deferred.resolve(compiled);
							});
						}, function trailGetRejection() {
							$http.get('models/org.clothocad.trails.LearningClotho.json').then(function (data) {
								Trails.compile(data.data).then(function (compiled) {
									deferred.resolve(compiled);
								});
							});
						});
						return deferred.promise;
					}]
				},
				hotkeys: [
					['alt+left', 'Previous page of Trail', 'prev()'],
					['alt+right', 'Next page of Trail', 'next()']
				]
			})
			.otherwise({
				redirectTo: '/'
			});

	})
	.run(function ($rootScope, $route, $window, interfaceConfig) {

		/****** Config *****/

		$rootScope.$on('$routeChangeSuccess', function (event, current, previous) {
			var title = angular.isDefined(current.$$route) ? current.$$route.title : null;
			//can't use interpolation in document title because ng-app is within body
			$window.document.title = 'Clotho Trails' + (angular.isDefined(title) ? ' | ' + title : '');
		});

		$rootScope.interfaceConfig = interfaceConfig;

	});

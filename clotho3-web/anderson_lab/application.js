'use strict';

angular.module('andersonLab', ['ngRoute', 'clotho.foundation'])
	.config(function ($routeProvider) {
		$routeProvider
			.when('/', {
				templateUrl: 'partials/home.html'
			})
			.when('/people', {
				templateUrl: 'partials/people.html',
				resolve: {
					faculty: ['$q', '$route', '$http', 'Clotho', function ($q, $route, $http, Clotho) {

						//testing
						//return $http.get("assets/labroster.json").then(function(data) { return data.data });


						var deferred = $q.defer();
						Clotho.query({"schema" : "org.clothocad.model.LabPerson", "lab" : "Anderson"}).then(function (result) {
							deferred.resolve(result);
						});
						return deferred.promise;

					}]
				}
			})
			.when('/people/:id', {
				templateUrl: 'partials/person.html',
				resolve: {
					individual: ['$q', '$route', '$http', 'Clotho', function ($q, $route, $http, Clotho) {
						//testing
						//return $http.get("assets/labroster.json").then(function(data) { return data.data[data.data.length-1] });

						var deferred = $q.defer();

						Clotho.query({
							"schema" : "org.clothocad.model.LabPerson",
							"lab" : "Anderson",
							"name" : $route.current.params.id
						}).then(function (result) {
							deferred.resolve(result[0]);
						});

						return deferred.promise;

					}]
				}
			})
			.when('/research', {
				templateUrl: 'partials/research.html'
			})
			.when('/software', {
				templateUrl: 'partials/software.html'
			})
			.otherwise({
				redirectTo: '/'
			})
	})
	.run(function ($timeout) {
		//note - need to load and run these scripts after bootstrapping so that the elements are present
		$timeout(function () {
			$script(["scripts/slides.min.jquery.js", "scripts/jquery.columnizer.min.js"], "Busse-components", function() {
				$script("scripts/cl02.js", "clotho");
			});
		});
	})
	.controller('PeopleCtrl', function ($scope, $route, $location) {

		var faculty = $route.current.locals.faculty;

		$scope.faculty = {};
		$scope.faculty.management = [];
		$scope.faculty.current = [];
		$scope.faculty.past = [];
		
		//stuff below is to parse faculty into bins so don't need to use filters on the page, since won't be changing past controller instantiation
		/*
		//end up with this format
		$scope.faculty.current = [
			{ "name" : "postdocs"          ,      "members"  :    [] },
			{ "name" : "grads"             ,      "members"  :    [] },
			{ "name" : "programmers"       ,      "members"  :    [] },
			{ "name" : "undergrads"        ,      "members"  :    [] },
			{ "name" : "igem"              ,      "members"  :    [] },
			{ "name" : "collaborators"     ,      "members"  :    [] }
		];
		*/
		var facultyBin = {};
		facultyBin.current = {
			"postdocs": [],
			"grads": [],
			"programmers": [],
			"undergrads": [],
			"igem": [],
			"collaborators": []
		};
		facultyBin.past = angular.copy(facultyBin.current);

		angular.forEach(faculty, function(fac) {
			if (!!fac.isManagement) {
				$scope.faculty.management.push(fac);
				return;
			}

			var time = fac.current ? 'current' : 'past';

			switch (fac.role) {
				case "Postdoc" : {
					facultyBin[time].postdocs.push(fac);
					return;
				}
				case "Programmer" : {
					facultyBin[time].programmers.push(fac);
					return;
				}
				case "Graduate Student" : {
					facultyBin[time].grads.push(fac);
					return;
				}
				case "Undergraduate Researcher" : {
					facultyBin[time].undergrads.push(fac);
					return;
				}
				default : {
					if ((/iGEM/ig).test(fac.role)) {
						facultyBin[time].igem.push(fac);
						return;
					} else {
						//future - handle collaborators
						console.log(fac);
					}
				}
			}
		});

		var nameMap = [
			{ "short" : "postdocs"      ,   "screen" : "Postdocs"                     },
			{ "short" : "grads"         ,   "screen" : "Graduate Students"            },
			{ "short" : "programmers"   ,   "screen" : "Programmers"                  },
			{ "short" : "undergrads"    ,   "screen" : "Undergraduate Researchers"    },
			{ "short" : "igem"          ,   "screen" : "iGEM Participants"            },
			{ "short" : "collaborators" ,   "screen" : "Collaborators"                }
		];

		angular.forEach(facultyBin, function (items, time) {
			angular.forEach(nameMap, function (nameObj) {
				$scope.faculty[time].push({
					"name" : nameObj.screen ,
					"members"  : items[nameObj.short]
				})
			})
		});

	})
	.controller('PersonCtrl', function ($scope, $route) {

		//inherit from routeProvider.resolve(), promise fulfilled before render
		$scope.indiv = $route.current.locals.individual;

		console.log($scope.indiv);

		$scope.id = $route.current.params.id;

	})
	.directive('slideshow', function () {

		//todo ??

	});
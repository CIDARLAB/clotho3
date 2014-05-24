'use strict';

angular.module('clotho.webapp')
  .controller('TestContstructiontrailCtrl', function ($scope, $route, $timeout, Clotho, Trails, $location) {

		//inherited from $routeProvider.resolve clause in application.js
		$scope.id = $route.current.params.id;
		$scope.trail = $route.current.locals.trail;

		//initial position check
		if ($route.current.params.position) {
			$scope.current = $route.current.params.position;
			$scope.currentPage = Trails.extractPage($scope.trail, $scope.current);
		}

		$scope.activate = Trails.activate;

		$scope.favorite = function() {
			//future - once this concept exists, do initial check etc.
			Trails.favorite($scope.id);
		};

		$scope.share = function() {
			Trails.share($scope.id)
		};

		$scope.next = function() {
			$scope.activate(Trails.calcNextPage($scope.trail, $scope.current));
		};

		$scope.prev = function() {
			$scope.activate(Trails.calcPrevPage($scope.trail, $scope.current));
		};

		$scope.mapIcon = Trails.mapIcon;

		//listen for position changes
		$scope.$on('$routeUpdate', function(scope, next, current) {
			//don't activate already active one
			if ($scope.current == next.params.position) return;

			$scope.current = next.params.position || null;
			$scope.currentPage = Trails.extractPage($scope.trail, $scope.current);

			//if page doesn't exist, re-route to start
			if (angular.isEmpty($scope.currentPage)) {
				$scope.activate('0-0');
			}
		});

		$scope.$on('$destroy', function (scope, next, current) {
			$location.search('position', null).replace();
			$location.search('id', null).replace();
		});

		//kickoff
		$scope.activate($scope.current);


	});


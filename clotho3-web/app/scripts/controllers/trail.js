angular.module('clotho.webapp').controller('TrailCtrl', function($scope, $route, $timeout, Clotho, Trails, $keypress, $location) {

	//inherited from $routeProvider.resolve clause in application.js
	$scope.id = $route.current.params.id;
	$scope.trail = $route.current.locals.trail;

	$scope.activate = function(indices) {
		//if passed nothing
		if (!indices || !angular.isString(indices)) return;

		//don't activate already active one
		if ($scope.current == indices) return;

		$scope.current = indices;
		$scope.currentPage = Trails.extractPage($scope.trail, $scope.current);

		//if page doesn't exist, re-route to start
		if (angular.isEmpty($scope.currentPage)) {
			$scope.activate('0-0');
		}

		$location.search('position', $scope.current);
	};

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

	$keypress.on('keydown', {'alt-right' : 'next()', 'alt-left' : 'prev()'}, $scope);

	//listen for position changes
	$scope.$on('$routeUpdate', function(scope, next, current) {
		$scope.activate(next.params.position);
	});

	$scope.$on('$destroy', function (scope, next, current) {
		$location.search('position', null).replace();
	});

	//kickoff
	$scope.activate($route.current.params.position || '0-0');
});

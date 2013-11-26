angular.module('clotho.webapp').controller('TrailCtrl', function($scope, $route, $timeout, Clotho, Trails,$keypress) {

	//inherited from $routeProvider.resolve clause in application.js
	$scope.id = $route.current.params.id;
	$scope.trail = $route.current.locals.trail;


	$scope.activate = function(indices) {

		console.log(indices);

		//if passed nothing
		if (!indices || !angular.isString(indices)) return;

		//don't activate already active one
		if ($scope.current == indices) return;

		$scope.current = indices;
		$scope.currentPage = Trails.extractPage($scope.trail, $scope.current);
		//$scope.content = Trails.loadPage(page, $scope);
	};

	$scope.home = function() {
		$scope.content = $scope.trail.description;
		$scope.current = undefined;
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

	$scope.base64icon = base64icon;

	//kickoff
	$scope.activate('0-0');
});

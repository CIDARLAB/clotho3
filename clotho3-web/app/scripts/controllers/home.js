angular.module('clotho.webapp').controller('HomeCtrl', function ($scope, $location) {
	$scope.modalContent = '<p>Welcome to Clotho!</p>' +
		'<p>Clotho is a platform for automating your genetic engineering projects. Learn how to use Clotho by starting the trail below!</p>';

	$scope.enterClotho = function() {
		$location.path('/trails/org.clothocad.trails.LearningClotho');
	};

	$scope.enterEugene = function() {
		$location.path('/trails/org.clothocad.trails.EugeneCADIntro');
	};

	$scope.enterRaven = function() {
		$location.path('/trails/org.clothocad.trails.RavenCADIntro');
	};
});

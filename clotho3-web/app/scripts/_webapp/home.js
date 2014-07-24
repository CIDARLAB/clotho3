angular.module('clotho.webapp').controller('HomeCtrl', function ($scope, Clotho) {
	$scope.modalContent = '<p>Welcome to Clotho!</p>' +
		'<p>Clotho is a platform for automating your genetic engineering projects. Learn how to use Clotho by starting the trail below!</p>';

	$scope.enterClotho = function() {
		Clotho.startTrail('org.clothocad.trails.LearningClotho');
	};

	$scope.enterEugene = function() {
		Clotho.startTrail('org.clothocad.trails.EugeneCADIntro');
	};

	$scope.enterRaven = function() {
		Clotho.startTrail('org.clothocad.trails.RavenCADIntro');
	};
});

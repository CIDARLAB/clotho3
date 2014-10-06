angular.module('clotho.webapp').controller('HomeCtrl', function ($scope, Clotho, hotkeys) {
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

  //register listener in controller as callback instead of expression so doesn't get triggered in inputs
  hotkeys.bindTo($scope)
    .add({
      combo: 'h',
      description: 'Show Intro Modal',
      callback: function() {
        $scope.showHelp = !$scope.showHelp;
      }
    })
});

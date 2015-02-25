angular.module('clotho.webapp').controller('HomeCtrl', function ($scope, Clotho, hotkeys) {
    $scope.modalContent = '<p>Welcome to Clotho!</p>' +
        '<p>Clotho is a platform for automating your genetic engineering projects. Learn how to use Clotho by starting the trail below!</p>';

    $scope.downloadClotho = function(){
        window.open("https://github.com/CIDARLAB/clotho3", '_blank');
    };

    $scope.browseApps = function(){
        window.open("https://github.com/CIDARLAB/clotho3apps", '_blank');
    };

    $scope.enterClotho = function() {
        Clotho.startTrail('org.clothocad.trails.LearningClotho');
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

'use strict';

angular.module('clotho.trails')
  .controller('TestTrailBrowserCtrl', function ($scope, Clotho) {
    Clotho.query({"schema" : "org.clothocad.model.Trail"}).then(function (result) {
	    $scope.trails = result;
    });

		$scope.startTrail = Clotho.startTrail;
  });

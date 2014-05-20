'use strict';

angular.module('clotho.trails')
  .controller('TestTrailBrowserCtrl', function ($scope, Clotho) {
    Clotho.query({"schema" : "Trail"}).then(function (result) {
	    $scope.trails = result;
    });

		$scope.startTrail = Clotho.startTrail;
  });

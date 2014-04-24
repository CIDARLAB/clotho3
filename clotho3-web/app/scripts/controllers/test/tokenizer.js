'use strict';

angular.module('clotho.webapp')
  .controller('TestTokenizerCtrl', function ($scope) {

		$scope.myModel = [];
		
		$scope.$watch('myModel', function (newval, oldval) {

			//todo - model is updating, but this isn't firing for some reason past first call

			console.log('CTRL model update', newval, oldval);
		});

});

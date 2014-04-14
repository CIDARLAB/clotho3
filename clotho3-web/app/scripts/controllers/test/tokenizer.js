'use strict';

angular.module('clotho.webapp')
  .controller('TestTokenizerCtrl', function ($scope) {

		$scope.myModel = [];
		
		$scope.$watch('myModel', function (newval, oldval) {
			console.log(newval, oldval);
		});

});

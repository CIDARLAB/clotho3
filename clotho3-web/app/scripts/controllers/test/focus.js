'use strict';

angular.module('clotho.webapp')
  .controller('TestFocusCtrl', function ($scope, CommandBar, $focus, $parse) {

    $scope.setInput = function () {
	    CommandBar.setInput('revcomp aacc');
    };

		$scope.type = function () {
			$focus.typeOutSearch('revcomp aacc');
		};
  });

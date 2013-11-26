'use strict';

angular.module('clotho.webapp').controller('AboutCtrl', function ($scope) {
	$scope.joinMailingList = function(from) {
		alert('this doesnt do anything right now');
		$scope.emailToAdd = null;
	}
});

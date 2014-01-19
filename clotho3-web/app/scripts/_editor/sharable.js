angular.module('clotho.editor').controller('Editor_SharableCtrl', function($scope, $compile, $q, Clotho) {

	$scope.addNewField = function () {
		if ($scope.newFieldKey && $scope.newFieldVal) {
			$scope.sharable[$scope.newFieldKey] = $scope.newFieldVal;
			$scope.newFieldKey = null;
			$scope.newFieldVal = null;
			console.log($scope.sharable);
		}
	};

});
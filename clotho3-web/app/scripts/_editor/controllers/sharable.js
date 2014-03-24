angular.module('clotho.editor').controller('Editor_SharableCtrl', function($scope) {

	$scope.addNewField = function () {
		if ($scope.newFieldKey && $scope.newFieldVal) {
			$scope.sharable[$scope.newFieldKey] = $scope.newFieldVal;
			$scope.newFieldKey = null;
			$scope.newFieldVal = null;
		}
	};

});
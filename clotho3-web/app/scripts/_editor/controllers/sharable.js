angular.module('clotho.editor').controller('Editor_SharableCtrl', function($scope) {

	$scope.$watch('sharable', function (newval, oldval) {
		console.log('sharable updated ', newval);
	});

	$scope.addNewField = function () {
		if ($scope.newFieldKey && $scope.newFieldVal) {
			$scope.sharable[$scope.newFieldKey] = $scope.newFieldVal;
			$scope.newFieldKey = null;
			$scope.newFieldVal = null;
		}
	};

});
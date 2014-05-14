angular.module('clotho.editor').controller('Editor_SharableCtrl', function($scope) {

	$scope.addNewField = function () {
		if ($scope.newFieldKey && $scope.newFieldVal) {
			$scope.sharable[$scope.newFieldKey] = $scope.newFieldVal;
			$scope.newFieldKey = null;
			$scope.newFieldVal = null;
		}
	};

	//check if the new field key will overwrite another
	//directive may make more sense...
	$scope.checkNewFieldName = function () {
		if (angular.isDefined($scope.sharable[$scope.newFieldKey])) {
			$scope.newField.$setInvalid();
		}
	}

})

/**
 * checks on changes to ngModel (newFieldKey) to see if is in sharable already (as attr to this)
 *
 * @example (sharable editor)
 *  <input type="text" ng-model="newFieldKey" sharable-check-field-exists="sharable">
 */
.directive('sharableCheckFieldExists', function ($parse) {
		return {
			restrict : 'A',
			require : 'ngModel',
			link : function (scope, element, attrs, ngModelCtrl) {
				ngModelCtrl.$parsers.push(function (viewValue) {
					if (viewValue) {
						var sharable = $parse(attrs['sharableCheckFieldExists'])(scope);
						var isDup = angular.isDefined(sharable[viewValue]);
						ngModelCtrl.$setValidity('duplicateField', !isDup);
						return viewValue;
					}
				});
			}
		}
	});
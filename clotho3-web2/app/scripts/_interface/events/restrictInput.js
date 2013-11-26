//alternative to ng-pattern, which will allow all input but if model not acceptable then not set -- this will not allow invalid values to propagate to model, or be visible in the input field
//todo - update
angular.module('clotho.interface').directive('restrictInput', function() {
	return {
		restrict: 'A',
		require: 'ngModel',
		link: function(scope, iElement, iAttrs, ngModel) {
			var fn = function(input) {
				//in function so remains dynamic
				var regexp = scope.$eval(iAttrs.restrictInput);

				var transformedInput = input.replace(regexp, '');
				if (transformedInput != input) {
					ngModel.$setViewValue(transformedInput);
					ngModel.$render();
				}
				return transformedInput;
			};

			ngModel.$parsers.unshift( fn );
		}
	}
});
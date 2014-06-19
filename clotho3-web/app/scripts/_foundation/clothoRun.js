'use strict';

/**
 *
 * @ngdoc directive
 * @name clotho-run
 *
 * @description
 * This is meant for read-only modifications to the model - Can't use ngModel $formatters with promises (angular-1.1.5). Calling $setViewValue affects model and propagates, which is often undesired (e.g. if revcomp a sequence for display, don't want to change model).
 *
 * note : form of ngModel and functions of higher arity --
 * it is assumed that if ngModel is an array, it is in the correct format as if being passed to a function.apply()....
 * if it is not an array, it is assumed the function is of single arity, and ngModel is wrapped in an array as the only value
 * (i.e. [ngModel.$modelValue] )
 *
 * note : updating parent scope
 * There are two ways to do this.
 * (1) Because an isolate scope is created, you could pass in $parent.<model> which should do normal angular binding.
 * (2) If you want to update the model only with the run function, add the attribute tag clotho-run-update-model="true" ... this will update the model with the result of the run function once it is complete.
 *
 * //todo - allow passage of specific model to update - clotho-run-update=""  - using $parse().assign
 *
 * @example
 *
 * Assuming `$scope.myModel = 'HEY THERE';`
 *
 * `<p clotho-run="lowercase" ng-model="myModel" clotho-run-update="updateMe"></p>`
 *
 * will output
 *
 * `<p clotho-run="lowercase" ng-model="myModel">hey there</p>`
 *
 * and update the model updateMe to 'hey there' as well
 *
 */
angular.module('clotho.clothoDirectives')
.directive('clothoRun', function(Clotho, $parse) {

	var inputsVal = {input: true, textarea : true, select: true};

	return {
		restrict : 'A',
		require : 'ngModel',
		scope : true,
		link: function (scope, element, attrs, ngModel) {

			//config
			var useVal = !!inputsVal[angular.lowercase(element[0].nodeName)];
			if (useVal) {
				//avoid flicker
				ngModel.$render = angular.noop;
			}

			var updateModel = $parse(attrs.clothoRunUpdate);

			//command, args
			scope.$watch(function() {
				return attrs.clothoRun
			}, function(newval, oldval) {
				if (!!newval) {
					runFunction(ngModel.$modelValue);
				}
			});

			//model changes
			scope.$watch(function() {
				return ngModel.$modelValue
			}, function(newval, oldval) {
				runFunction(newval);
			});

			function updateParentModel (newModel) {
				angular.isDefined(updateModel.assign) && updateModel.assign(scope, newModel);
			}

			function updateElement (newval) {
				var method = useVal ? 'val' : 'text';
				element[method](newval);
			}

			//form array out of arguments if not an array
			function arrayify (input) {
				return angular.isArray(input) ? input : [input];
			}

			var runFunction = function(input) {
				input = arrayify(input);

				return Clotho.run(attrs.clothoRun, input).then(function(result) {
					updateParentModel(result);
					updateElement(result);
				});
			};
		}
	}
});
//adds Clotho and $safeApply to $rootScope
angular.module('clotho.setup', [])
	.run(function ($rootScope, Clotho) {

		//extend scope with Clotho API so don't need to do this in each controller.
		$rootScope.Clotho = Clotho;

		/**
		 @name $rootScope.$safeApply
		 @adaptedFrom: https://github.com/yearofmoo/AngularJS-Scope.SafeApply
		 @description Particularly for 3rd party apps, when need to force digest or apply safely. Alternative is $timeout, but can cause flicker.

		 @usage
		 $scope.$safeApply();
		 $rootScope.$safeApply($scope);
		 $scope.$safeApply(<function>);
		 $rootScope.$safeApply($scope, <function>, <force>);
		 */
		$rootScope.$safeApply = function() {
			var $scope, fn, force = false;
			if(arguments.length == 1) {
				var arg = arguments[0];
				if(typeof arg == 'function') {
					fn = arg;
				}
				else {
					$scope = arg;
				}
			}
			else {
				$scope = arguments[0];
				fn = arguments[1];
				if(arguments.length == 3) {
					force = !!arguments[2];
				}
			}
			$scope = $scope || this;
			fn = fn || function() { };
			if(force || !$scope.$$phase) {
				$scope.$apply ? $scope.$apply(fn) : $scope.apply(fn);
			}
			else {
				fn();
			}
		};
	});
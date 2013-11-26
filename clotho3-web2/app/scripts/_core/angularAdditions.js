//angular components needed, add a couple methods to prototype
angular.module('clotho.angularAdditions', [
		'ngCookies', 'ngSanitize', 'ngRoute'
	])
	.config(function() {
		//angular function extensions
		var ext = {};
		ext.isEmpty = function(value) {
			return angular.isUndefined(value) || value === '' || value === null || value !== value;
		};
		ext.isScope = function(obj) {
			return obj && obj.$evalAsync && obj.$watch;
		};
		angular.extend(angular, ext);
	})
	.run(function($rootScope) {
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
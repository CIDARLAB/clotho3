//angular components needed, add a couple methods to prototype
angular.module('clotho.angularAdditions', [])
	.config(function() {
		//angular function extensions
		var ext = {};
		/**
		 * @name angular.isEmpty
		 * @description Determines whether object is empty, delegating to lodash if present
		 * @param {*} value
		 * @returns {*|boolean|Boolean}
		 */
		ext.isEmpty = function(value) {
			return (angular.isDefined(_) && angular.isObject(value) && _.isEmpty(value)) || (angular.isUndefined(value) || value === '' || value === null || value !== value);
		};
		/**
		 * @name angular.isScope
		 * @description Determines whether an object is an angular $scope
		 * @param {*} obj
		 * @returns {*|$evalAsync|$watch|Function}
		 */
		ext.isScope = function(obj) {
			return obj && obj.$evalAsync && obj.$watch;
		};
		/**
		 * @name angular.once
		 * @description Creates a function that is restricted to execute `func` once. Repeat calls to the function will return the value of the first call. The `func` is executed with the `this` binding of the created function.
		 * @param {Function} func The function to restrict.
		 * @returns {Function} Returns the new restricted function.
		 */
		ext.once = function (func) {
			var ran,
				result;

			if (!angular.isFunction(func)) {
				throw new TypeError;
			}
			return function() {
				if (ran) {
					return result;
				}
				ran = true;
				result = func.apply(this, arguments);

				func = null;
				return result;
			};
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
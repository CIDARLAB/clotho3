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
		ext.isEmpty = function isEmpty(value) {
			//can't check if _ defined with angular because throws ReferenceError
			try {
				return _.isEmpty(value)
			} catch (err) {
				if (angular.isNumber(value)) {
					return value === value;
				}
				else if (angular.isObject(value)) {
					if (value.length === 0) {
						return true;
					} else {
						for (var key in value) {
							if (value.hasOwnProperty(key)) {
								return false;
							}
						}
						return true;
					}
				}
				else if (angular.isUndefined(value) || value === null || value.length == 0) {
					return true;
				}
				return false;
			}
		};
		/**
		 * @name angular.isScope
		 * @description Determines whether an object is an angular $scope
		 * @param {*} obj
		 * @returns {Boolean}
		 */
		ext.isScope = function isScope(obj) {
			return !!obj && angular.isFunction(obj.$evalAsync) && angular.isFunction(obj.$watch);
		};
		/**
		 * @name angular.once
		 * @description Creates a function that is restricted to execute `func` once. Repeat calls to the function will return the value of the first call. The `func` is executed with the `this` binding of the created function.
		 * @param {Function} func The function to restrict.
		 * @returns {Function} Returns the new restricted function.
		 */
		ext.once = function once(func) {
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

		/**
		 * Removes all elements from an array that the callback returns truey for and returns an array of removed elements.
		 * The callback is invoked with three arguments; (value, index, array).
		 * @param {Array} array The array to modify
		 * @param {Function} callback The function called per iteration
		 * @param {*} context `this`
		 * @returns {Array} Returns a new array of removed elements.
		 */
		ext.remove = function remove(array, callback, context) {
			var index = -1,
				length = array ? array.length : 0,
				result = [];

			context = context || null;

			while (++index < length) {
				var value = array[index];
				if (angular.isUndefined(value) || callback.call(context, value, index, array)) {
					angular.isDefined(value) && result.push(value);
					array.splice(index--, 1);
					length--;
				}
			}
			return result;
		};

		ext.map = function map(obj, iterator, context) {
			var results = [];
			angular.forEach(obj, function(value, index, list) {
				if (angular.isFunction(iterator)) {
					results.push(iterator.call(context, value, index, list));
				}
			});
			return results;
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
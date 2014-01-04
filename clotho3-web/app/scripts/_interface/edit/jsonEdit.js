'use strict';

angular.module('clotho.interface')
	.directive('jsonEdit', function () {
		return {
			restrict: 'A',
			require: 'ngModel',
			scope: {
				model : '=ngModel'
			},
			link: function(scope, element, attrs, ngModelCtrl) {

				function fromUser(text) {

					console.log(text, typeof text);
					try {
						return angular.fromJson(text);
					} catch (err) {
						ngModelCtrl.$setValidity(false);
						return text;
					}
				}

				function toUser(object) {
					// better than JSON.stringify(), because it formats + filters $$hashKey etc.
					return angular.toJson(object, true);
				}

				function isValidJson (model) {
					var flag = true;
					try {
						angular.fromJson(model);
					} catch (err) {
						flag = false;
					}
					return flag;
				}

				// push() if faster than unshift(), and avail. in IE8 and earlier (unshift isn't)
				ngModelCtrl.$parsers.push(fromUser);
				ngModelCtrl.$formatters.push(toUser);

				scope.$watch('model', function(newValue, oldValue) {
					if (newValue != oldValue && isValidJson(newValue)) {
						ngModelCtrl.$setViewValue(toUser(newValue));
						ngModelCtrl.$render();
					}
				}, true); // MUST use objectEquality (true) here, for some reason..
			}
		};
	});

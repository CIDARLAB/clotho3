'use strict';

/*
 example usage: <textarea json-edit="myObject" rows="8" class="form-control"></textarea>

 jsonEditing is a string which we edit in a textarea. we try parsing to JSON with each change. when it is valid, propagate model changes via ngModelCtrl use isolate scope to prevent model propagation when invalid - will update manually. cannot replace with template, or will override ngModelCtrl, and not hide behind facade will override element type to textarea and add own attribute ngModel tied to jsonEditing

 As far as I know, there is currently no way to achieve this using $parsers (other than one of the function errors and kills the pipeline)

 //todo - refactor so can use ngModel... but can't create element inside.
 */

angular.module('clotho.editor')
	.directive('jsonEdit', function () {
		return {
			restrict: 'A',
			require: 'ngModel',
			template: '<textarea ng-model="jsonEditing"></textarea>',
			replace : true,
			scope: {
				model: '=jsonEdit'
			},
			link: function (scope, element, attrs, ngModelCtrl) {

				function string2JSON(text) {
					try {
						return angular.fromJson(text);
					} catch (err) {
						setInvalid();
						return text;
					}
				}

				function JSON2String(object) {
					// better than JSON.stringify(), because it formats + filters $$hashKey etc.
					// NOTE that this will remove all $-prefixed values
					return angular.toJson(object, true);
				}

				function setEditing (value) {
					scope.jsonEditing = JSON2String(value);
				}

				function updateModel (value) {
					scope.model = string2JSON(value);
				}

				function setValid() {
					ngModelCtrl.$setValidity('json', true);
				}

				function setInvalid () {
					ngModelCtrl.$setValidity('json', false);
				}

				function isValidJson(model) {
					var flag = true;
					try {
						angular.fromJson(model);
					} catch (err) {
						flag = false;
					}
					return flag;
				}

				//init
				setEditing(scope.model);

				//check for changes going out
				scope.$watch('jsonEditing', function (newval, oldval) {
					if (newval != oldval) {
						if (isValidJson(newval)) {
							setValid();
							updateModel(newval);
						} else {
							setInvalid();
						}
					}
				}, true);

				//check for changes coming in
				scope.$watch('model', function (newval, oldval) {
					if (newval != oldval) {
						setEditing(newval);
					}
				}, true);

			}
		};
	});
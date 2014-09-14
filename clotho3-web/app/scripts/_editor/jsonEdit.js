'use strict';

/*
 example usage: <textarea json-edit="myObject" rows="8" class="form-control"></textarea>

 jsonEditing is a string which we edit in a textarea. we try parsing to JSON with each change. when it is valid, propagate model changes via ngModelCtrl use isolate scope to prevent model propagation when invalid - will update manually. cannot replace with template, or will override ngModelCtrl, and not hide behind facade will override element type to textarea and add own attribute ngModel tied to jsonEditing

 As far as I know, there is currently no way to achieve this using $parsers (other than one of the function errors and kills the pipeline)

 //todo - refactor so can use ngModel... but can't create element inside.
 //todo - NEED TO REMOVE REPLACE - see GH#412
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
					// NOTE that angular.toJson will remove all $$-prefixed values
					// alternatively, use JSON.stringify(object, null, 2);
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
					var flag = angular.isDefined(model);
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
          if (isValidJson(newval)) {
            setValid();
            updateModel(newval);
          } else {
            setInvalid();
          }
				}, true);

				//check for changes coming in
				scope.$watch('model', function (newval, oldval) {
          setEditing(newval);
				}, true);

			}
		};
	})

/*

When used on a textarea (or input), allows the direct editing of a JSON object (any valid JSON - string, array, boolean, etc.). Invalid JSON will set the value to undefined. It is recommended you do not allow saving etc. while this field is invalid.

Handles validation under the `json` attribute.

example usage:

 <form name="myForm">
  <textarea json-editor ng-model="myObject" rows="8" name="myFormElement" class="form-control"></textarea>
  <p ng-show="myForm.myFormElement.$error.json">JSON is invalid!</p>
 </form>

todo - NEED TO ADD SUPPORT FOR PRIMITIVES (doesn't work for strings because try to serialize)

 */

.directive('jsonEditor', function () {
	return {
		restrict: 'A',
		require: 'ngModel',
		link: function (scope, element, attrs, ngModelCtrl) {

			function isValidJson(model) {
        var flag = angular.isDefined(model);
				try {
					angular.fromJson(model);
				} catch (err) {
					flag = false;
				}
				return flag;
			}

			function string2JSON(text) {
				try {
					return angular.fromJson(text);
				} catch (err) {
					//returning undefined results in a parser error as of angular-1.3-rc.0, and will not go through $validators
					//return undefined
					return text;
				}
			}

			function JSON2String(object) {
				// NOTE that angular.toJson will remove all $$-prefixed values
				// alternatively, use JSON.stringify(object, null, 2);
				return angular.toJson(object, true);
			}

      /*
      todo - only want to update view (formatters) if objects differ
      if (!angular.equals(object, angular.fromJson(ngModelCtrl.$viewValue))) {

      }
      */

			//$validators is an object, where key is the error
			ngModelCtrl.$validators.json = isValidJson;

			//array pipelines
			ngModelCtrl.$parsers.push(string2JSON);
			ngModelCtrl.$formatters.push(JSON2String);
		}
	}
});


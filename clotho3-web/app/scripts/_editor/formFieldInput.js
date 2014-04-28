'use strict';

//note - not good enough to work. doesn't handle transclusion well (e.g. form-field directive). works ok standalone.

angular.module('clotho.editor')
/**
 * @name formFieldInput
 *
 * @description
 * Given a type and model, generates an input field to handle data of that type and binds it to the model. Element will be replaced, so you can define normal input attributes.
 *
 */
	.directive('formFieldInput', function (ClothoSchemas, $compile, $parse) {

		return {
			restrict: 'A',
			link: function (scope, element, attrs) {

				function replaceElement(passedType) {
					var javascriptType = ClothoSchemas.formTypeMap[passedType] || false;

					var inputElement;
					if (javascriptType) {
						//if we have an input field, just create an input and extend with attrs given
						if (javascriptType['input']) {
							inputElement = angular.element('<input>');
							inputElement.attr(javascriptType['input']);
							inputElement.attr({
								'ng-model': attrs.formFieldModel
							});
						}
						//otherwise, need to handle a bit differently (probably an object)
						else {
							inputElement = angular.element('<textarea>');
							inputElement.attr({
								'json-edit': attrs.formFieldModel,
								rows: 3
							});
						}
					}
					else {
						//didn't map, handle as default, allow specification in UI
						inputElement = angular.element('<input>');
						inputElement.attr({
							type: 'text',
							'ng-model': attrs.formFieldModel
						});
					}


					//simplified version of template merge, for use in compile (template will automatically handle) to reapply old attributes
					angular.forEach(attrs, function (val, key) {
						if (key == 'class') {
							inputElement.addClass(val); //won't handle SVG (read-only)
						}
						else if (key.charAt(0) != '$' && !inputElement.attr(key)) {
							inputElement.attr(attrs.$attr[key], val);
						}
					});

					element.attr('type', null);
					element.attr('ng-model', null);
					element.attr('json-edit', null);

					//todo - ngModel collisions for json-edit because in element.data()
					element.replaceWith($compile(inputElement)(scope));
				}

				//todo - won't render on first run
				var watcher;
				watcher = scope.$watch(function () {
					return $parse(attrs.formFieldType)(scope)
				}, function (newval, oldval) {
					if (newval != oldval) {
						replaceElement(newval);
						watcher();
					}
				});
			}
		}
	});

/*
		return {
			restrict: 'A',
			replace: true,
			*/
/*
			 //don't want to create isolate scope
			 scope : {
			 model : '@formFieldModel',
			 type : '@formFieldType',
			 },
			  note - we would use template because angular handles the template merge well. Create an angular element because its cleaner then just use the outerHTML

			template: function (tElement, tAttrs) {

				var javascriptType = ClothoSchemas.formTypeMap[tAttrs.formFieldType] || false;

				var childElement;
				if (javascriptType) {
					//if we have an input field, just create an input and extend with attrs given
					if (javascriptType['input']) {
						childElement = angular.element('<input>');
						childElement.attr(javascriptType['input']);
						childElement.attr({
							'ng-model': tAttrs.formFieldModel
						});
					}
					//otherwise, need to handle a bit differently (probably an object)
					else {
						childElement = angular.element('<textarea>');
						childElement.attr({
							'json-edit': tAttrs.formFieldModel,
							rows: 3
						});
					}
				}
				else {
					//didn't map, handle as default, allow specification in UI
					childElement = angular.element('<input>');
					childElement.attr({
						type: 'text',
						'ng-model': tAttrs.formFieldModel
					});
				}

				console.log(childElement[0].outerHTML);

				return childElement[0].outerHTML;
			},
			link : function (scope, element, attrs) {

				scope.$watch(function () {
					return attrs.formFieldType;
				}, function (newval, oldval) {
					console.log(newval, oldval);
					if (newval != oldval) {

							element.attr('form-field-type', newval);
							element.attr('type', null);
							console.log(element[0].outerHTML);
							element.replaceWith($compile(element[0].outerHTML)(scope));
					}
				});
			}
			 compile: function (tElement, tAttrs) {

			 var javascriptType = ClothoSchemas.formTypeMap[tAttrs.formFieldType] || false;

			 var childElement;
			 if (javascriptType) {
			 //if we have an input field, just create an input and extend with attrs given
			 if (javascriptType['input']) {
			 childElement = angular.element('<input>');
			 childElement.attr(javascriptType['input']);
			 childElement.attr({
			 'ng-model' : tAttrs.formFieldModel
			 });
			 }
			 //otherwise, need to handle a bit differently (probably an object)
			 else {
			 childElement = angular.element('<textarea>');
			 childElement.attr({
			 'json-edit': tAttrs.formFieldModel,
			 rows: 3
			 });
			 }
			 }
			 else {
			 //didn't map, handle as default, allow specification in UI
			 childElement = angular.element('<input>');
			 childElement.attr({
			 type : 'text',
			 'ng-model' : tAttrs.formFieldModel
			 });
			 }


			 //simplified version of template merge, for use in compile (template will automatically handle) to reapply old attributes
			 angular.forEach(tAttrs, function (val, key) {
			 if (key == 'class') {
			 childElement.addClass(val); //won't handle SVG (read-only)
			 }
			 else if (key.charAt(0) != '$' && !childElement.attr(key)) {
			 childElement.attr(tAttrs.$attr[key], val);
			 }
			 });

			 tElement.replaceWith(childElement);


			 return function (scope, element, attrs) {

			 */
/*
			 function replaceElement () {
			 var javascriptType = ClothoSchemas.formTypeMap[attrs.formFieldType] || false;

			 var childElement;
			 if (javascriptType) {
			 //if we have an input field, just create an input and extend with attrs given
			 if (javascriptType['input']) {
			 childElement = angular.element('<input>');
			 childElement.attr(javascriptType['input']);
			 childElement.attr({
			 'ng-model' : attrs.formFieldModel
			 });
			 }
			 //otherwise, need to handle a bit differently (probably an object)
			 else {
			 childElement = angular.element('<textarea>');
			 childElement.attr({
			 'json-edit': attrs.formFieldModel,
			 rows: 3
			 });
			 }
			 }
			 else {
			 //didn't map, handle as default, allow specification in UI
			 childElement = angular.element('<input>');
			 childElement.attr({
			 type : 'text',
			 'ng-model' : attrs.formFieldModel
			 });
			 }


			 //simplified version of template merge, for use in compile (template will automatically handle) to reapply old attributes
			 angular.forEach(attrs, function (val, key) {
			 if (key == 'class') {
			 childElement.addClass(val); //won't handle SVG (read-only)
			 }
			 else if (key.charAt(0) != '$' && !childElement.attr(key)) {
			 childElement.attr(attrs.$attr[key], val);
			 }
			 });

			 element.replaceWith($compile(childElement)(scope));

			 }

			 scope.$watch(function () {
			 return attrs.formFieldType;
			 }, function (newval, oldval) {
			 console.log(newval, oldval);
			 if (newval) {
			 $timeout(function () {
			 //replaceElement();
			 });
			 }
			 });
			 }


		}
	});
*/

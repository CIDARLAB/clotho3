'use strict';

angular.module('clotho.editor')
/**
 * @name formFieldEnumeration
 *
 * @description
 * you can pass fields explicitly, and/or a sharable with a schema. Only one is required (probably sharable).
 *
 * Each type is given it's own section (fields, schema, sharable), where sharable fields not matching the schema will appear in the sharable section
 *
 * - fields: pass fields to generate form fields for a specified list of fields. Fields are given priority over schema and sharable
 * - sharable: pass a sharable (containing a schema) to generate a form containing its fields, and the fields of its schema
 */
	.directive('formFieldEnumeration', function (Clotho, ClothoSchemas, $q, $compile) {

		//todo - handle ID differently
		//future - incorporate select and radios

		function generateDynamicFields (fields) {
			var allFields = angular.element('<div ng-form>');

			angular.forEach(fields, function(field) {

				var formField = angular.element('<form-field>');
				formField.attr({
					name : field.name,
					horizontal : 'formHorizontal'
				});

				/*
				//this doesn't compile properly in transcluded scope
				var formFieldInput = angular.element('<div>');
				formFieldInput.attr({
					'form-field-input' : '',
					'form-field-type' : field.type,
					'form-field-model' : 'sharable.' + field.name,
					'name' : field.name,
					'ng-required' : !!field.required,
					'placeholder' : field.description
				});
				*/

				var javascriptType = ClothoSchemas.formTypeMap[field.type] || false;

				var inputElement;
				if (javascriptType) {
					//if we have an input field, just create an input and extend with attrs given
					if (javascriptType['input']) {
						inputElement = angular.element('<input>');
						inputElement.attr(javascriptType['input']);
						inputElement.attr({
							'ng-model': 'sharable.' + field.name
						});
					}
					//otherwise, e.g. object
					else {
						inputElement = angular.element('<textarea>');
						inputElement.attr({
							'json-edit': 'sharable.' + field.name,
							rows: 3
						});
					}
				}
				else {
					//didn't map, handle as default, allow specification via JSON
					inputElement = angular.element('<textarea>');
					inputElement.attr({
						'json-edit': 'sharable.' + field.name,
						rows: 1,
						placeholder: "Edit JSON directly, use quotes for strings"
					});
				}

				inputElement.attr({
					'placeholder' : field.description
				});

				if (field.required) {
					inputElement.attr({
						'ng-required' : true
					});
				}

				formField.append(inputElement);
				allFields.append(formField);
			});

			return allFields;
		}

		return {
			restrict: 'A',
			scope: {
				fields: '=?',
				sharable: '=?',
				stripBasicFields : '@?'
			},
			controller: function($scope, $element, $attrs) {},
			link: function (scope, element, attrs) {

				if (angular.isUndefined(attrs.fields) &&
						angular.isUndefined(attrs.sharable)) {
					element.html('no schema information passed');
					return;
				}

				//styling pass-through
				scope.formHorizontal = scope.$parent.formHorizontal;

				//container elements for each set of fields
				var schemaFieldsElement = angular.element('<div>'),
					sharableFieldsElement = angular.element('<div>'),
					definedFieldsElement = angular.element('<div>');

				//regenerate if fields change
				scope.$watch('fields', function (newval) {
					if (!!newval) {
						definedFieldsElement.replaceWith($compile(generateDynamicFields(newval))(scope));
					}
				}, true);

				//update sharable keys change - schema and sharable specific fields
				scope.$watch(function () {
					return _.keys(scope.sharable)
				}, function(newkeys) {
					if (scope.sharable && scope.sharable.schema) {

						Clotho.get(scope.sharable.schema)
						.then(function (schema) {

							ClothoSchemas.downloadSchemaDependencies(schema)
							.then(function(compiledSchema) {

								var schemaFieldNames = _.pluck(compiledSchema.fields, 'name');
								var sharableFieldNames = _.keys(scope.sharable);
								var sharableOnlyFieldNames = _.difference(sharableFieldNames, schemaFieldNames);

								//todo - handle field types if defined --- where to define?
								var sharableOnlyFields = _.map(sharableOnlyFieldNames, function (fieldName) {
									return {
										name : fieldName
									}
								});

								if (scope.stripBasicFields) {
									_.remove(compiledSchema.fields, function (field) {
										return angular.isDefined(ClothoSchemas.sharableBasicFields[field.name]);
									});
									_.remove(sharableOnlyFields, function (field) {
										return angular.isDefined(ClothoSchemas.sharableBasicFields[field.name]);
									});
								}

								schemaFieldsElement = $compile(generateDynamicFields(compiledSchema.fields))(scope);

								sharableFieldsElement = $compile(generateDynamicFields(sharableOnlyFields))(scope);

								replaceFieldsView();
							});
						});
					}
					//no schema...
					else {
						schemaFieldsElement = angular.element('<div class="alert alert-warning">Sharable has no schema...</div>');

						var strippedFields = newkeys;
						_.remove(strippedFields, function (field) {
							return angular.isDefined(ClothoSchemas.sharableBasicFields[field]);
						});

						var sharableFieldNames = _.map(strippedFields, function (fieldName) {
							return {
								name : fieldName
							}
						});
						sharableFieldsElement = $compile(generateDynamicFields(sharableFieldNames))(scope);
						replaceFieldsView();
					}
				}, true);

				function replaceFieldsView () {
					element.empty();
					//add elements
					element.prepend(sharableFieldsElement);
					element.prepend(schemaFieldsElement);
					element.prepend(definedFieldsElement);
				}


			}
		}
	});
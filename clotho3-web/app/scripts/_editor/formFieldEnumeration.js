angular.module('clotho.editor')
/**
 * @name formFieldEnumeration
 *
 * you can pass fields > schema > sharable (in terms of priority given to fields generated). Only one is requried.
 * - fields: pass fields to generate form fields for a specified list of fields. Fields are given priority over schema and sharable
 * - schema: pass a schema to generate a list of fields from the schema
 * - sharable: pass a sharable (containing a schema) to generate a form containing its fields, and the fields of its schema
 */
	.directive('formFieldEnumeration', function (Clotho, ClothoUtils, $q, $compile) {

		function generateDynamicForm (fields) {
			var fulltext = "";

			//todo - need BSON types from server (GH #99)
			//todo - options for radio and select?
			/*
			 string -> test vs. textarea
			 id not editable
			 */
			var typeMap = {
				'string' : 'text',
				'boolean' : 'checkbox',
				'select' : 'select',
				'radio' : 'radio',
				'object' : 'object'
			};

			angular.forEach(fields, function(field) {

				//fixme -- don't all convert to text if not in typemap
				var type = typeMap[field.type] || 'text';
				if (type == '?') field.type == 'text';
				var required = field.required ? "required" : "";

				var htmlText_pre = '<form-field name="' + field.name + '" removable-field="' + field.name + '">';
				var htmlText_post = '</form-field>';
				var inputText;

				switch (type) {
					case 'object' : {
						inputText = '<textarea json-edit rows="3" ' + required + ' ng-model="sharable.'+field.name+'"></textarea>';
						break;
					}
					case 'radio' : {
						angular.forEach(field.options, function(value, key) {
							inputText += '<input type="radio" value="'+value+'" ng-model="sharable.'+field.name+'>' + value;
						});
						break;
					}
					case 'select': {
						var optionsText = "";
						//todo - use ng-options (what are options) + attach array to scope
						angular.forEach(field.options, function(value, key) {
							optionsText = optionsText + '<option value="'+value+'">'+ value + '</option>';
						});

						inputText = '<select ' + required + ' ng-model="sharable.'+field.name+'">' + optionsText + '</select>';
						break;
					}
					default: {
						inputText = '<input type="' + type + '" ' + required + 'ng-model="sharable.'+field.name+'" >';
						break;
					}

				}

				fulltext += htmlText_pre + inputText + htmlText_post;
			});

			return fulltext;
		}

		return {
			restrict: 'A',
			scope: {
				fields: '=?',
				schema: '=?',
				sharable: '=?'
			},
			controller: function($scope, $element, $attrs) {

				$scope.generateFields = function (opts) {

					var fields = opts.fields || $scope.fields || [];
					var schema = opts.schema || $scope.schema || {fields: []};
					var sharable = opts.sharable || $scope.sharable || {schema : {fields: []}};

					console.log(fields, schema, sharable);

					var concatenatedFields = fields.concat(schema.fields);

					//add sharable fields
					_.forEach(sharable, function (val, key) {
						concatenatedFields.push({
							name : key
						});
					});

					console.log(_.pluck(concatenatedFields, 'name'));

					var uniqueFields = _.uniq(concatenatedFields, function (field) {
						return field.name;
					});

					$element.html(generateDynamicForm(uniqueFields));
					$compile($element.contents())($scope);

				};

			},
			link: function (scope, element, attrs) {

				if (angular.isUndefined(scope.fields) && angular.isUndefined(scope.schema) && angular.isUndefined(scope.sharable)) {
					element.html('no information passed');
					return;
				}

				function updateWithSchema (schema) {
					return ClothoUtils.downloadSchemaDependencies(schema)
					.then(function(compiledSchema) {
						scope.generateFields({schema : compiledSchema});
					});
				}

				//regenerate if fields change
				scope.$watch('fields', function (newval) {
					if (!!newval) {
						console.log('updating fields');

						scope.generateFields({fields: newval});
					}
				}, true);

				//update schema with superClass fields
				scope.$watch('schema', function (newval) {
					if (!!newval && !!newval.fields && newval.fields.length) {
						console.log('updating schema');

						updateWithSchema(newval);
					}
				}, true);

				//update sharable schema based on superClass
				scope.$watch(function () {
					return _.keys(scope.sharable)
				}, function(newkeys) {
					if (!!scope.sharable) {
						console.log('updating sharable');

						//todo - don't get schema and just pass fields

						Clotho.get(scope.sharable.schema)
						.then(function(retrievedSchema) {
							updateWithSchema(retrievedSchema);
						})
					}
				}, true);

			}
		}
	});
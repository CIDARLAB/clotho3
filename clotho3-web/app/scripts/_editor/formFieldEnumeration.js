angular.module('clotho.editor')
/**
 * @name formFieldEnumeration
 *
 * you can pass fields > schema > sharable (in terms of priority given to fields generated). Only one is requried.
 * - fields: pass fields to generate form fields for a specified list of fields. Fields are given priority over schema and sharable
 * - schema: pass a schema to generate a list of fields from the schema
 * - sharable: pass a sharable (containing a schema) to generate a form containing its fields, and the fields of its schema
 */
	.directive('formFieldEnumeration', function (Clotho, $q, $compile) {

		function generateDynamicForm (fields) {
			var fulltext = "";

			//todo - need BSON types from server (GH #99)
			/*
			 string -> test vs. textarea
			 when use select?
			 what to do for objects?
			 id not editable
			 */
			var typeMap = {
				'string' : 'text',
				'boolean' : 'checkbox'
			};

			angular.forEach(fields, function(field) {

				var type = field.type || 'text';
				if (type == '?') field.type == 'text';
				var required = field.required ? "required" : "";

				var htmlText_pre = '<form-field name="' + field.name + '">';
				var htmlText_post = '</form-field>';
				var inputText;

				switch (type) {
					case "textarea": {
						inputText = '<textarea rows="2" ' + required + ' ng-model="sharable.'+field.name+'"></textarea>';
						break;
					}
					case "select": {
						var optionsText = "";
						//todo(low) - use ng-options + attach array to scope
						angular.forEach(field.options, function(value, key) {
							optionsText = optionsText + '<option value="'+value+'">'+ value + '</option>';
						});

						inputText = '<select ' + required + ' ng-model="sharable.'+field.name+'">' + optionsText + '</select>';
						break;
					}
					case "sharable": {
					}
					//todo - add filedrop support, and radio. checkbox works.
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

				//todo - clean up
				//will update schema as reference (whether sharable.schema or just a schema object)
				$scope.getSuperClasses = function (schema) {

					//initial check
					if (!schema.superClass) {
						return schema
					}

					var finalSchema = angular.copy(schema);
					var promiseChain = $q.when();
					var reachedBottom = $q.defer();

					function getSuperClass (passedSchema) {
						if (passedSchema.superClass) {
							promiseChain.then(function () {
								//testing console.log('retriving ' + passedSchema.superClass);
								return Clotho.get(passedSchema.superClass)
								.then(function (retrieved) {
									//testing console.log('retrieved ' + retrieved.id + ' - ' + retrieved.name, _.pluck(retrieved.fields, 'name'), retrieved);
									finalSchema.fields = finalSchema.fields.concat(retrieved.fields);
									//testing console.log('finalSchema now', _.pluck(finalSchema.fields, 'name'));
									return getSuperClass(retrieved)
								});
							});
						} else {
							reachedBottom.resolve();
						}
					}

					getSuperClass(schema);

					return reachedBottom.promise.then(function() {
						return promiseChain
					})
					.then(function (chain) {
						return finalSchema;
					});
				};

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

				//regenerate if fields change
				scope.$watch('fields', function (newval) {
					if (!!newval) {
						console.log('updating fields');
						scope.generateFields();
					}
				}, true);

				//update schema with superClass fields
				scope.$watch('schema', function (newval) {
					if (!!newval && !!newval.fields && newval.fields.length) {
						console.log('updating schema');
						scope.getSuperClasses(newval)
							.then(function() {
								scope.generateFields();
							});
					}
				}, true);

				//update sharable schema based on superClass
				scope.$watch('sharable', function(newval) {
					if (!!newval) {
						console.log('updating sharable');
						Clotho.get(newval.schema)
						.then(function(retrievedSchema) {
							return scope.getSuperClasses(retrievedSchema)
						})
						.then(function(compiledSchema) {
							console.log(compiledSchema);
							scope.generateFields({schema : compiledSchema});
						});
					}
				}, true);

			}
		}
	});
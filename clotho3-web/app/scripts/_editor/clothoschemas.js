'use strict';

angular.module('clotho.editor')
  .service('ClothoSchemas', function EditorSchemas(Clotho, $q) {

		var SCHEMA_BUILTIN = 'org.clothocad.core.schema.BuiltInSchema';
		var SCHEMA_CLOTHOSCHEMA = 'org.clothocad.core.schema.ClothoSchema';

		var retrievedSchemas = $q.defer();

		var sharableTypes = {
			"Instance" : {
				readable : "Instance",
				editor_template_url : 'views/_editor/sharable.html',
				scaffold : {
					schema: false,
					language: "JSONSCHEMA"
				}
			},
			"Function": {
				readable : "Function",
				editor_template_url : 'views/_editor/function.html',
				scaffold : {
					schema: "org.clothocad.core.datums.Function",
					language: "JSONSCHEMA"
				}
			},
			"Schema": {
				readable : "Schema",
				editor_template_url : 'views/_editor/schema.html',
				scaffold : {
					schema: SCHEMA_CLOTHOSCHEMA,
					language: "JSONSCHEMA"
				}
			},
			"View" :  {
				readable : "View",
				editor_template_url : 'views/_editor/view.html',
				scaffold : {
					schema: "org.clothocad.core.datums.View",
					language: "JSONSCHEMA"
				}
			}
		};

		var accessTypes = [
			{name:'Public', value:'PUBLIC'},
			{name:'Private', value:'PRIVATE'},
			{name:'Read Only', value:'READONLY'}
		];

		var constraintTypes = [
			{name:'RegExp', value:'regex'},
			{name: 'Not Null', value: 'notnull'}
		];

		var primitiveToJava = {
			"string" : "java.lang.String",
			"number" : "java.lang.Long",
			"boolean" : "java.lang.Boolean",
			"object" : "java.util.HashMap",
			"array" : "java.util.List"
		};

		/* QUERIES */

		//hold our query promises
		var schemaQueries = {};

		//query for built-in schemas, removing instances in sharableTypes
		schemaQueries.builtIn = Clotho.query({"schema": SCHEMA_BUILTIN})
		.then(function (resultSchemas) {
			_.remove(resultSchemas, function (schema) {
				return !!sharableTypes[schema.name];
			});
			return resultSchemas;
		});

		//query for clothoSchemas
		schemaQueries.clothoSchemas = Clotho.query({"schema": SCHEMA_CLOTHOSCHEMA});

		//get them all, resolve promise
		$q.all(schemaQueries).then(function (results) {
			var allSchemas = [];
			angular.forEach(results, function (schemas, type) {
				allSchemas.concat(schemas);
			});
			retrievedSchemas.resolve(allSchemas);
		});

		/* FUNCTIONALITY */

		//returns schema of a sharable, or null
		function determineSchema (sharable) {
			return sharable.schema || null;
		}

		//determine whether a sharable is a schema
		function isSchema (sharable) {
			var sharableSchema = determineSchema(sharable);
			return sharableSchema == SCHEMA_BUILTIN || sharableSchema == SCHEMA_CLOTHOSCHEMA;
		}

		function isFunction (sharable) {
			return determineSchema(sharable) == sharableTypes.Function.scaffold.schema;
		}

		function isView (sharable) {
			return determineSchema(sharable) == sharableTypes.View.scaffold.schema;
		}

		//determines main type: Instance, Function, View, Schema
		function determineType (sharable) {
			if (isSchema(sharable)) {
				return "Schema";
			} else if (isFunction(sharable)) {
				return 'Function';
			} else if (isView(sharable)) {
				return 'View';
			} else {
				return 'Instance';
			}
		}

		//creates scaffold given SchemaName, does not create on server
		function createScaffold (schemaName) {
			var scaffold;
			if (sharableTypes[schemaName]) {
				scaffold = sharableTypes[schemaName].scaffold;
			} else {
				scaffold = angular.extend({schema : schemaName}, sharableTypes.Instance.scaffold)
			}
			return scaffold;
		}

		return {
			retrievedSchemas : retrievedSchemas.promise,
			sharableTypes : sharableTypes,
			accessTypes : accessTypes,
			constraintTypes : constraintTypes,
			primitiveToJava : primitiveToJava,

			isSchema : isSchema,
			determineType : determineType,
			determineSchema : determineSchema,
			createScaffold : createScaffold
		}
  });

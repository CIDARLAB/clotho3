'use strict';

angular.module('clotho.foundation')
  .service('ClothoSchemas', function EditorSchemas(Clotho, $q) {

		var SCHEMA_ALL = 'org.clothocad.core.schema.Schema';
		var SCHEMA_BUILTIN = 'org.clothocad.core.schema.BuiltInSchema';
		var SCHEMA_CLOTHOSCHEMA = 'org.clothocad.core.schema.ClothoSchema';

		var retrievedSchemas = $q.defer();

		//Main types of instances, and basic scaffold for creating each
		var sharableTypes = {
			"Instance" : {
				readable : "Instance",
				editor_template_url : 'views/_editor/sharable.html',
				schema: false,
				scaffold : {
					language: "JSONSCHEMA"
				},
				color : '#4488cc'
			},
			"Function": {
				readable : "Function",
				editor_template_url : 'views/_editor/function.html',
				schema: "org.clothocad.core.datums.Function",
				scaffold : {
					language: "JSONSCHEMA"
				},
				color : '#66bb66'
			},
			"Schema": {
				readable : "Schema",
				editor_template_url : 'views/_editor/schema.html',
				schema: SCHEMA_CLOTHOSCHEMA,
				scaffold : {
					language: "JSONSCHEMA"
				},
				color : '#eeaa55'
			},
			"View" :  {
				readable : "View",
				editor_template_url : 'views/_editor/view.html',
				schema: "org.clothocad.core.datums.View",
				scaffold : {
					language: "JSONSCHEMA"
				},
				color : '#dd33dd'
			}
		};

		//Basic fields to show for an instance for basic distinguishing
		var sharableBasicFields = {
			"name" : {
				name : 'name',
				type : 'string',
				description : 'Name of the Instance, given by author'
			},
			"id" : {
				name : 'id',
				type : 'string',
				description : 'Unique ID referring to this object'
			},
			"schema" : {
				name : 'schema',
				type : 'string',
				description : 'Pattern describing contents and organization of instance data'
			},
			"description" : {
				name : 'description',
				type : 'string',
				description : 'Description of object, written by author'
			},
			"author" : {
				name: 'author',
				type: 'string',
				description : 'User who created this object'
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
			"boolean" : "java.lang.Boolean",
			"number" : "java.lang.Long",
			"object" : "java.util.HashMap",
			"array" : "java.util.List",
			"string" : "java.lang.String"
		};

		//determine a JSON field / javascript variable's type
		function determineFieldType (value) {
			if (value === null) {
				return 'null';
			} else if (angular.isBoolean(value)) {
				return 'boolean';
			} else if (angular.isNumber(value)) {
				return 'number';
			} else if (angular.isArray(value)) {
				return 'array';
			} else if (angular.isObject(value)) {
				return 'object';
			} else {
				return 'string';
			}
		}

		function colorByType (type) {
			if (sharableTypes[type]) {
				return sharableTypes[type].color;
			} else {
				return '#31b0d5';
			}
		}

		/* QUERIES */

		Clotho.query({"schema" : SCHEMA_ALL})
		.then(function (resultSchemas) {
			_.remove(resultSchemas, function (schema) {
				return !!sharableTypes[schema.name];
			});
				retrievedSchemas.resolve(resultSchemas)
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

		//determine if schemas is a built in schema
		function isBuiltIn (sharable) {
			return determineSchema(sharable) == SCHEMA_BUILTIN;
		}

		//determine if a schema is a Clotho written schema
		function isClothoSchema (sharable) {
			return determineSchema(sharable) == SCHEMA_CLOTHOSCHEMA;
		}

		function isFunction (sharable) {
			return determineSchema(sharable) == sharableTypes.Function.schema;
		}

		function isView (sharable) {
			return determineSchema(sharable) == sharableTypes.View.schema;
		}

		//determines main type: Instance, Function, View, Schema
		function determineInstanceType (sharable) {
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
				scaffold = angular.extend({schema : sharableTypes[schemaName].schema }, sharableTypes[schemaName].scaffold);
			} else {
				scaffold = angular.extend({schema : schemaName}, sharableTypes.Instance.scaffold)
			}
			return scaffold;
		}

		/**
		 * @description Returns a new object containing just the basic fields of a sharable, defined by sharableBasicFields
		 * @param {Object} sharable Sharable to prune
		 * @returns {Object} New object with limited fields
		 */
		function pruneToBasicFields (sharable) {
			return angular.map(sharableBasicFields, function (obj, key) {
				return sharable[key]
			});
		}

		return {
			retrievedSchemas : retrievedSchemas.promise,
			sharableTypes : sharableTypes,
			accessTypes : accessTypes,
			constraintTypes : constraintTypes,
			primitiveToJava : primitiveToJava,

			isSchema : isSchema,
			isBuiltIn : isBuiltIn,
			isClothoSchema : isClothoSchema,

			determineInstanceType : determineInstanceType,
			determineFieldType : determineFieldType,
			determineSchema : determineSchema,
			createScaffold : createScaffold,
			colorByType : colorByType,

			sharableBasicFields : sharableBasicFields,
			pruneToBasicFields : pruneToBasicFields
		}
  });

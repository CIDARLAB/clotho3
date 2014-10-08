'use strict';

angular.module('clotho.foundation')
  .service('ClothoSchemas', function EditorSchemas(Clotho, Debug, $filter, $q) {

		var Debugger = new Debug('clothoSchemas', '#992299');

		var SCHEMA_ALL = 'org.clothocad.core.schema.Schema';
		var SCHEMA_BUILTIN = 'org.clothocad.core.schema.BuiltInSchema';
		var SCHEMA_CLOTHOSCHEMA = 'org.clothocad.core.schema.ClothoSchema';

		var retrievedSchemas = $q.defer();
		//map from schema ID to schema name, ready once retrievedSchemas resolves
		var schemaNameMap = {};

		//Main types of instances, and basic scaffold for creating each
		var sharableTypes = {
			"Instance" : {
				readable : "Instance",
				editor_template_url : 'views/_editor/sharable.html',
				schema: false,
				scaffold : {
					language: "JSONSCHEMA"
				},
				class : 'info' // blue
			},
			"Function": {
				readable : "Function",
				editor_template_url : 'views/_editor/function.html',
				schema: "org.clothocad.core.datums.Function",
				scaffold : {},
				class : 'success' //green
			},
			"Schema": {
				readable : "Schema",
				editor_template_url : 'views/_editor/schema.html',
				schema: SCHEMA_CLOTHOSCHEMA,
				scaffold : {
					language: "JSONSCHEMA"
				},
				class : 'danger' //orange
			},
			"View" :  {
				readable : "View",
				editor_template_url : 'views/_editor/view.html',
				schema: "org.clothocad.core.datums.View",
				scaffold : {},
				class : 'warning' //pink
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


		function escapeRegexp (reg) {
			return (reg.toString()).replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&');
		}
		//matches spaces unless between single or double quotes
		var spaceRegexp = /[^\s"']+|"([^"]*)"|'([^']*)'/;
		var spaceRegexpEscaped =  spaceRegexp.toString();

		var formTypeMap = {
			"boolean" : {
				type : "boolean",
				input : {
					type : 'checkbox'
				}
			},
			"number" : {
				type : "number",
				input : {
					type : 'number'
				}
			},
			"string" : {
				type : "string",
				input : {
					type : 'text'
				}
			},
			"date" :{
				type : "date",
				input : {
					type : 'date'
				}
			},
			"array" : {
				type : "array"
			},
			"org.bson.types.ObjectId" : {
				type : "id",
				input : {
					type : 'text'
				}
			},
			"object" : {
				type : "object"
			},
			//suggest do not just use this catch, but allow further specification in UI or something
			"default" : {
				type : "string",
				input : {
					type: "text"
				}
			}
		};

		// has object 'input' if can be displayed as input. add all as attrs.
		// otherwise will have to handle another way (e.g. json-edit directive)
		var javaToJavascript = {
			"java.lang.Boolean" : {
				type : "boolean",
				input : {
					type : 'checkbox'
				}
			},
			"java.lang.Number" : {
				type : "number",
				input : {
					type : 'number'
				}
			},
			"java.lang.Long" : {
				type : "number",
				input : {
					type : 'number'
				}
			},
			"java.lang.Integer" : {
				type : "number",
				input : {
					type : 'number',
					"ng-pattern" : "/^[-+]?\\d+$/"
				}
			},
			"java.lang.Short" : {
				type : "number",
				input : {
					type : 'number',
					min : -32767,
					max : 32767
				}
			},
			"java.lang.String" : {
				type : "string",
				input : {
					type : 'text'
				}
			},
			"java.awt.Color" :{
				type : "color",
				input : {
					type : 'color'
				}
			},
			"java.util.Date" :{
				type : "date",
				input : {
					type : 'date'
				}
			},
			"java.util.Set" : {
				type : "array"
			},
			"java.util.List" : {
				type : "array"
			},
			"org.bson.types.ObjectId" : {
				type : "id",
				input : {
					type : 'text'
				}
			},
			"java.util.Map" : {
				type : "object"
			},
			"java.util.HashMap" : {
				type: "object"
			},
			//suggest do not just use this catch, but allow further specification in UI or something (json-edit directive)
			"default" : {
				type : "string",
				input : {
					type: "text"
				}
			}
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

		function typeToColorClass (type) {
			if (sharableTypes[type]) {
				return sharableTypes[type].class;
			} else {
				return 'default';
			}
		}

		/* QUERIES */

		Clotho.query({"schema" : SCHEMA_ALL}, {mute : true})
		.then(function (resultSchemas) {
			//populate name map
			angular.forEach(resultSchemas, function (schema) {
				schemaNameMap[schema.id] = schema.name;
			});

			retrievedSchemas.resolve(resultSchemas);
		});

		//depends on resolution of retrievedSchemas
		function mapSchemaIdToName (id) {
			return schemaNameMap[id] || '';
		}

		function retrieveSuperclass (schemaSharable) {

			var superclassId = schemaSharable.superClass;

			if (!superclassId) {
				return $q.when(null);
			}

			return retrievedSchemas.promise.then(function (schemas) {
				var inPromise = $filter('filter')(schemas, {id : superclassId})[0];
				if (angular.isEmpty(inPromise)) {
					return $q.when(inPromise);
				} else {
					return Clotho.get(superclassId);
				}
			});
		}

		var getParentSchemaIds = function (schemaId) {
			return Clotho.run('clotho.functions.schema.getParents', [schemaId], {mute : true});
		};

		var getSuperclassFields = function (schemaSharable) {

			var fields = [];

			//initial checks
			if (angular.isUndefined(schemaSharable) || !schemaSharable.superClass) {
				return $q.when(fields);
			}

			var promiseChain = $q.when();
			var reachedBottom = $q.defer();

			//todo - retrieveSuperclass returns null -- loop until don't return null instead of crazy promises (or use getParentSchemaIds)
			function getSuperClass (passedSchema) {
				if (passedSchema.superClass) {
					promiseChain.then(function () {
						return retrieveSuperclass(passedSchema).then(function (retrieved) {
							fields = fields.concat(retrieved.fields);
							return getSuperClass(retrieved);
						}, function (err) {
							Debugger.warn('couldnt get parent schema ' + passedSchema.superClass);
							reachedBottom.reject(fields);
						});
					});
				} else {
					reachedBottom.resolve();
				}
			}

			getSuperClass(schemaSharable);

			return reachedBottom.promise.then(function() {
				return promiseChain;
			})
			.then(function (chain) {
				return fields;
			});

		};

		//fixme - DEPRECATE
		var downloadSchemaDependencies = function (schemaSharable) {

			//initial checks
			if (angular.isUndefined(schemaSharable)) {
				return $q.when();
			}
			if (!schemaSharable.superClass) {
				return $q.when(schemaSharable);
			}

			var finalSchema = angular.copy(schemaSharable);
			var promiseChain = $q.when();
			var reachedBottom = $q.defer();

			function getSuperClass (passedSchema) {
				if (passedSchema.superClass) {
					promiseChain.then(function () {
						return Clotho.get(passedSchema.superClass)
							.then(function (retrieved) {
								//testing console.log('retrieved ' + retrieved.id + ' - ' + retrieved.name, _.pluck(retrieved.fields, 'name'), retrieved);
								finalSchema.fields = finalSchema.fields.concat(retrieved.fields);
								//testing console.log('finalSchema now', _.pluck(finalSchema.fields, 'name'));
								return getSuperClass(retrieved)
							}, function (err) {
								Debugger.warn('couldnt get parent schema ' + passedSchema.superClass);
								//couldn't get parent schema
								reachedBottom.reject(finalSchema);
							});
					});
				} else {
					reachedBottom.resolve();
				}
			}

			getSuperClass(schemaSharable);

			return reachedBottom.promise.then(function() {
				return promiseChain;
			})
			.then(function (chain) {
				return finalSchema;
			});
		};

		/* FUNCTIONALITY */

		//forgiving check to see if potential sharable
		function isSharable (sharable) {
			return angular.isObject(sharable) && !angular.isEmpty(sharable) && angular.isDefined(sharable.id) && angular.isDefined(sharable.schema);
		}

		//returns schema of a sharable, or null
		function determineSchema (sharable) {
			return ( !angular.isEmpty(sharable) && angular.isDefined(sharable.schema) ) ? sharable.schema : null;
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
			var sharableSchema = determineSchema(sharable);
			return sharableSchema == sharableTypes.Function.schema || sharableSchema == 'org.clothocad.core.datums.Module';
		}

		function isView (sharable) {
			return determineSchema(sharable) == sharableTypes.View.schema;
		}

		//determines main type: Instance, Function, View, Schema -- ASYNC so kinda slow...
		function determineSharableType (sharable) {
			return Clotho.run('clotho.functions.schema.determineSharableType', [sharable.id], {mute : true});
		}

		//this is a synchronous version for things like autocomplete where we're running a lot of these
    //returns empty string if does not appear to be a sharable
		function dirtyDetermineType (sharable) {
			if (isSchema(sharable)) {
				return 'Schema';
			} else if (isFunction(sharable)) {
				return 'Function';
			} else if (isView(sharable)) {
				return 'View';
			} else if (isSharable(sharable)) {
				return 'Instance';
			} else {
        return '';
      }
		}

		//todo - currently this only checks name, but should check for child as well
		//hack
		function isInstanceOfSchema (sharable, schema) {

			if (angular.isUndefined(schema)) {
				return false;
			}
			//hack - handle NucSeqs
			if (schema == 'org.clothocad.model.NucSeq') {
				return /clotho.demo.sequence.*/.test(sharable.id) || _.indexOf(['org.clothocad.model.NucSeq', 'Part', 'Vector'], determineSchema(sharable)) > 0;
			}

			return determineSchema(sharable) == schema;
		}

		var sharableIconMap = {
			Instance : "glyphicon glyphicon-file",
			Function : "glyphicon glyphicon-play-circle",
			View : "glyphicon glyphicon-picture",
			Schema : "glyphicon glyphicon-cog",
			default : "glyphicon glyphicon-file"
		};

		//returns class for icon of sharable type given
		function determineSharableIcon (type) {
			return sharableIconMap[type] || sharableIconMap["default"];
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
			downloadSchemaDependencies : downloadSchemaDependencies,
			getSuperclassFields : getSuperclassFields,
			getParentSchemaIds : getParentSchemaIds,

			sharableTypes : sharableTypes,
			//checks if a key is a basic sharable field (name, author, id, etc.)
			isBasicField : function (fieldType) {
				return angular.isDefined(sharableBasicFields[angular.lowercase(fieldType)]);
			},
			isPrimitiveField : function (fieldType) {
				return angular.isDefined(primitiveToJava[angular.lowercase(fieldType)]);
			},
			accessTypes : accessTypes,
			constraintTypes : constraintTypes,
			primitiveToJava : primitiveToJava,
			formTypeMap : formTypeMap,
			javaToJavascript : javaToJavascript,

			isSharable : isSharable,
			isFunction : isFunction,
			isView : isView,
			isSchema : isSchema,
			isInstance : function (sharable) {
				return dirtyDetermineType(sharable) == 'Instance';
			},
			isBuiltIn : isBuiltIn,
			isClothoSchema : isClothoSchema,
			mapSchemaIdToName : mapSchemaIdToName,
			isInstanceOfSchema : isInstanceOfSchema,

			dirtyDetermineType : dirtyDetermineType,
			determineSharableType : determineSharableType,
			determineSharableIcon : determineSharableIcon,
			determineFieldType : determineFieldType,
			determineSchema : determineSchema,
			createScaffold : createScaffold,
			typeToColorClass : typeToColorClass,

			sharableBasicFields : sharableBasicFields,
			pruneToBasicFields : pruneToBasicFields
		}
  });

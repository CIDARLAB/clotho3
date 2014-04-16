angular.module('clotho.utils', ['clotho.core'])
	.service('ClothoUtils', function($q, $http, Clotho) {

		/**
		 * @name verifyUUID
		 * @param uuid
		 * @description Verify that a UUID string is valid
		 */
		var validUUID = function (uuid) {
			return angular.isString(uuid) && uuid.length == 16 && (/[a-zA-Z0-9]{16}/).test(uuid);
		};

		/***********
		VIEWS
		************/

		//generates absolute url for file in view
		function generateWidgetUrl (viewId, url) {
			return 'widgets/' + viewId + (!!url ? '/' + url : '');
		}

		//TEMPORARY client polyfill for Clotho.get() to retrieve a view for testing
		var clientGetView = function(viewId) {
			return $http.get(generateWidgetUrl(viewId) + '/model.json')
				.then(function(data){
					return data.data
				});
		};

		/*
		 Download view dependencies recursively.

		 Will go into imported views sequentially, download *their* dependencies, then bubble up to current view.

		 Returns promise that is fulfilled when all dependencies downloaded
		 */
		var downloadViewDependencies = function (view) {
			//create array of promises of nested dependencies from importedView
			var nestedDeps = [];
			_.forEach(view.importedViews, function (id, alias) {
				//return Clotho.get(id).then(function(retrievedView) {    //when server
				nestedDeps.push(clientGetView(id)               //testing
					.then(function (retrievedView) {
						return downloadViewDependencies(retrievedView);
					})
				);
			});

			//download nested dependencies, from deepest import bubbling up
			return $q.all(nestedDeps)
				//after imported dependencies downloaded, mixin current view's dependencies
				.then(function() {
					var relativeDeps = [];
					_.forEach(view.dependencies, function (dep) {
						relativeDeps.push(generateWidgetUrl(view.id, dep));
					});

					return $clotho.extensions.mixin(relativeDeps);
				})
				.then(function() {
					return view;
				});
		};

		/***********
		 SCHEMAS
		 ************/

		//todo - clean up
		var downloadSchemaDependencies = function (schema) {

			//initial checks
			if (angular.isUndefined(schema)) {
				return $q.when();
			}
			if (!schema.superClass) {
				return $q.when(schema);
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
				return promiseChain;
			})
			.then(function (chain) {
				return finalSchema;
			});
		};

		return {
			validUUID : validUUID,

			downloadViewDependencies : downloadViewDependencies,
			generateWidgetUrl : generateWidgetUrl,

			downloadSchemaDependencies : downloadSchemaDependencies
		}
	});

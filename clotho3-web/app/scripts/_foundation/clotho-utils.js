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

    var getViewInfo = clientGetView;

		/*
		 Download view dependencies recursively.

		 Will go into imported views sequentially, download *their* dependencies, then bubble up to current view.

		 Returns promise that is fulfilled when all dependencies downloaded
		 */
		var downloadViewDependencies = function (view) {
			//create array of promises of nested dependencies from importedView
			var nestedDeps = [];
			angular.forEach(view.importedViews, function (id, alias) {
				nestedDeps.push(getViewInfo(id)
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
					angular.forEach(view.dependencies, function (dep) {
						relativeDeps.push(generateWidgetUrl(view.id, dep));
					});
					return $clotho.extensions.mixin(relativeDeps);
				})
				.then(function() {
					return view;
				});
		};

		return {
			validUUID : validUUID,

      getViewInfo : getViewInfo,
			downloadViewDependencies : downloadViewDependencies,
			generateWidgetUrl : generateWidgetUrl
		}
	});

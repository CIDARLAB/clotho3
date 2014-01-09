angular.module('clotho.clothoDirectives', [])
/**
 * @name clotho-show
 *
 * @usage <div clotho-show="VIEW_ID"></div>
 */
.directive('clothoShow', function ($q, $http, $timeout, $browser, $rootScope, Clotho, PubSub) {

	function generateWidgetUrl (viewId, url) {
		return 'widgets/' + viewId + (!!url ? '/' + url : '');
	}

	//simple polyfill for Clotho.get() to retrieve a view
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
	var downloadDependencies = function (view) {

		//create array of promises of nested dependencies from importedView
		var nestedDeps = [];
		_.forEach(view.importedViews, function (id, alias) {
			//return Clotho.get(id).then(function(retrievedView) {    //when server
			nestedDeps.push($q.when(clientGetView(id))                //testing
				.then(function (retrievedView) {
					return downloadDependencies(retrievedView);
				})
			);
		});

		//download nested dependencies, from deepest import bubbling up
		return $q.all(nestedDeps)
			//after imported dependencies downloaded, mixin current view's dependencies
			.then(function() {
				console.log('downloaded nested dependencies of view ' + view.id);
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




	return {
		terminal: true,
		restrict: 'E',
		scope: {
			id: '@',
			callback: '=?'
		},
		controller: function ($scope, $element, $attrs) {

		},
		link: function linkFunction (scope, element, attrs) {

			if (!scope.id) return;

			//todo - basic check to make sure proper form

			//config element
			element.addClass('clothoWidget');

			//retrieve view
			$q.when(clientGetView('123456789'))           //testing
				//Clotho.get(scope.id)                      //when server handles
				.then(function(view){
					return downloadDependencies(view);
				})
				.then(function (view) {

					console.log(view);

					//configure dictionary
					view.dictionary = angular.extend({}, view.dictionary, view.importedViews);
					view.dictionary.id = view.id;



					if (view.bootstrap) {
						//hack-y creating custom module so we can set some stuff up without taking the module creation out of the user's control
						var customModuleName = view.id + '-additions';
						var dependencies = view.includeExtensionsModule === true ? ['clotho.extensions'] : [];
						angular.module(customModuleName, dependencies)
							.run(function($rootScope) {
								//extend scope with dictionary
								angular.extend($rootScope, view.dictionary);

								/**
								 * @name prefixUrl
								 *
								 * @description Function which will prefix partial URLs appropriately, e.g. in ng-include
								 * @param url {string} URL of partial, relative to View root
								 * @param specifyView {string} ID of view. Pass nothing to default to this view's id
								 */
									//todo - handle specifying view by name
								$rootScope.prefixUrl = function (url, specifyView) {
									return generateWidgetUrl(specifyView ? specifyView : view.id, url)
								}
							});

						//Modules : overwrite some default services
						var modules = [];
						modules.push(function($provide) {
							$provide.value('$anchorScroll', angular.noop);
							$provide.value('$browser', $browser);
						});
						//Modules : add declared dependencies. May also list in module definition
						modules = modules.concat(view.bootstrap.modules, customModuleName);

						//assumes module includes a template index.html. Otherwise nothing will be shown.
						element.html('<div ng-include="prefixUrl(\'index.html\')"></div>');

						//if pass controller, attach it to DOM for compilation
						if (!!view.controller) {
							attrs.$set('ng-controller', view.controller)
						}

						//BOOTSTRAP
						element.data('$injector', null);
						angular.bootstrap(element, modules);
					}
					else {
						//handle according to type


					}




					//CALLBACK
					//can also use run clause of module, but not passed the element
					$timeout(function() {
						angular.isFunction(scope.callback) && scope.callback(element);
						PubSub.trigger('clothoShow:' + scope.id, [scope.id, element, view])
					});
				})
		}
	};
})

	//todo
	.directive('clothoModal', function () {

	});

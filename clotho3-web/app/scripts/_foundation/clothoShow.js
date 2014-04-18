angular.module('clotho.clothoDirectives')
/**
 * @name clotho-show
 *
 * @usage <div clotho-show="VIEW_ID"></div>
 */
.directive('clothoShow', function ($q, $http, $timeout, $browser, $rootScope, $compile, Clotho, PubSub, ClothoUtils) {

	var generateWidgetUrl = ClothoUtils.generateWidgetUrl;

	//client polyfill for Clotho.get() to retrieve a view for testing
	var clientGetView = function(viewId) {
		return $http.get(generateWidgetUrl(viewId) + '/model.json')
		.then(function(data){
			return data.data
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

			console.log('directive linked');

			scope.$watch('id', function (newval, oldval) {
				if (!!newval) {

					console.log(newval);

					//todo - basic check to make sure proper form

					//config element
					element.addClass('clothoWidget');

					//retrieve view
					$q.when(clientGetView(scope.id))              //testing
						//Clotho.get(scope.id)                      //when server handles
						.then(function(view){
							return ClothoUtils.downloadViewDependencies(view);
						})
						.then(function (view) {

							//configure dictionary
							view.dictionary = angular.extend({}, view.dictionary, view.importedViews);
							view.dictionary.id = view.id;

							if (view.bootstrap) {
								//creating custom module so we can set some stuff up without taking the module creation out of the user's control
								var customModuleName = view.id + '-additions';
								var dependencies = view.bootstrap.excludeExtensionsModule !== false ? ['clotho.extensions'] : [];
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
										$rootScope.prefixUrl = function (url, specifyView) {
											return generateWidgetUrl(specifyView ? specifyView : view.id, url)
										}
									});

								//Modules : overwrite some default services
								var modules = [];
								modules.push(function($provide, $compileProvider) {
									$provide.value('$anchorScroll', angular.noop);
									$provide.value('$browser', $browser);

									//service to handle inline scripts
									$provide.service('lazyScripts', ['$q', '$timeout', '$document', function ($q, $timeout, $document) {

										var promises = [];

										this.register = function (url) {
											promises.push($clotho.extensions.mixin(url));
										};

										$timeout(function() {
											$q.all(promises).then(function() {
												//broadcast event
												$document.triggerHandler('WidgetContentLoaded');
											})
										});
									}]);

									//override script directive to handle relative assets, register with lazyScripts
									/*
									 Angular provides was to dynamically load templates with dynamic names via ng-include. We needed to download scripts, relative to the path of the .html partial calling them. (i.e. we had a directory with a file to include, and wanted the .html file to declare the scripts etc. it needed itself).

									 Inline JS and CSS will load fine, but not scripts with src like the one below.

									 This directive will patch angular so that you can use ng-src and pass in an interpolated src, or one generated by a function. Here, we are assuming that prefixUrl() is a method on the $rootScope (or accessible some other way). For example, prefixUrl('x.js') may generate 'namespaced/path/x.js'

									 You can extend the functionality of the message broadcasted above to only run after all tags are downloaded... Just register all the scripts you need to download and emit the message after $q.all()
									 */
									$compileProvider.directive('script', ['$parse', '$rootScope', 'lazyScripts', function($parse, $rootScope, lazyScripts) {
										return {
											restrict: 'E',
											terminal: true,
											compile: function(element, attr) {
												if (attr.ngSrc) {
													var scriptUrl = $parse(attr.ngSrc)($rootScope);
													lazyScripts.register(scriptUrl);
												}
											}
										};
									}]);
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

							//if don't pass bootstrap clause, handle appropriately
							else {

								//this is what will be inserted
								var htmlString;

								//if pass index.html, include template
								if (_.indexOf(view.files, 'index.html') >= 0) {

									htmlString = '<div ';

									if (!!view.controller) {
										htmlString += 'ng-controller="' + view.controller + '" ';
									}

									htmlString += 'ng-include="prefixUrl(\'index.html\')"></div>';

								}
								//if only one file, handle according to type
								else if (view.files.length == 1) {

									//todo - other types?

									var reg_image = /([a-z\-_0-9\/\:\.]*\.(jpg|jpeg|png|gif))/i;

									//if image, create img tag
									if (reg_image.test(view.files[0])) {
										htmlString = '<img src="' +
											generateWidgetUrl(view.id, view.files[0]) + '"' +
											'alt="view '+view.id+'" />';
									}
								}

								angular.extend(scope, view.dictionary);
								element.html($compile(htmlString)(scope));

							}




							//CALLBACK
							//can also use run clause of module, but not passed the element
							$timeout(function() {
								angular.isFunction(scope.callback) && scope.callback(element);
								PubSub.trigger('clothoShow:' + scope.id, [scope.id, element, view])
							});
						});
				}
			});
		}
	};
});
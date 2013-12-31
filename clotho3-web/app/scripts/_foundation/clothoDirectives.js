'use strict';

/*
 This is meant for read-only modifications to the model - Can't use ngModel $formatters with promises (angular-1.1.5). Calling $setViewValue affects model and propagates, which is often undesired (e.g. if revcomp a sequence for display, don't want to change model).

 note : form of ngModel and functions of higher arity --
 it is assumed that if ngModel is an array, it is in the correct format as if being passed to a function.apply().... if it is not an array, it is assumed the function is of single arity, and ngModel is wrapped in an array as the only value (i.e. [ngModel.$modelValue] )

 note : updating parent scope
 There are two ways to do this.
 (1) Because an isolate scope is created, you could pass in $parent.<model> which should do normal angular binding.
 (2) If you want to update the model only with the run function, add the attribute tag clotho-run-update-model="true" ... this will update the model with the result of the run function once it is complete.

 @example
 <p clotho-run="lowercase" ng-model="'HEY THERE'"></p> will output <p>hey there</p>
 */
angular.module('clotho.clothoDirectives', [])
.directive('clothoRun', function(Clotho) {

	var inputsVal = {input: true, textarea : true, select: true};

	return {
		restrict : 'A',
		require : 'ngModel',
		scope : true,
		link: function (scope, element, attrs, ngModel) {

			//config
			var useVal = !!inputsVal[angular.lowercase(element[0].nodeName)];
			if (useVal) {
				//avoid flicker
				ngModel.$render = angular.noop;
			}
			var updateParent = false;

			//command, args
			scope.$watch(function() {
				return attrs.clothoRun
			}, function(newval, oldval) {
				if (!!newval) runFunction(ngModel.$modelValue);
			});

			//model changes
			scope.$watch(function() {
				return ngModel.$modelValue
			}, function(newval, oldval) {
				//console.log(newval);
				runFunction(newval);
			});

			//update model?
			scope.$watch(function() {
				return attrs.clothoRunUpdateModel
			}, function(newval, oldval) {
				updateParent = !!newval;
			});

			//form array out of arguments if not an array
			var parseInput = function(input) {
				return angular.isArray(input) ? input : [input];
			};

			var updateParentModel = function(newModel) {
				if (updateParent) {
					//todo - make sure passes up to $parent
					ngModel.$setViewValue(newModel);
				}
			};

			var updateElement = function(newval) {
				var method = useVal ? 'val' : 'text';
				element[method](newval);
			};

			var runFunction = function(input) {
				input = parseInput(input);

				return Clotho.run(attrs.clothoRun, input).then(function(result) {
					console.log(result);
					updateParentModel(result);
					updateElement(result);
				});
			};


		}
	}
})

/**
 * @name clotho-show
 *
 * @usage <div clotho-show="VIEW_ID"></div>
 */
	.directive('clothoShow', function ($q, $timeout, $browser, $rootScope, Clotho, PubSub) {

		function generateWidgetUrl (url, viewId) {
			return 'widgets/' + viewId + '/' + url;
		}

		/* VIEW OBJECT - what is expected from server	*/
		var exampleReturnedView = {
			//id of the view, used
			id: "123456789",

			//views declared as dependencies so can use filenames locally
			importedViews : {
				"otherView" : "987654321"
			},

			//ATM no need to pass this to the client... these can be requested lazily
			files: [
				'index.html', //it is expected a view to be displayed will have a partial named index.html. This is the template that will be used when added to the DOM
				'lazyPartial.html'
			],

			//files to download, URLs passed by server (likely to be namespaced by id)
			dependencies: [
				'widgets/123456789/external-module.js',
				'widgets/123456789/widgetModule.js',
				'widgets/123456789/widgetController.js'
			],

			//to extend the scope. client will add in id and imported views
			dictionary: {
				"dictString" : "My String",
				"dictObject" : {
					"myKey" : "myValue"
				}
			},

			//if true, will include module clotho.extensions so can add components easily
			includeExtensionsModule : true,

			//attach a controller for the whole widget - use angular.module('clotho.extensions') if not using own module and use boolean to include module
			controller: '123456789_parentController',


			//information for bootstrap
			bootstrap: {
				modules: ['123456789']
			}

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

				//config element
				element.addClass('clothoWidget');

				//retrieve view
				$q.when(exampleReturnedView)                //testing
					//Clotho.get(scope.id)                      //when server handles
					.then(function (view) {

						//todo - get importedView dependencies

						$clotho.extensions.mixin(view.dependencies).then(function() {

							//configure dictionary
							angular.extend(view.dictionary, view.importedViews);
							view.dictionary.id = view.id;



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
										return generateWidgetUrl(url, specifyView ? specifyView : view.id)
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


							//CALLBACK
							//can also use run clause of module, but not passed the element
							$timeout(function() {
								angular.isFunction(scope.callback) && scope.callback(element);
								PubSub.trigger('clothoShow:' + scope.id, [scope.id, element, view])
							});
						})
					});
			}
		};
	});

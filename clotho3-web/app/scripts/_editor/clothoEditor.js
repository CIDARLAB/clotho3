'use strict';

/**
 * @name clotho-editor
 *
 * @param sharable {string|object} If String, ID of sharable. Otherwise the sharable directly.
 * @param editMode {*} pass false, null, undefined to start in view mode. Pass anything that passes !!var to start in edit mode
 *
 *
 */
angular.module('clotho.editor').directive('clothoEditor', function (Clotho, $compile, $parse, $http, $templateCache, $filter, $q, Debug, ClothoSchemas) {

	var Debugger = new Debug('Editor', '#dd44dd');

	return {
		restrict: 'A',
		require: ['ngModel', '^form'],
		scope: {
			inputSharable: '=ngModel',
			editMode: '=?'
		},
		controller: function ($scope, $element, $attrs) {
			/******
			 GENERAL
			 ******/

				//add novalidate to the form to bypass browser validation
			$attrs.$set('novalidate', "novalidate");

			//if no name is set, set so can use for form validation
			if (!angular.isString($attrs.name)) {
				$attrs.$set('name', 'sharableEditor')
			}

			//for testing, easily log scope
			$scope.logScope = function () {
				Debugger.log($scope);
			};

			//scoped variables for including templates
			$scope._objBaseFields = 'views/_editor/_baseFields.html';
			$scope._formActions = 'views/_editor/_formActions.html';

			/** model **/

			var sharableId;

			/*********
			 COMPILATION
			 **********/

			// should check for custom template and use if exists. otherwise, do the form generation for the generic
			//note - assumes scope already has sharable bound for compilation
			$scope.getPartialAndCompile = function (type) {
				$scope.showJsonEditor = false;

				//todo - handle instance-specific templates

				if (angular.isEmpty(type)) {
					$element.html('<p class="text-center">Please select an Object</p>');
					return;
				}

				$http.get(ClothoSchemas.sharableTypes[type].editor_template_url, {cache: $templateCache})
				.success(function (data) {
					var el = $compile(data)($scope);
					$element.html(el);
				})
				.error(function (data) {
					Debugger.error('Could not retrieve template for type ' + type);
					$element.html('<p class="text-center">Please select an Object</p>');
				});
			};

			//process input sharable to see if object, or id string
			$scope.processInputSharable = function (sharable) {

				if (angular.isEmpty(sharable)) {
					Debugger.log('sharable is undefined / empty');
					$scope.getPartialAndCompile();
					return;
				}

				$scope.sharable = {};

				Debugger.debug('processing sharable: ', sharable);

				if (angular.isObject(sharable) && !angular.isEmpty(sharable)) {
					Debugger.log('sharable object passed');
					sharableId = sharable.id | '';
					$scope.sharable = sharable;
				}
				else if (angular.isString(sharable)) {
					//if its a string, call clotho.get()
					Debugger.log('passed a string, assuming a valid ID: ' + sharable);
					sharableId = sharable;
				} else {
					Debugger.warn('sharable must be an object or a string, you passed: ', sharable);
					sharableId = '';
					$scope.sharable = {};
				}

				//if its a string, it's an ID so get it, otherwise just return a wrapper promise for the object
				var getObj = (angular.isEmpty($scope.sharable) && sharableId != '') ? Clotho.get(sharableId) : $q.when($scope.sharable);

				getObj.then(function (result) {
					$scope.sharable = result;
					$scope.getPartialAndCompile(ClothoSchemas.determineInstanceType(result));
				});

			};

			//init with input sharable
			$scope.processInputSharable($scope.inputSharable);

			//controller API
			return {
				removeField: function (name) {
					Debugger.warn('(DEPRECATED) - removing fields from sharables is not supported');
					Debugger.log('removing key ' + name);
					delete $scope.sharable[name];
				}
			}
		},
		compile: function compile(tElement, tAttrs, transclude) {

			return {
				pre: function preLink(scope, iElement, iAttrs, controllers) {

				},
				post: function postLink(scope, iElement, iAttrs, controllers) {

					scope.isValidJson = function (model) {
						var flag = true;
						try {
							angular.fromJson(model);
						} catch (err) {
							flag = false;
						}
						return flag;
					};

					/* config */

					scope.formHorizontal = scope.$eval(iAttrs.formHorizontal);

					scope.formCtrl = controllers[1];

					scope.edit = function () {
						scope.editMode = true;
					};

					scope.toggleJsonEdit = function () {
						scope.showJsonEditor = !scope.showJsonEditor;
					};

					scope.reset = function () {
						Clotho.get(scope.sharable.id).then(function (result) {
							scope.sharable = result;
							scope.formCtrl.$setPristine();
						});
					};

					scope.save = function () {
						if (!scope.isValidJson(scope.sharable)) {
							//this shouldn't happen often... JSON editor should not propagate to model until valid
							Clotho.alert('Your JSON is invalid... please fix it');
						} else {
							Clotho.set(scope.sharable).then(function (id) {
								if (!scope.sharable.id) {
									console.log('adding ID', id);
									scope.sharable.id = id;
								}
							});
							scope.editMode = false;
						}
					};

					scope.discard = function () {
						scope.reset();
						scope.editMode = false;
					};

					scope.destroy = function () {
						Clotho.destroy(scope.sharable.id).then(function () {
							scope.editMode = false;
							scope.sharable = null;
							scope.id = null;
							scope.processInputSharable({});
						});
					};

					/* watchers */

					//listen for collector_reset and reget & recompile
					Clotho.listen('collector_reset', function Editor_onCollectorReset() {
						scope.processInputSharable();
					}, scope);


					//watch for internal PubSub Changes
					var sharableWatcher = angular.noop;
					scope.$watch('sharable', function (newval, oldval) {
						if (!!newval && newval.id && (!oldval || newval.id != oldval.id)) {
							sharableWatcher();
							sharableWatcher = Clotho.watch(newval.id, function (newObj) {
								console.log('\n\n\n\nclothoEditor CLOTHO WATCH', newObj);
								scope.sharable = newObj;
							}, scope);
						}
					});

					scope.$watch('inputSharable', function (newval, oldval) {
						if (!oldval || (!!newval && !!oldval && newval.id != oldval.id)) {
							console.log(newval, oldval);
							scope.editMode = false;
							scope.processInputSharable(newval);
						}
					});

					scope.$watch('editMode', function (newval, oldval) {
						Debugger.log('edit mode: ', newval);
					});
				}
			}
		}
	}
});
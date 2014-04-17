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

			/*********
			 COMPILATION
			 **********/

			// should check for custom template and use if exists. otherwise, do the form generation for the generic
			$scope.getPartialAndCompile = function (type, obj) {
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

				$scope.sharable = {};

				if (angular.isEmpty(sharable)) {
					Debugger.log('sharable is undefined / empty');
					$scope.getPartialAndCompile();
					return;
				}

				Debugger.debug('processing sharable: ', sharable);

				if (angular.isObject(sharable) && !angular.isEmpty(sharable)) {
					Debugger.log('sharable object passed');
					$scope.id = sharable.id | '';
					$scope.sharable = sharable;
				}
				else if (angular.isString(sharable)) {
					//if its a string, call clotho.get()
					Debugger.log('passed a string, assuming a valid ID: ' + sharable);
					$scope.id = sharable;
				} else {
					Debugger.warn('sharable must be an object or a string, you passed: ', sharable);
					//todo - better exit if nothing present
					$scope.id = '';
					$scope.sharable = {};
				}


				//if its a string, it's an ID so get it, otherwise just return a wrapper promise for the object
				var getObj = (angular.isEmpty($scope.sharable) && $scope.id != '') ? Clotho.get($scope.id) : $q.when($scope.sharable);

				getObj.then(function (result) {
					$scope.sharable = result;
					$scope.getPartialAndCompile(ClothoSchemas.determineType(result), result);
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

					/* config */

					scope.formCtrl = controllers[1];

					scope.edit = function () {
						scope.editMode = true;
					};

					scope.toggleJsonEdit = function () {
						scope.showJsonEditor = !scope.showJsonEditor;
					};

					scope.reset = function () {
						scope.formCtrl.$setPristine();
						Clotho.get(scope.id).then(function (result) {
							scope.sharable = result;
						});
					};

					scope.save = function () {
						Clotho.set(scope.sharable);
						scope.editMode = false;
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
					Clotho.watch(scope.id, scope.sharable, scope);

					scope.$watch('inputSharable', function (newval, oldval) {
						if (!oldval || !!newval && !!oldval && newval.id != oldval.id) {
							scope.editMode = false;
						}
						scope.processInputSharable(newval);
					});

					scope.$watch('editMode', function (newval, oldval) {
						Debugger.log('edit mode: ', newval);
					});

					scope.$watch('sharable', function (newval, oldval) {
						Debugger.log('sharable updated ', newval);
					});

				}
			}
		}
	}
});
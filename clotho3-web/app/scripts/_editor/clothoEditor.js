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
		require: '^form',
		scope: {
			inputSharable: '=sharable',
			editMode: '=?'
		},
		controller: function ($scope, $element, $attrs) {
			/******
			 GENERAL
			 ******/

			//scoped variables for including templates, set up immediately
			$scope._objBaseFields = 'views/_editor/_baseFields.html';
			$scope._formActions = 'views/_editor/_formActions.html';

			var passedId,
				passedModel,
				retrievedModel;

			/*********
			 COMPILATION
			 **********/

			//process input sharable to see if object, or id string
			$scope.processInputSharable = function (sharable) {

				if (angular.isEmpty(sharable)) {
					Debugger.warn('sharable is undefined / empty');
					$scope.getPartialAndCompile();
					return;
				}

				Debugger.log('processing input');


				//reset what existed before
				$scope.sharable = {}; //can be reset because re-assign in compilation
				retrievedModel = null;
				passedModel = null;
				passedId = '';

				Debugger.log('processing sharable: ', sharable);

				if (angular.isObject(sharable) && !angular.isEmpty(sharable)) {
					Debugger.log('sharable object passed');
					passedModel = sharable;
					passedId = passedModel.id | '';
				} else if (angular.isString(sharable)) {
					//if its a string, call clotho.get()
					Debugger.log('passed a string, assuming a valid ID: ' + sharable);
					passedId = sharable;
				} else {
					Debugger.warn('sharable must be an object or a string, you passed: ', sharable);
					$scope.getPartialAndCompile();
					return;
				}

				//if its a string, it's an ID so get it, otherwise just return a wrapper promise for the object
				var getObj = (angular.isEmpty(passedModel) && passedId != '') ?
					Clotho.get(passedId) :
					$q.when(passedModel);

				getObj.then(function (result) {
					console.log('got object', result);
					retrievedModel = result;
					$scope.getPartialAndCompile(ClothoSchemas.determineInstanceType(result));
				});
			};

			//gets partial for type of sharable and compiles editor HTML
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

						$scope.sharable = retrievedModel;
						$scope.panelClass = ClothoSchemas.sharableTypes[type].class || 'default';

						console.log('compiling', type);
						console.log(passedModel == retrievedModel, retrievedModel == $scope.sharable, passedModel, retrievedModel, $scope.sharable);

						var el = $compile(data)($scope);
						$element.html(el);
					})
					.error(function (data) {
						Debugger.error('Could not retrieve template for type ' + type);
						$element.html('<p class="text-center">Please select an Object</p>');
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
		link: function postLink(scope, element, attrs, formCtrl) {

			/* config */

			//add novalidate to the form to bypass browser validation
			attrs.$set('novalidate', 'novalidate');

			//if no name is set, set so can use for form validation
			if (!angular.isString(attrs.name)) {
				attrs.$set('name', 'sharableEditor');
			}

			scope.formHorizontal = scope.$eval(attrs.formHorizontal);

			/* functionality */

			scope.isValidJson = function (model) {
				var flag = true;
				try {
					angular.fromJson(model);
				} catch (err) {
					flag = false;
				}
				return flag;
			};

			scope.edit = function () {
				scope.editMode = true;
			};

			scope.toggleJsonEdit = function () {
				scope.showJsonEditor = !scope.showJsonEditor;
			};

			scope.reset = function () {
				Clotho.get(scope.sharable.id).then(function (result) {
					scope.sharable = result;
					formCtrl.$setPristine();
				});
			};

			scope.save = function () {
				if (!scope.isValidJson(scope.sharable)) {
					//this shouldn't happen often... JSON editor should not propagate to model until valid
					Clotho.alert('Your JSON is invalid... please fix it');
				} else {
					Clotho.set(scope.sharable).then(function (id) {
						if (!scope.sharable.id) {
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


			//watch for internal PubSub Changes for given sharable
			var sharableWatcher = angular.noop;
			scope.$watch('sharable', function (newval, oldval) {
				if (!angular.isEmpty(newval) && ( !oldval || (newval.id && newval.id != oldval.id))) {
					Debugger.log('sharable is different', newval, oldval);
					//cancel the last watcher
					sharableWatcher();
					sharableWatcher = Clotho.watch(newval.id, function (newObj) {
						scope.sharable = newObj;
					}, scope);
				}
			});

			scope.$watch('inputSharable', function (newval, oldval) {
				Debugger.log('inputSharable: ', newval, oldval);
				//no old value, or new and old but id's don't match, or no new id field
				if (!oldval || (!!newval && !!oldval && (!newval.id || newval.id != oldval.id))) {
					scope.processInputSharable(newval);
				}
			});

			scope.$watch('editMode', function (newval, oldval) {
				Debugger.log('edit mode: ', newval);
			});
		}
	}
});
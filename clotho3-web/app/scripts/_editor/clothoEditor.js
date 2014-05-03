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

	var baseFieldTemplate = 'views/_editor/_baseFields.html',
		formActionsTemplate = 'views/_editor/_formActions.html';

	return {
		restrict: 'A',
		require: '^form',
		scope: {
			sharable: '=',
			editMode: '=?',
			formHorizontal : '=?'
		},
		link: function postLink(scope, element, attrs, formCtrl) {

			/* config */

			//add novalidate to the form to bypass browser validation
			attrs.$set('novalidate', 'novalidate');

			//if no name is set, set so can use for form validation
			if (!angular.isString(attrs.name)) {
				attrs.$set('name', 'sharableEditor');
			}

			scope.$watch('formHorizontal', function (newval) {
				element.toggleClass('form-horizontal', newval);
			});

			//scoped variables for including templates
			scope._objBaseFields = baseFieldTemplate;
			scope._formActions = formActionsTemplate;

			//DEPRECATED
			scope.removeField = function (name) {
				Debugger.warn('(DEPRECATED) - removing fields from sharables is not supported');
				Debugger.log('removing key ' + name);
				delete scope.sharable[name];
			};

			/* compilation */

			var retrievedModel;

			//process input sharable to see if object, or id string
			scope.processInputSharable = function (sharable) {

				if (angular.isEmpty(sharable)) {
					Debugger.warn('sharable is undefined / empty');
					scope.getPartialAndCompile();
					return;
				}

				//will store retrieved value as promise
				var getObj;

				//reset what existed before
				scope.sharable = {}; //can be reset because re-assign in compilation
				retrievedModel = null;

				//if its a string, it's an ID so get it, otherwise just return a wrapper promise for the object
				if (angular.isObject(sharable) && !angular.isEmpty(sharable)) {
					Debugger.log('sharable object passed');
					getObj = $q.when( sharable );
				} else if (angular.isString(sharable)) {
					//if its a string, call clotho.get()
					Debugger.log('passed a string, assuming a valid ID: ' + sharable);
					getObj = Clotho.get(sharable);
				} else {
					Debugger.warn('sharable must be an object or a string, you passed: ', sharable);
					scope.getPartialAndCompile(null);
					return;
				}

				getObj.then(function (result) {
					retrievedModel = result;
					scope.getPartialAndCompile(ClothoSchemas.determineInstanceType(result));
				});
			};

			//gets partial for type of sharable and compiles editor HTML
			//note - assumes scope already has sharable bound for compilation
			scope.getPartialAndCompile = function (type) {
				scope.showJsonEditor = false;

				if (angular.isEmpty(type)) {
					element.html('<p class="text-center">Please select an Object</p>');
					return;
				}

				//todo - handle instance-specific templates

				$http.get(ClothoSchemas.sharableTypes[type].editor_template_url, {cache: $templateCache})
					.success(function (data) {
						scope.sharable = retrievedModel;
						scope.panelClass = ClothoSchemas.sharableTypes[type].class || 'default';

						var el = $compile(data)(scope);
						element.html(el);
					})
					.error(function (data) {
						Debugger.error('Could not retrieve template for type ' + type);
						element.html('<p class="text-center">Please select an Object</p>');
					});
			};

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
					scope.processInputSharable(null);
				});
			};

			/* watchers */

			//listen for collector_reset and reget & recompile
			Clotho.listen('collector_reset', function Editor_onCollectorReset() {
				scope.processInputSharable(scope.editable);
			}, scope);


			//watch for internal PubSub Changes for given sharable
			var sharableWatcher = angular.noop;
			scope.$watch('sharable', function (newval, oldval) {
				//not empty, and either: no old value, is string and different, or object and id is different
				if (!angular.isEmpty(newval) &&
					 ( !oldval ||
						 ( angular.isString(newval) && newval != oldval ) ||
						 ( newval.id && newval.id != oldval.id)
					 )
				) {
					Debugger.log('sharable is different', newval, oldval);

					//cancel the last watcher, and set up a new one
					sharableWatcher();
					sharableWatcher = Clotho.watch(newval.id, function (newObj) {
						scope.sharable = newObj;
					}, scope);

					//re-compile the actual editor
					scope.processInputSharable(newval);
				}
			});

			/* init */

			//init with input sharable
			scope.processInputSharable(scope.sharable);
		}
	}
});
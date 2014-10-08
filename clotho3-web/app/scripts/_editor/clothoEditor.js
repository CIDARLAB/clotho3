'use strict';

/**
 * @name clotho-editor
 *
 * @param sharable {object} The sharable to be edited. Preferred way of passing in data.
 * @param sharableId {string} ID of sharable to edit. will be retrieved via Clotho.get(), and overwrite sharable passed (if id is passed)
 * @param editMode {*} pass false, null, undefined to start in view mode. Pass anything that passes !!var to start in edit mode
 *
 * @description
 * Will watch for changes to sharable.id
 */
angular.module('clotho.editor')
.directive('clothoEditor', function (Clotho, $compile, $parse, $http, $templateCache, $filter, $q, Debug, ClothoSchemas) {

	var Debugger = new Debug('Editor', '#dd44dd');

	var baseFieldTemplate = 'views/_editor/_baseFields.html',
		formActionsTemplate = 'views/_editor/_formActions.html';

	return {
		restrict: 'A',
		require: '^form',
		scope: {
			sharable: '=?',
			sharableId : '=?',
			editMode: '=?',
			formHorizontal : '=?'
		},
		link: function postLink(scope, element, attrs, formCtrl) {

			if (angular.isUndefined(attrs.sharable) && angular.isUndefined(attrs.sharableId)) {
				Debugger.warn('editor not compiling - did not pass either sharable or sharable-id');
				return;
			}

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
			/*
			scope.removeField = function (name) {
				Debugger.warn('(DEPRECATED) - removing fields from sharables is not supported');
				Debugger.log('removing key ' + name);
				delete scope.sharable[name];
			};
			*/

			/* shared dataTypes */

			ClothoSchemas.retrievedSchemas.then(function (schemas) {
				scope.clothoSchemas = schemas;
			});

			/* compilation */

			function createEmptyEditor () {
				element.html('<p class="text-center">Please select an Object</p>');
			}

			function createEditorElement (sharable) {

				console.log(sharable);

				if (angular.isUndefined(sharable) || !angular.isObject(sharable) || angular.isEmpty(sharable)) {
					Debugger.warn('editor must be created from object, was given: ', sharable);
					createEmptyEditor();
					return;
				}

				ClothoSchemas.determineSharableType(sharable).then(function (type) {
					return type;
				}, function (err) {
					//hack - need to update server version to handle scaffolds (i.e. not in DB)
					//if return null, do dirty check ourselves
					return ClothoSchemas.dirtyDetermineType(sharable)
				})
				.then(function (type) {
					var templateUrl;

					console.warn('sharable type is', type);

					// if it's an instance, check for a more specific template
					if (type == 'Instance') {
						//todo - handle instance-specific templates
						templateUrl = ClothoSchemas.sharableTypes['Instance'].editor_template_url
					} else {
						templateUrl = ClothoSchemas.sharableTypes[type].editor_template_url
					}

					//gets partial for type of sharable and compiles editor HTML
					$http.get(templateUrl, {cache: $templateCache})
						.success(function (data) {
							scope.showJsonEditor = false;
							scope.panelClass = ClothoSchemas.sharableTypes[type].class || 'default';

							var el = $compile(data)(scope);
							element.html(el);
						})
						.error(function (data) {
							Debugger.error('Could not retrieve template: ' + templateUrl);
							element.html('<p class="text-center">Please select an Object</p>');
						});
				});
			}

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
				Clotho.get(scope.sharable.id, {mute: true}).then(function (result) {
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
						scope.sharable.id = id;
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
					createEmptyEditor();
				});
			};

			/*** watchers ***

			 When the sharable changes, we may not need to do anything beyond updating the models.
			 When the id changes, set up a new PubSub watch
			 When the schema changes, the editor needs to recompile. It is assumed the schema is selected via a dropdown, or it will be recompiled on each keystroke.

			 When sharableId changes, that means it has been triggered externally, so we get the ID and set scope.sharable.
			 When the collector resets, we don't want to keep a stale model, so reset the editor.

			*/

			//listen for collector_reset and reget & recompile
			Clotho.listen('collector_reset', function Editor_onCollectorReset() {
				createEmptyEditor();
			}, scope);

			//watch for internal PubSub Changes, update when input sharable id changes
			var sharableWatcher = angular.noop;
			function setNewWatch (id) {
				//cancel the last watcher, and set up a new one
				sharableWatcher();
				sharableWatcher = Clotho.watch(id, function (newObj) {
					if (!angular.equals(scope.sharable, newObj)) {
						scope.sharable = newObj;
					}
				}, scope);
			}

			//note - to just listen to changes from parent coming in, you could use attrs.$observe instead of a scope watch

			//watch the input id if attr is present
			scope.$watch('sharableId', function (newval, oldval) {
				Debugger.log('new id input', newval, oldval);

				newval && Clotho.get(newval, {mute: true}).then(function (result) {
					scope.sharable = result;
				});
			});

			//when schema changes, update
			scope.$watch('sharable.schema', function (newval, oldval) {
				Debugger.log('sharable schema has changed', newval, oldval);
				createEditorElement(scope.sharable);
			});

			//watch the sharable itself for changes to the id
			scope.$watch('sharable.id', function (newval, oldval) {
				Debugger.log('sharable id has changed', newval, oldval);
				setNewWatch(newval);
				if (!angular.isEmpty(scope.sharable) && angular.isUndefined(scope.sharable.schema)) {
					createEditorElement(scope.sharable);
				}
			});

		}
	}
});

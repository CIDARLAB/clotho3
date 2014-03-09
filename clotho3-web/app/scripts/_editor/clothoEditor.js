'use strict';

/**
 * @name clotho-editor
 * 
 * @param sharable {string|object} If String, ID of sharable. Otherwise the sharable directly.
 * @param editMode {*} pass false, null, undefined to start in view mode. Pass anything that passes !!var to start in edit mode
 * 
 * 
 */
angular.module('clotho.editor').directive('clothoEditor', function(Clotho, $compile, $parse, $http, $templateCache, $filter, $q) {

    return {
        restrict: 'A',
        require: ['ngModel', '^form'],
        scope: {
          inputSharable: '=ngModel',
	        editMode: '=?'
        },
        controller: function($scope, $element, $attrs) {
            /******
            GENERAL
            ******/

            //add novalidate to the form to bypass browser validation
            $attrs.$set('novalidate', "novalidate");

	          //if no name is set, set so can use for form validation
	          if (!angular.isString($attrs.name)) {
		          $attrs.$set('name', 'sharableEditor')
	          }

	          //determine whether to start in edit mode
	          $scope.editMode = !!$scope.editMode;

	          //for testing, easily log scope
	          $scope.logScope = function() { console.log($scope); };


            /** config **/
            $scope.schema = {}; //todo - remove

	          //scoped variables for including templates
	          $scope._objBaseFields = 'views/_editor/_baseFields.html';
	          $scope._formActions = 'views/_editor/_formActions.html';

            /** model **/

            //process input sharable to see if object, or id string
	          $scope.processInputSharable = function(sharable) {

		          $scope.sharable = {};
		          console.log(sharable);

		          if (angular.isObject(sharable) && !angular.isEmpty(sharable)) {
			          console.log('sharable object passed');
			          $scope.id = sharable.id;
			          $scope.sharable = sharable;
		          }
		          else if (angular.isString(sharable)) {
			          //if its a string, call clotho.get()
			          console.log('passed a string, assuming a valid ID: ' + sharable);
			          $scope.id = sharable;
		          } else {
			          console.log('sharable must be an object or a string, you passed: ', sharable);
			          //todo - better exit if nothing present
			          $scope.id = '';
			          $scope.sharable = {};
		          }


		          //if its a string, it's an ID so get it, otherwise just return a wrapper promise for the object
		          var getObj = (angular.isEmpty($scope.sharable) && $scope.id != '') ? Clotho.get($scope.id) : $q.when($scope.sharable);

		          getObj.then(function(result) {

			          $scope.sharable = result;

			          $scope.type = $scope.determineType(result);
			          $scope.getPartialAndCompile($scope.type, result);
		          });

	          };

	        //init with input sharable
	        $scope.processInputSharable($scope.inputSharable);

            /*********
            COMPILATION
            **********/

	          //todo - deprecate. needs to be generic to handle all schemas so don't need a switch for each schema name
	          //(or rename all templates and remove this method)
            $scope.determineType = function(obj) {
                function endsWith(str, suffix) {
                    return str.indexOf(suffix, str.length - suffix.length) !== -1;
                }

                if (angular.isUndefined(obj) || angular.isEmpty(obj) || angular.isUndefined(obj.schema)) { return 'undefined' }

                var schema = angular.lowercase(obj.schema);
	            
	              console.log(obj, obj.schema, schema);

                if (endsWith(schema, 'schema'))
                    return 'schema';

                if (endsWith(schema, 'function')) {
                    return 'function';
                }

                //default, even if empty
                return 'sharable';
            };

	          // should check for custom template and use if exists. otherwise, do the form generation for the generic
            $scope.getPartialAndCompile = function(type, obj) {
                $scope.editMode = false;
	              $scope.showJsonEditor = false;

	              //todo - just $http the suspected path based on schema, use error clause to default to generic sharable editor

                switch (angular.lowercase(type)) {
                    case 'function' : {

                        $http.get('views/_editor/function.html',  {cache: $templateCache})
                            .success(function(data) {
                                var el = $compile(data)($scope);
                                $element.html(el);
                            });

                        break;
                    }
                    case 'schema' : {

                        $http.get('views/_editor/schema.html',  {cache: $templateCache})
                            .success(function(data) {
                                var el = $compile(data)($scope);
                                $element.html(el);
                            });

                        break;
                    }
		                case 'sharable' : {

			                $http.get('views/_editor/sharable.html', {cache: $templateCache})
				                .then(function(result) {
					                var el = $compile(result.data)($scope);
					                $element.html(el);
				                });

			                break;

		                }
		                default : {
			                console.log('unknown type');
			                $element.html('<p class="text-center">Please select an Object</p>');

		                }
                }
            };

	          //controller API
	          return {
		          removeField : function (name) {
			          console.log('removing key ' + name);
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

                    scope.edit = function() {
                        scope.editMode = true;
                    };

		                scope.toggleJsonEdit = function() {
			                scope.showJsonEditor = !scope.showJsonEditor;
		                };

                    scope.reset = function() {
                        scope.formCtrl.$setPristine();
                        Clotho.get(scope.id).then(function(result) {
                            scope.sharable = result;
                        });
                    };

                    scope.save = function() {
                        Clotho.set(scope.sharable);
                        scope.editMode = false;
                    };

                    scope.discard = function() {
                        scope.reset();
                        scope.editMode = false;
                    };

                    scope.destroy = function()  {
                        Clotho.destroy(scope.sharable.id).then(function() {
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

                    scope.$watch('inputSharable', function(newval, oldval) {
                        scope.processInputSharable(newval);
                    });

                }
            }
        }
    }
});
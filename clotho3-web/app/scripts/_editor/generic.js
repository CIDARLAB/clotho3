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
        require: '^form',
        scope: {
          sharable: '=?',
	        editMode: '='
        },
        controller: function($scope, $element, $attrs) {
            /******
            GENERAL
            ******/

            console.log($scope.sharable); //testing

            //add novalidate to the form to bypass browser validation
            $attrs.$set('novalidate', "novalidate");

	          //determine whether to start in edit mode
	          $scope.editMode = !!$scope.editMode;

	          //for testing, easily log scope
	          $scope.logScope = function() { console.log($scope); };


            /** config **/
            $scope.schema = {}; //todo - remove
            $scope.formDirty = false; // todo - remove


            /** model **/

            //process input sharable to see if object, or id string
	          $scope.processInputSharable = function(sharable) {
		          $scope.sharable = sharable || $scope.sharable;

		          if (angular.isObject($scope.sharable)) {
			          console.log('sharable object passed');
			          $scope.id = $scope.sharable.id;
		          }
		          else if (angular.isString($scope.sharable)) {
			          //if its a string, call clotho.get()
			          $scope.id = $scope.sharable;
		          } else {
			          console.log('sharable must be an object or a string, you passed: ', $scope.sharable);
			          //todo - better exit if nothing present
		          }
	          };

	          $scope.processInputSharable();

            /*********
            COMPILATION
            **********/


            //wrapper function
            $scope.compileEditor = function() {
                var getObj = ($scope.id == 'id') ? $q.when() : Clotho.get($scope.id);

                getObj.then(function(result) {
                    $scope.sharable = result;

                    $scope.type = $scope.determineType(result);
                    $scope.getPartialAndCompile($scope.type, result);
                });
            };


	          //todo - needs to be generic to handle all schemas so don't need a switch for each schema name
	          //(or rename all templates and remove this method)
            $scope.determineType = function(obj) {
                function endsWith(str, suffix) {
                    return str.indexOf(suffix, str.length - suffix.length) !== -1;
                }

                if (angular.isUndefined(obj) || angular.isEmpty(obj)) { return 'undefined' }

                var schema = angular.lowercase(obj.schema);

                if (endsWith(schema, 'schema'))
                    return 'schema';

                if (endsWith(schema, 'function')) {
                    return 'function';
                }

                //default, even if empty
                return 'sharable';
            };

	          //todo - handle types that do not require additional processing but have custom template (e.g. trail)
	          // should check for custom template and use if exists. otherwise, do the form generation for the generic
            $scope.getPartialAndCompile = function(type, obj) {
                $scope.editMode = false;

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
        },
        compile: function compile(tElement, tAttrs, transclude) {

            return {
                pre: function preLink(scope, iElement, iAttrs, controller) {

                    //note - separation at this point into pre is not hugely important as nothing is linked to form controller
                    scope.compileEditor();

                },
                post: function postLink(scope, iElement, iAttrs, controller) {

                    /* config */

                    scope.form = $parse(iAttrs.name)(scope);

                    scope.edit = function() {
                        scope.editMode = true;
                    };

                    scope.reset = function() {
                        scope.form.$setPristine();
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
                            scope.compileEditor();
                        });
                    };

                    /* watchers */

		                //listen for collector_reset and reget & recompile
		                Clotho.listen('collector_reset', function Editor_onCollectorReset() {
			                scope.compileEditor();
		                }, scope);

		                //watch for internal PubSub Changes
		                Clotho.watch(scope.id, scope.sharable, scope);


	                  //todo - can probably join these two watchers
	                  //watch scope id
                    scope.$watch('id', function(newval, oldval) {
                        if (!!newval && newval != oldval) {
                            scope.compileEditor();
                        }
                    });

                    scope.$watch('sharable', function(newval, oldval) {
                        scope.processInputSharable(newval);
                    });



                }
            }
        }
    }
});
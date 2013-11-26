'use strict';

//todo - integrate ngFormController --- nested ngForms so validation works
// then in template: ng-class="{error: myForm.name.$invalid}"

//todo - avoid compile if object doesn't exist

angular.module('clotho.editor').directive('clothoEditor', function(Clotho, $compile, $parse, $http, $templateCache, $filter, $q) {

    return {
        restrict: 'A',
        require: '^form',
        scope: {
            uuid: '=',
            editable: '=?'
        },
        controller: function($scope, $element, $attrs) {
            /******
            GENERAL
            ******/
            console.log($scope.editable);

            /** config **/

            $scope.schema = {};
            $scope.editMode = false;
            $scope.formDirty = false;
            $scope.logScope = function() { console.log($scope); };

            /** model **/

            //if editable is defined, ignore passed id
            if (angular.isDefined($scope.editable) && angular.isObject($scope.editable)) {
                console.log('editable passed');

                $scope.uuid = $scope.editable.id;

            } else {
                //if we're not linking, just pull it from the attrs -- i.e. put in a string
                if (typeof $scope.uuid == 'undefined') {
                    console.log("id [" + $scope.uuid + "] not in controller, pulling");
                    $scope.uuid = $attrs.uuid;
                }

            }

            /** watchers **/

            //todo - don't assume same schema necessarily
            Clotho.listen('collector_reset', function Editor_onCollectorReset() {
                Clotho.get($scope.uuid).then(function(result) {
                    $scope.editable = result;
                });

            }, $scope);

            Clotho.watch2($scope.uuid, $scope, 'editable', $scope);



            /*********
            COMPILATION
            **********/

            /**
             * @description Generates HTML for a form provided a properly formed schema. Implicit parameters are scope, with $scope.schema defined
             * @returns {string}
             */
            function generateDynamicForm () {
                var fulltext = "";

                //todo
                /*
                string -> test vs. textarea
                when use select?
                what to do for objects?
                id not editable
                 */
                var typeMap = {
                    'string' : 'text',
                    'boolean' : 'checkbox'
                };

                angular.forEach($scope.schema.fields, function(field) {


                    var type = field.type || 'text';
                    if (type == '?') field.type == 'text';
                    var required = field.required ? "required='required'" : "";

                    var htmlText_pre = '<div class="control-group">' +
                        '<label class="control-label" for="' + field.name + '">' + field.name + '</label>' +
                        '<div class="controls">';
                    var htmlText_post = '</div>' +
                        '</div>';
                    var inputText;

                    switch (type) {
                        case "textarea": {
                            inputText = '<textarea class="input-large" id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-model="editable.'+field.name+'" ng-disabled="!editMode"></textarea>';
                            break;
                        }
                        case "select": {
                            var optionsText = "";
                            //todo - use ng-options
                            angular.forEach(field.options, function(value, key) {
                                optionsText = optionsText + '<option value="'+value+'">'+ value + '</option>';
                            });

                            inputText = '<select id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-disabled="!editMode" ng-model="editable.'+field.name+'">' + optionsText + '</select>';
                            break;
                        }
                        case "sharable": {
                        }
                        //todo - add filedrop support, and radio. checkbox works.
                        default: {
                            inputText = '<input type="' + type + '" class="input-large" id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-disabled="!editMode" ng-model="editable.'+field.name+'" >';
                            break;
                        }

                    }

                    fulltext += htmlText_pre + inputText + htmlText_post;
                });
                return fulltext;
            }

            //wrapper function
            $scope.compileEditor = function() {
                var getObj = ($scope.uuid == 'id') ? $q.when() : Clotho.get($scope.uuid);

                getObj.then(function(result) {
                    $scope.editable = result;

                    $scope.type = $scope.determineType(result);
                    $scope.getPartialAndCompile($scope.type, result);
                });
            };


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

            $scope.getPartialAndCompile = function(type, obj) {
                $scope.editMode = false;

                switch (angular.lowercase(type)) {
                    case 'sharable' : {

                        $http.get('views/_editor/sharable.html', {cache: $templateCache})
                            .then(function(result) {
                                var el = $compile(result.data)($scope);
                                $element.html(el);
                            })
                            .then(function() {
                                //todo - fallthrough
                                $scope.schemaName = obj.schema;

                                Clotho.get($scope.schemaName).then(function(result) {
                                    $scope.schema = result;
                                    $scope.schema_custom = result.custom;

                                    var insert = $element.find('insert-fields').html(generateDynamicForm($scope));
                                    $compile(insert.contents())($scope);
                                });
                            });

                        break;

                    }
                    case 'function' : {

                        $http.get('views/_editor/function.html',  {cache: $templateCache})
                            .success(function(data) {
                                var el = $compile(data)($scope);
                                $element.html(el);
                            });

                        break;
                    }
                    case 'schema' : {


	                    var getSuperClass = ($scope.editable.superClass) ?
                            Clotho.get($scope.editable.superClass).then(function(result) {$scope.superClassObj = result})
                            : $q.when($scope.superClassObj = null);

                        getSuperClass.then(function() {
                            $http.get('views/_editor/schema.html',  {cache: $templateCache})
                                .success(function(data) {
                                    var el = $compile(data)($scope);
                                    $element.html(el);
                                });
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
                        Clotho.get(scope.uuid).then(function(result) {
                            scope.editable = result;
                        });
                    };

                    scope.save = function() {
                        Clotho.set(scope.editable);
                        scope.editMode = false;
                    };

                    scope.discard = function() {
                        scope.reset();
                        scope.editMode = false;
                    };

                    scope.destroy = function()  {
                        Clotho.destroy(scope.editable.id).then(function() {
                            scope.editMode = false;
                            scope.editable = null;
                            scope.uuid = null;
                            scope.compileEditor();
                        });
                    };

                    /* watchers */

                    scope.$watch('uuid', function(newval, oldval) {
                        if (!!newval && newval != oldval) {
                            scope.compileEditor();
                        }
                    });


                    scope.$watch('editable', function(newval, oldval) {
                        console.log(newval);

                        if (!!newval && !!newval.id) {
                            scope.uuid = newval.id;
                        } else if (!!newval) {
                            console.log('object doesnt have id -- assume has type');
                            scope.getPartialAndCompile(scope.determineType(newval), newval);
                        } else {
	                        //todo
                        }
                    });

                }
            }
        }
    }
});
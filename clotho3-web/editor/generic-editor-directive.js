'use strict';

//todo - integrate ngFormController --- nested ngForms so validation works
// then in template: ng-class="{error: myForm.name.$invalid}"

Application.Editor.directive('clothoEditor', ['Clotho', '$compile', '$parse', '$http', '$templateCache', '$filter', '$q', function(Clotho, $compile, $parse, $http, $templateCache, $filter, $q) {

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


            /******
             SHARABLE
             ******/



            /******
             FUNCTION
             ******/

            $scope.langTypes = [
                {name:'JavaScript', value:'JAVASCRIPT'},
                {name:'Java', value:'JAVA'},
                {name:'Python', value:'PYTHON'},
                {name:'Groovy', value:'GROOVY'}
            ];

            $scope.outputTypes = [
                {name:'Value', value:'VALUE'},
                {name:'Reference', value:'REFERENCE'}
            ];

            $scope.simpleTypes = {
                "object" : true,
                "string" : true,
                "number" : true,
                "boolean" : true,
                "array" : true
            };

            $scope.paramTypes = [
                {name:'object', type: "object", category:'Primitive', javaType : "java.util.HashMap", reference: false},
                {name:'array', type : "array", category:'Primitive', javaType : "java.util.Arrays", reference: false},
                {name:'string', type : "string", category:'Primitive', javaType : "java.lang.String", reference: false},
                {name:'number', type : "number", category:'Primitive', javaType : "java.lang.Long", reference: false},
                {name:'boolean', type : "boolean", category:'Primitive', javaType : "java.lang.Boolean", reference: false}
            ];
            Clotho.query({schema:"Schema"}).then(function(data){
                angular.forEach(data, function(schema){
                    $scope.paramTypes.push(angular.extend(schema, {category:'Schema'}));
                });
            });

            $scope.clothoFunctions = [];
            Clotho.query({schema: "Function"}).then(function(result) {
                $scope.clothoFunctions = result;
            });



            $scope.addArg = function() {
                if (angular.isEmpty($scope.editable.args)) {$scope.editable.args = [];}
                $scope.editable.args.push({"type" : "", "name" : ""});
            };

            $scope.addDep = function() {
                if (angular.isEmpty($scope.editable.dependencies)) {$scope.editable.dependencies = [];}
                $scope.editable.dependencies.push({"id" : "", "name" : ""});
            };

            $scope.addTest = function() {
                if (angular.isEmpty($scope.editable.tests)) {$scope.editable.tests = [];}
                $scope.editable.tests.push({"args" : [], "output" : {"value" : "", "type" : ""}});
            };

            $scope.testResults = {};
            $scope.singleTest = function(index) {

                var data = {};
                data.id = $scope.editable.id;
                if (angular.isEmpty($scope.editable.tests)) {
                    data.args = [];
                }
                else {
                    data.args = $scope.editable.tests[index].args;
                }

                Clotho.run(data.id, data.args).then(function (result){
                    console.log(result);
                    console.log(result == $scope.editable.tests[index].output.value);
                    $scope.testResults[index] = (result == $scope.editable.tests[index].output.value);
                   /* if (result == angular.fromJson($scope.editable.testResult)) {
                        ClientAPI.say({text:"test success!"});
                    } else {
                        ClientAPI.say({text:"test failed!"});
                    }*/
                });
            };

            $scope.runAllTests = function() {
                for (var i = 0; i < $scope.editable.tests.length; i++) {
                    $scope.singleTest(i);
                }
            };

            $scope.resetTests = function() {
                $scope.testResults = {};
            };

            $scope.querySchemaWrapper = function(schemaType) {
                return Clotho.query({schema: schemaType}).then(function (result) {
                    return $filter('limitTo')(result, 10);
                })
            };


            /*********
             SCHEMA
             **********/

            $scope.schemas = [];
            Clotho.query({"schema": "Schema"}).then(function(data) {
                $scope.schemas = data;
            });

            $scope.accessTypes = [
                {name:'Public', value:'PUBLIC'},
                {name:'Private', value:'PRIVATE'},
                {name:'Read Only', value:'READONLY'}
            ];

            $scope.constraintTypes = [
                {name:'RegExp', value:'regex'},
                {name: 'Not Null', value: 'notnull'}
            ];

            $scope.primitiveToJava = {
                "string" : "java.lang.String",
                "number" : "java.lang.Long",
                "boolean" : "java.lang.Boolean",
                "object" : "java.util.HashMap",
                "array" : "java.util.List"
            };

            $scope.findSpacesRegExp = /\s/ig;

            $scope.parseField = function(field) {
                //only passed field.value so model maps onto options properly in html
                if ($scope.simpleTypes[field.type]) {
                    field.javaType = $scope.primitiveToJava[field.type];
                    field.reference = false;
                } else {

                    field.reference = true;
                }
            };

            $scope.newMethod = function() {
                return ""
            };

            $scope.addMethod = function(method) {
                if (angular.isEmpty($scope.editable.methods)) {$scope.editable.methods = [];}
                $scope.editable.methods.push(method);
            };

            $scope.addNewMethod = function() {
                if (angular.isEmpty($scope.newMethodObj)) return;

                $scope.addMethod($scope.newMethodObj);
                $scope.newMethodObj = $scope.newMethod();
            };

            $scope.newField = function() {
                return {
                    name: "",
                    type: "",
                    description: "",
                    example: "",
                    constraints: null,
                    access: "PUBLIC"
                }
            };

            $scope.addField = function() {
                if (angular.isEmpty($scope.editable.fields)) {$scope.editable.fields = [];}
                $scope.editable.fields.push($scope.newField());
            };

            //note - constraints are processed in saveSchema (link function)
            $scope.newConstraint = function() {
                return {
                    type: "",
                    value: ""
                };
            };

            $scope.addConstraint = function(index) {
                if (angular.isEmpty($scope.editable.fields[index].constraints))
                    $scope.editable.fields[index].constraints = [];
                $scope.editable.fields[index].constraints.push($scope.newConstraint())
            };




            /*********
            COMPILATION
            **********/

            /**
             * @description Generates HTML for a form provided a properly formed schema. Implicit parameters are scope, with $scope.schema defined
             * @returns {string}
             */
            function generateDynamicForm () {
                var fulltext = "";

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

                //todo - better fallthrough
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

                        $http.get('/editor/sharable-partial.html', {cache: $templateCache})
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

                        $http.get('/editor/function-partial.html',  {cache: $templateCache})
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
                            $http.get('/editor/schema-partial.html',  {cache: $templateCache})
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

                    //todo - verify works
                    //iAttrs.$set('novalidate', true);

                    //note - separation at this point into pre is not important as nothing is linked to form
                    scope.compileEditor();

                },
                post: function postLink(scope, iElement, iAttrs, controller) {

                    /* config */

                    scope.form = $parse(iAttrs.name)(scope);
                    scope.resetTests();

                    scope.edit = function() {
                        scope.editMode = true;
                    };

                    scope.createSchema = function(schema) {
                        //create and set to edit mode, changing buttons

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

                    scope.saveSchema = function () {
                        //process constraints to object
                        angular.forEach(scope.editable.fields, function(field) {
                            console.log(field);
                            if (field.constraints) {
                                var constraintsArray = field.constraints;
                                field.constraints = {};
                                angular.forEach(constraintsArray, function(constraint) {
                                    field.constraints[constraint.type] = constraint.value;
                                });
                            } else {
                                field.constraints = null;
                            }
                        });
                        
                        console.log(scope.editable);

                        scope.save();
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
                        } else {
                            console.log('object doesnt have id -- assume has type');
                            scope.getPartialAndCompile(scope.determineType(newval), newval);
                        }
                    });

                }
            }
        }
    }
}]);
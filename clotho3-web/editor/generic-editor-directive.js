'use strict';

//todo - integrate ngFormController --- nested ngForms so validation works
// then in template: ng-class="{error: myForm.name.$invalid}"

Application.Editor.directive('clothoEditor', ['Clotho', '$compile', '$parse', '$http', '$templateCache', '$filter', function(Clotho, $compile, $parse, $http, $templateCache, $filter) {

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

            //todo - map editable to ngModel if present and don't require uuid

            //if we're not linking, just pull it from the attrs -- i.e. put in a string
            if (typeof $scope.uuid == 'undefined') {
                console.log("id [" + $scope.uuid + "] not in controller, pulling");
                $scope.uuid = $attrs.uuid;
            }


            $scope.schema = {};
            $scope.editable = $scope.editable || {};

            $scope.editMode = false;
            $scope.formDirty = false;

            //todo - don't assume same schema necessarily
            Clotho.listen('collector_reset', function Editor_onCollectorReset() {
                Clotho.get($scope.uuid).then(function(result) {
                    $scope.editable = result;
                });

            }, $scope);

            Clotho.watch2($scope.uuid, $scope, 'editable', $scope);

<<<<<<< HEAD
=======
            $scope.logScope = function() { console.log($scope); };

>>>>>>> refactoring-overhaul

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

<<<<<<< HEAD
=======
            $scope.outputTypes = [
                {name:'Value', value:'VALUE'},
                {name:'Reference', value:'REFERENCE'}
            ];
>>>>>>> refactoring-overhaul

            $scope.simpleTypes = {
                "Object" : true,
                "String" : true,
                "Integer" : true,
                "Boolean" : true
            };

            $scope.paramTypes = [
                {name:'Object', type:'Type'},
                {name:'String', type:'Type'},
                {name:'Integer', type:'Type'},
                {name:'Boolean', type:'Type'}
            ];
            Clotho.query({schema:"Schema"}).then(function(data){
                angular.forEach(data, function(schema){
                    $scope.paramTypes.push({name:schema.name, type:'Schema'});
                });
            });



            $scope.addArg = function() {
                if (angular.isEmpty($scope.editable.args)) {$scope.editable.args = [];}
<<<<<<< HEAD
                $scope.editable.args.push({"type" : "", "name" : "", "test" : {"uuid" : ""}});
=======
                $scope.editable.args.push({"type" : "", "name" : ""});
>>>>>>> refactoring-overhaul
            };

            $scope.addDep = function() {
                if (angular.isEmpty($scope.editable.dependencies)) {$scope.editable.dependencies = [];}
                $scope.editable.dependencies.push({"id" : "", "name" : ""});
            };

<<<<<<< HEAD
            $scope.testFunction = function() {

                var data = {};
                data.id = $scope.editable.id;
                if (angular.isEmpty($scope.editable.args)) {$scope.editable.args = [];}
                data.args = $scope.editable.args.map(function (param){
                    return param.test.uuid ? param.test.uuid : param.test.value;
                });

                Clotho.run(data.id, data.args).then(function (result){
                   Clotho.say({text: result});
=======
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
>>>>>>> refactoring-overhaul
                   /* if (result == angular.fromJson($scope.editable.testResult)) {
                        ClientAPI.say({text:"test success!"});
                    } else {
                        ClientAPI.say({text:"test failed!"});
                    }*/
                });
            };

<<<<<<< HEAD
=======
            $scope.runAllTests = function() {
                for (var i = 0; i < $scope.editable.tests.length; i++) {
                    $scope.singleTest(i);
                }
            };

>>>>>>> refactoring-overhaul
            $scope.queryWrapper = function(schemaType) {
                return Clotho.query({schema: schemaType}).then(function (result) {
                    return $filter('limitTo')(result, 10);
                })
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

<<<<<<< HEAD
                angular.forEach($scope.schema, function(field) {

                    var type = field.type || 'text';
                    var required = field.required ? "required='required'" : "";

                    var htmlText_pre = '<div class="control-group">' +
                        '<label class="control-label" for="' + field.name + '">' + field.readable + '</label>' +
=======
                angular.forEach($scope.schema.fields, function(field) {

                    var type = field.type || 'text';
                    if (type == '?') field.type == 'text';
                    var required = field.required ? "required='required'" : "";

                    var htmlText_pre = '<div class="control-group">' +
                        '<label class="control-label" for="' + field.name + '">' + field.name + '</label>' +
>>>>>>> refactoring-overhaul
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
<<<<<<< HEAD
=======
                        case "sharable": {
                        }
>>>>>>> refactoring-overhaul
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

            $scope.compileEditor = function() {

                Clotho.get($scope.uuid).then(function(result) {
                    $scope.editMode = false;
                    $scope.editable = result;
<<<<<<< HEAD
                    $scope.type = result.type;

                    if (angular.lowercase(result.schema).indexOf('function', ((result.schema).length - 8)) !== -1) {
                        $scope.type = 'function';
                    }

=======

                    //todo - rewrite
                    //$scope.type = result.type;
                    var suffix = 'function';
                    var str = angular.lowercase(result.schema);
                    var endswith = str.indexOf(suffix, str.length - suffix.length) !== -1
                    if (endswith){
                        $scope.type = 'function';
                    } else {
                        $scope.type = 'sharable';
                    }


                    if (angular.lowercase(result.schema).indexOf('function', ((result.schema).length - 8)) !== -1)
                        $scope.type = 'function';

>>>>>>> refactoring-overhaul
                    switch (angular.lowercase($scope.type)) {
                        case 'sharable' : {

                            $scope.schemaName = result.schema;

                            $http.get('/editor/sharable-partial.html', {cache: $templateCache})
                                .then(function(result) {
                                    var el = $compile(result.data)($scope);
                                    $element.html(el);
                                })
                                .then(function() {
                                    Clotho.get($scope.schemaName).then(function(result) {
<<<<<<< HEAD
                                        $scope.schema = result.schema;
=======
                                        $scope.schema = result;
>>>>>>> refactoring-overhaul
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
<<<<<<< HEAD
                                    $element.html(data);
                                    $compile($element.contents())($scope);
=======
                                    var el = $compile(data)($scope);
                                    $element.html(el);
>>>>>>> refactoring-overhaul
                                });

                            break;
                        }
                        case 'schema' : {


                            break;
                        }
                        default : {
                            console.log('unknown type');

                        }
                    }
                });
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

                    scope.form = $parse(iAttrs.name)(scope);

                    scope.$watch('uuid', function(newval, oldval) {
                        if (!!newval && newval != oldval) {
                            scope.compileEditor();
                        }
                    });

                    //e.g. scope.formConst.$setPristine()
                    scope.formConst = $parse(iAttrs.name)(scope);
                    //switch to 'edit' mode
                    scope.edit = function() {
                        scope.editMode = true;
                    };

                    //discard edits
                    scope.reset = function() {
                        scope.formConst.$setPristine();
                        Clotho.get(scope.uuid).then(function(result) {
                            scope.editable = result;
                        });
                    };

                    //save edits, switch to 'view'
                    scope.save = function() {
                        Clotho.set(scope.editable);
                        scope.editMode = false;
                    };

                    //discard edits, switch to 'view'
                    scope.discard = function() {
                        scope.reset();
                        scope.editMode = false;
                    };

                    scope.destroy = function()  {
                        Clotho.destroy(scope.editable.id).then(function() {
                            console.log('bam');
                            scope.editMode = false;
                            scope.editable = undefined;
                            scope.id = undefined;
                        });

                    };

                }
            }
        }
    }
}]);
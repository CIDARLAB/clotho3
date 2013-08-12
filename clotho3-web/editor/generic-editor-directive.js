'use strict';

//todo - integrate ngFormController --- nested ngForms so validation works
// then in template: ng-class="{error: myForm.name.$invalid}"

Application.Editor.directive('clothoEditor', ['Clotho', '$compile', '$parse', '$http', '$templateCache', '$filter', function(Clotho, $compile, $parse, $http, $templateCache, $filter) {

    return {
        restrict: 'A',
        require: '^form',
        scope: {
            uuid: '='
            //fixme - make not required, doesn't work otherwise
            //editable: '='
        },
        controller: function($scope, $element, $attrs) {
            /******
            GENERAL
            ******/

            //if we're not linking, just pull it from the attrs -- i.e. put in a string
            if (typeof $scope.uuid == 'undefined') {
                console.log("id [" + $scope.uuid + "] not in controller, pulling");
                $scope.uuid = $attrs.uuid;
            }

            $scope.schema = {};
            $scope.editable = $scope.editable || {};

            $scope.editMode = false;
            $scope.formDirty = false;

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
                $scope.editable.args.push({"type" : "", "name" : "", "test" : {"uuid" : ""}});
            };

            $scope.addDep = function() {
                if (angular.isEmpty($scope.editable.dependencies)) {$scope.editable.dependencies = [];}
                $scope.editable.dependencies.push({"id" : "", "name" : ""});
            };

            $scope.testFunction = function() {

                var data = {};
                data.id = $scope.editable.id;
                if (angular.isEmpty($scope.editable.args)) {$scope.editable.args = [];}
                data.args = $scope.editable.args.map(function (param){
                    return param.test.uuid ? param.test.uuid : param.test.value;
                });

                Clotho.run(data.id, data.args).then(function (result){
                   Clotho.say({text: result});
                   /* if (result == angular.fromJson($scope.editable.testResult)) {
                        ClientAPI.say({text:"test success!"});
                    } else {
                        ClientAPI.say({text:"test failed!"});
                    }*/
                });
            };

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

                angular.forEach($scope.schema, function(field) {

                    var type = field.type || 'text';
                    var required = field.required ? "required='required'" : "";

                    var htmlText_pre = '<div class="control-group">' +
                        '<label class="control-label" for="' + field.name + '">' + field.readable + '</label>' +
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
                    $scope.type = result.type;

                    if (angular.lowercase(result.schema).indexOf('function', ((result.schema).length - 8)) !== -1)
                        $scope.type = 'function';

                    switch (angular.lowercase($scope.type)) {
                        case 'sharable' : {

                            $scope.schemaName = result.schema_id;

                            $http.get('/editor/sharable-partial.html', {cache: $templateCache})
                                .then(function(result) {
                                    var el = $compile(result.data)($scope);
                                    $element.append(el);
                                })
                                .then(function() {
                                    Clotho.get($scope.schemaName).then(function(result) {
                                        $scope.schema = result.schema;
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
                                    $element.html(data);
                                    $compile($element.contents())($scope);
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

                    //verify works
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

                }
            }
        }
    }
}]);
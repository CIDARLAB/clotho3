'use strict';

//todo - sharable bindable??
//decide if want to make sharable bindable also / instead of uuid (so can bind to models outside of directive)
// would need watch statement on sharable

Application.Editor.directive('sharableEditor', ['Clotho', '$compile', '$parse', function(Clotho, $compile, $parse) {

    return {
        restrict: 'A',
        replace: false,
        template: '<legend>{{schemaName}} View</legend>' +
            '<insert-fields></insert-fields>' +
            '<div class="form-actions" ng-hide="editMode">' +
            '   <button class="btn" type="button" ng-click="edit()">Edit</button>' +
            '   <button type="button" class="btn btn-inverse" ng-click="logScope()">Log</button>' +
            '</div>' +
            '<div class="form-actions" ng-show="editMode">' +
            '   <button type="submit" class="btn btn-primary" ng-click="save()" ng-disabled="sharableEditor.$invalid">Save</button>' +
            '   <button type="submit" class="btn btn-danger" ng-click="discard()">Discard</button>' +
            '   <button type="button" class="btn btn-warning" ng-click="reset()">Reset</button>' +
            '   <button type="button" class="btn btn-inverse" ng-click="logScope()">Log</button>' +
            '</div>',
        scope: {
            //sharable: '=',
            uuid: '='
        },
        controller: function($scope, $element, $attrs) {
            //if we're not linking, just pull it from the attrs
            if (typeof $scope.uuid == 'undefined') {
                console.log("uuid not in controller, pulling");
                $scope.uuid = $attrs.uuid;
            }

            $scope.schema = {};
            $scope.sharable = $scope.sharable || {};

            $scope.editMode = false;
            $scope.formDirty = false;

            // Listeners

            Clotho.listen('collector_reset', function SharableEditor_onCollectorReset() {

                //note - this uses promises in the default way, angular knows how to deal with it
                //$scope.sharable = Clotho.get($scope.uuid, true);

                //this is how people would ideally access things (normal promise pattern)
                Clotho.get($scope.uuid).then(function(result) {
                    $scope.sharable = result;
                });

                // testing - sync Clotho.get()
                //console.log(Clotho.get($scope.uuid, true));
                
            }, $scope);

            /*
            //note - alternate version - see also watch2 below
            Clotho.watch($scope.uuid, function (data) {
                $scope.sharable = data;
            }, $scope);
             */

            Clotho.watch2($scope.uuid, $scope, 'sharable', $scope);

            //future - use return {} syntax? i.e. for inheritable directive controllers
        },
        compile: function compile(tElement, tAttrs, transclude) {

            return {
                pre: function preLink(scope, iElement, iAttrs, controller) {

                    //todo - check out ngView and see how it compiles / links async content.
                    //todo - interface with ngForm
                    scope.generate_fields = function() {
                        var insert = iElement.find('insert-fields').html('');
                        var fulltext = "";

                        angular.forEach(scope.schema, function(field) {

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
                                    inputText = '<textarea class="input-large" id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-model="sharable.'+field.name+'" ng-disabled="!editMode"></textarea>';
                                    break;
                                }
                                case "select": {
                                    var optionsText = "";
                                    angular.forEach(field.options, function(value, key) {
                                        optionsText = optionsText + '<option value="'+value+'">'+ value + '</option>';
                                    });

                                    inputText = '<select id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-disabled="!editMode" ng-model="sharable.'+field.name+'">' + optionsText + '</select>';
                                    break;
                                }
                                default: {
                                    inputText = '<input type="' + type + '" class="input-large" id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-disabled="!editMode" ng-model="sharable.'+field.name+'" >';
                                    break;
                                }

                            }

                            fulltext += htmlText_pre + inputText + htmlText_post;
                        });
                        insert.html(fulltext);
                        $compile(insert.contents())(scope);
                    };


                    //get the sharable (which says which schema it needs)
                    //note variables not compiled yet in 'pre' (e.g. if use scope: {var : '@'} and should go through $parse)
                    Clotho.get(scope.uuid).then(function(result) {
                        scope.sharable = result;
                        scope.schemaName = result.$clotho.schema;

                        //get the schema
                        Clotho.get(scope.schemaName).then(function(result) {
                            scope.schema = result.schema;
                            scope.schema_custom = result.custom;
                            scope.generate_fields();
                        });
                    });
                },
                post: function postLink(scope, iElement, iAttrs, controller) {

                    //e.g. scope.formConst.$setPristine()
                    scope.formConst = $parse(iAttrs.name)(scope);

                    //todo - move to native angular code to handle this
                    //fixme - using directive to generate fields currently prevents ngForm from working
                    //avoid complex statements in watch statements (e.g. do we need the else if?)
                    //this statement only works for when form elements do not have new scope
                    scope.$watch(iAttrs.name + '.$dirty', function (newValue, oldValue) {
                        if (newValue != oldValue && newValue === true) {
                            scope.formDirty = true;
                        } else if (newValue != oldValue && newValue === false) {
                            scope.formDirty = false;
                        } else {}
                    });

                    /*
                     functions
                     note - not in directive controller because reinstantiated each change
                    */

                    //switch to 'edit' mode
                    scope.edit = function() {
                        scope.editMode = true;
                    };

                    //discard edits
                    scope.reset = function() {
                        scope.formConst.$setPristine();
                        Clotho.get(scope.uuid).then(function(result) {
                            scope.sharable = result;
                        });
                    };

                    //save edits, switch to 'view'
                    scope.save = function() {
                        Clotho.set(scope.uuid, scope.sharable);
                        scope.editMode = false;
                    };

                    //discard edits, switch to 'view'
                    scope.discard = function() {
                        scope.reset();
                        scope.editMode = false;
                    };

                    //testing
                    scope.logScope = function() {
                        console.log(scope);
                    };
                }
            }
        }
    }
}]);

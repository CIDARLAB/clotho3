'use strict';

/* TODO - pass in UUID, get schema and fields from within directive
    - rearchitect model binding (i.e. passing in of sharable)
 */

Application.Directives.directive('editorGenerator', ['$compile', function($compile) {
    return {
        restrict: 'E',
        replace: true,
        //future - choose which are necessary,
        scope: {
            sharable: '=',
            schema: '@',
            schemaName: '@',
            uuid: '=',
            editMode: '='
        },
        controller: function($scope, $element, $attrs, Collector) {
            $scope.logScope = function() {
                console.log($scope);
            }
        },
        link: function(scope, element, attrs, controller) {
            /* testing
            console.log("editor_generator");
            console.log(scope.schema);
            console.log(scope);
            */

            //todo - add ng-multiple for selects (http://docs.angularjs.org/api/ng.directive:ngMultiple)

            function generate_fields() {
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

                    var htmlText = $compile(htmlText_pre + inputText + htmlText_post)(scope);
                    element.append(htmlText);

                });
            }

            /* future -- check back for this being fixed
               todo - try to get this resolved via the resolve in router configuration
            - values are not interpolated before link function runs
            - i.e. promise not passed through, so schema is undefined when directive compiles
            - see: https://github.com/angular/angular.js/pull/1555
            - current solutions
                - use attrs.$observe('variable', function(changedValue){ //do stuff });
                    - see e.g. http://plnkr.co/edit/ICslXxkIuyB0Ok0olKpS?p=preview
                    - requires isolate scope is @ and {{}} in template
                - require promise resolved before directive called
                    - e.g. https://groups.google.com/forum/?fromgroups=#!searchin/angular/directive$20promise$20undefined/angular/VEtVYUptpPs/5DnnBwCZ0ZwJ
             */
            attrs.$observe('schema', function(changedValue){
                if (typeof changedValue !== "undefined" && changedValue != "")
                    scope.schema = JSON.parse(changedValue);
                generate_fields();
            });

            //testing
            var logScopeHTML = $compile('<a class="btn btn-info" ng-click="logScope()">LogScope_EditorGenerator</a>')(scope);
            element.append(logScopeHTML);

        }
    }
}]);
'use strict';

Application.Directives.directive('formInput', function($compile) {
    //FUTURE: to make more modular - set up so each has its own scope? could pull from some sort of schema resource service
    return {
        restrict: 'E',
        replace: true,
        /*
         NOTES
         - basic prototypical scope inheritance shouldn't cause problems (e.g. transclude = true, scope = true)
         - setting up isolate scope does cause problems, presumably because of the ng-repeat
         - ... will be apparent when have more components using this directive
         - may make sense to have a directive that ITSELF loops to create elements ...
         and to not have the loop outside the directive

         //transclude: true,

         //need to explicitly declare these e.g. <form-input field="{{field}}" ...>
         scope: {
         field: '=',
         model: '='
         },
         */
        controller: function($scope, $element, $attrs, Collector) {
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, iElement, iAttrs, controller) {

                },
                post: function postLink(scope, iElement, iAttrs, controller) {

                    //console.log(iAttrs);
                    //console.log(iElement);
                    console.log(scope);

                    var field = scope.field;
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
                    iElement.append(htmlText);
                }
            }
        }
    }
});